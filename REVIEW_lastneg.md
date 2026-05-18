# Code Review: `add-lastneg-likelihood` branch

**Date:** 2026-05-17
**Branch:** `add-lastneg-likelihood` vs `master`
**Scope:** All R and C++ changes adding last-negative test date support per host.

## Context

The branch threads `last.negative[i]` (date of host `i`'s last confirmed-negative test) into the model, with the intent that the inferred infection time `t_inf[i]` must satisfy `t_inf[i] > t_neg[i]`. The implementation adds:

- A new MCMC proposal `propose_tinf_lastneg` using rejection sampling.
- A Hastings-ratio correction `lastneg_logratio_correction` for the truncated proposal.
- NaN/Inf guards in MCMC acceptance steps.
- Plumbing of `last.negative` through `phybreakdata()`, `build_pbe`, and several `update_host_*` paths.
- A C++ negative-clamp in `likseqenv.cpp` and a coalescent-time shift heuristic in R.

Tests are assumed perfect (sensitivity 1, specificity 1) — this is not documented anywhere.

---

## Critical issues (likely invalidate the inference)

### C1. The constraint is not in the target distribution
**Files:** [R/logLik_phybreak.R:101-103](R/logLik_phybreak.R#L101-L103), [R/infect_distribution.R](R/infect_distribution.R)

`lik_sampletimes` is the unmodified Gamma density on `sampletime − inftime`. Nowhere in the joint log-posterior is there either:
- a hard indicator `1{t_inf > t_neg}`, or
- the truncated-density normalizer `−log(1 − F(t_neg − t_inf | sample params))`.

The branch only modifies proposals, so MCMC converges to neither the constrained nor the unconstrained target — detailed balance is broken. This also biases `sample.mean`/`sample.shape` toward shorter intervals.

**Fix:** add
```r
- sum(log1p(- pgamma(lastneg - inftimes, sample.shape, sample.shape / sample.mean)))
```
over hosts with a last-negative date, in `lik_sampletimes` (or a new `lik_lastneg`). Once present, the proposal and target agree.

### C2. Path I in `update_host_keepphylo` uses the wrong `tinf.prop`
**File:** [R/mcmc-updatehost-paths.R:127](R/mcmc-updatehost-paths.R#L127)

`tinf.prop <- v$inftimes[hostiorID]` overwrites the value drawn by `propose_tinf_lastneg` before `update_pathI` is dispatched. The downstream `lastneg_logratio_correction` at line 979 then receives the infector's existing infection time, not the proposal. Detailed balance is broken whenever `d$last.negative[hostID]` is set.

**Fix:** carry the original sampled `tinf.prop` separately and pass it to the correction, or skip the correction for path I.

### C3. Paths C and F have explicit "biased otherwise" TODOs
**File:** [R/mcmc-updatehost-paths.R:392-398](R/mcmc-updatehost-paths.R#L392-L398), [:549-554](R/mcmc-updatehost-paths.R#L549-L554)

The MH correction for index-swap proposals is not derived. The author's own comment says "Acceptable when last.negative is unset for hostID; biased otherwise."

**Fix (simplest):** early `return()` from `update_pathC` / `update_pathF` when either host involved in the swap has a finite `last.negative`.

### C4. Within-host rewire moves do not check the constraint
**Files:** [R/mcmc-rewire_paths_complete_edgewise.R:39](R/mcmc-rewire_paths_complete_edgewise.R#L39), [R/mcmc-rewire_paths_wide_edgewise.R:38](R/mcmc-rewire_paths_wide_edgewise.R#L38)

These set `pbe1$v$inftimes[ID] <- tinf` without enforcing `tinf > last.negative[ID]`.

**Fix:**
```r
if (!is.na(d$last.negative[ID]) && tinf < d$last.negative[ID]) {
  pbe1$logLiktoporatio <- -Inf
  return()
}
```

### C5. Likelihood components fabricated on NaN/Inf
**File:** [R/mcmc-environment-functions.R:154-167](R/mcmc-environment-functions.R#L154-L167)

When `logLikseq` is invalid, the code sets `logLiksam`, `logLikgen`, `logLikcoal` all to `-Inf` rather than computing them. Downstream acceptance ratios using these placeholders will reject everything and freeze the chain at the last accepted state.

**Fix:** detect once, log it with state info, halt the chain (or revert to last good state). Do not fabricate components.

### C6. Out-of-bounds read in C++ debug branch
**File:** [src/likseqenv.cpp:123](src/likseqenv.cpp#L123)

`likarray[(rootnode * nSNPs + j) * 4 + k]` is dereferenced after the `for` loop exits with `k == 4`, reading into the next SNP's slot.

**Fix:** drop the `+ k`, print the 4-vector explicitly, and switch from `std::cout` to `Rcpp::Rcout`.

### C7. `test_lik_sampletimes.R` signature mismatch
**File:** [tests/testthat/test_lik_sampletimes.R:14,20](tests/testthat/test_lik_sampletimes.R#L14)

Calls `lik_sampletimes(...)` with a 6th `lastneg` argument, but [R/logLik_phybreak.R:101](R/logLik_phybreak.R#L101) takes 5 args. The test errors immediately with "unused argument (lastneg)" and exercises nothing.

**Fix:** either remove the extra argument, or actually thread `lastneg` through `lik_sampletimes` (preferred, since C1 needs this anyway).

### C8. `update_mS` Gibbs draw assumes the untruncated sampling likelihood
**File:** [R/mcmc-updateparameters.R:37-69](R/mcmc-updateparameters.R#L37-L69) (`update_mS`), [:71-98](R/mcmc-updateparameters.R#L71-L98) (`update_mG`)

`update_mS` draws `sample.mean` from a Gamma full conditional that assumes the *untruncated* sampling-time likelihood, and only includes `logLikcoal` in the MH ratio (relying on Gibbs cancellation for `logLiksam`). Once the truncation normalizer `Π_i s(t_neg_i − t_inf_i; sample.shape, sample.mean)` is in the target, the Gibbs distribution is no longer the true conditional, and the MH correction omits the change in `Σ log s(·)`. `update_mG` does the same and then accepts unconditionally. The chain will drift `sample.mean` and `gen.mean` to biased values whenever any host has a last-negative date.

**Fix:** rewrite `update_mS` as a proper MH (drop the Gibbs draw and include `logLiksam_new − logLiksam_old + Σ log s_new − Σ log s_old` in `logaccprob`). `update_mG` must be MH at all.

### C9. `update_move_sampleedges()` silently skipped under last-negative
**File:** [R/mcmc-updatehost-paths.R:693-707](R/mcmc-updatehost-paths.R#L693-L707) (`update_move`)

`update_move_sampleedges()` was relocated *inside* the `if(logLiktoporatio > -Inf && logproposalratio > -Inf)` block AND inside the new `!is.na(logacceptanceprob)` guard. Previously it ran unconditionally after the proposal block. Within-host phylo updates are now silently skipped whenever the host-level transmission proposal is infeasible — which is very common when `last.negative` constrains `tinf`. Mixing of within-host moves collapses in last-negative-rich datasets.

**Fix:** hoist `if(which_protocol == "edgewise") update_move_sampleedges()` back out of the nested guards.

### C10. Path J also missing `lastneg_logratio_correction` for `hostiorID`
**File:** [R/mcmc-updatehost-paths.R:1085](R/mcmc-updatehost-paths.R#L1085)

Similar to C2: in path J, `hostiorID` receives a freshly Gamma-sampled `tinf2.prop` via `propose_tinf_lastneg`, but only `hostID`'s correction term is added to `logproposalratio`. If `hostiorID` has `last.negative` set, the corresponding correction is missing.

**Fix:** add `+ lastneg_logratio_correction(hostiorID, pbe0$v$inftimes[hostiorID], tinf2.prop, p, d)` to `logproposalratio` in path J. Audit path I similarly for the second host.

### C11. Path I `lastneg_logratio_correction` missing normalizer
**File:** [R/mcmc-updatehost-paths.R:961-979](R/mcmc-updatehost-paths.R#L961-L979)

The rejection-sampler proposal density is `q(t) ∝ Gamma(nodetime − t; shape=2/3·sample.shape) · [1 − F(t_neg − t; shape=sample.shape)]`, normalized by `Z(nodetime)`. For paths A/B/D/E/G/H/J where the forward and reverse proposals share the same `nodetime`, `Z` cancels. In path I, the reverse proposal uses `pbe0$v$nodetimes[hostiorID]` (a different `nodetime`), so `Z` does NOT cancel and the current correction is missing `log Z(nodetime_hostID) − log Z(nodetime_hostiorID)`.

**Fix:** compute `Z` numerically (1D integral, e.g., `integrate()` or a closed-form via incomplete gammas) per proposal in path I (and audit path L analogously), or restrict `last.negative` to never apply to hosts that can hit path I/L.

### C12. C++ negative-clamp `else` branch does not actually clamp
**File:** [src/likseqenv.cpp:85-93](src/likseqenv.cpp#L85-L93)

`if (val > -1e-3) val = 0.0;` else { prints via `Rcpp::Rcout` but `val` is **not** clamped }. The subsequent `likarray[...] *= val` writes a negative product. `log(SNPsums[j])` then produces `NaN`, which trips the new pbe1$logLikseq NaN-handler (C5) and corrupts the chain.

**Fix:** in the `else` branch, also assign `val = 0.0` or return `-Inf` cleanly. The print-only branch silently injects NaN downstream.

### C13. `lik_func` loop runs on a state already known to be corrupt
**File:** [R/mcmc-environment-functions.R:197-201](R/mcmc-environment-functions.R#L197-L201)

When `logLikseq` is invalid, the branch sets `logLiksam/logLikgen/logLikcoal <- -Inf` (see C5) but then unconditionally runs the user-supplied `lik_func` block at line 197. Custom likelihood modules then evaluate on a known-corrupt state, may themselves return `NaN`, get copied into `pbe0`, and corrupt downstream `prepare_pbe` checks.

**Fix:** move the `lik_func` loop inside the `else` branch (skip when `logLikseq` is invalid).

### C14. Kinetic boundary leakage across `t_neg`
**Files:** [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R), [R/mcmc-rewire_paths_*](R/)

Since `lik_sampletimes` does NOT contain the truncation factor (C1), the target is currently smooth across `t_neg`. The truncation is enforced only via the rejection sampler in `propose_tinf_lastneg`. Other moves (within-host paths K, path I/L when normalizer is wrong, any path bypassing the correction) can propose `tinf_old < t_neg ≤ tinf_new` and accept happily. Combined with the missing normalizer in the target, the sampler conditions the proposal on `t_inf > t_neg` but does not penalize positions where it isn't — an asymmetric kernel that violates detailed balance. The chain can leak across the boundary and never bounce back.

**Fix:** automatically resolved once C1 is in place — but every path's constraint check must then be audited (C4, C8 above).

---

## Important issues (probable bugs)

### I1. MH guards silently mask `+Inf` and `NaN`
**Files:** [R/mcmc-updateparameters.R](R/mcmc-updateparameters.R), [R/modules.R](R/modules.R), [R/swap_heats.R](R/swap_heats.R)

`if (!is.na(logaccprob) && !is.nan(logaccprob) && !is.infinite(logaccprob))` silently no-ops on `+Inf` (which should always accept) and `NaN` (which signals a bug). `-Inf` is already handled correctly by `runif(1) < exp(logaccprob)`.

**Fix:** accept on `+Inf`, log a warning on `NaN`. Do not silently no-op.

### I2. Rejection sampler silent bailout
**File:** [R/mcmc-updatehost-paths.R:34-41](R/mcmc-updatehost-paths.R#L34-L41)

`propose_tinf_lastneg` returns `NA` after `tinf.lastneg.maxtries = 1000`; the caller silently `return()`s. For hosts where `nodetime − lastneg` is small vs `sample.mean`, mixing collapses with no diagnostic.

**Fix:** replace rejection sampling with one-shot inverse-CDF draw from the truncated Gamma:
```r
F0 <- pgamma(nodetime - lastneg, shape, rate)
qgamma(F0 + runif(1) * (1 - F0), shape, rate)
```

### I3. `admission.times` plumbing is broken
**File:** [R/phybreak.R:168](R/phybreak.R#L168)

`# dataslot$admission.times <- dataset$admission.times` is commented out. The check at line 256 is therefore always FALSE, and the admission constraint in [mcmc-updatehost-paths.R:194-195](R/mcmc-updatehost-paths.R#L194-L195) never fires.

**Fix:** uncomment line 168.

### I4. No initialization rejection for last-negative
**File:** [R/sample_phybreak.R](R/sample_phybreak.R)

Once C1 is fixed, initial states with `inftimes[i] < last.negative[i]` will have `-Inf` likelihood.

**Fix:** resample violating values from the truncated Gamma at init, then refresh downstream node times.

### I5. Dead date-to-numeric conversion in `propose_pbe`
**File:** [R/mcmc-environment-functions.R:323-326](R/mcmc-environment-functions.R#L323-L326)

`build_pbe` already converted `d$last.negative` to numeric at lines 53-56, so the `inherits(d$last.negative, "Date")` check is always FALSE. Even if it fired, the local reassignment is not propagated.

**Fix:** remove the block; rely on the single conversion in `build_pbe`.

### I6. Wrong-class comparisons in `phybreakdata`
**File:** [R/phybreakdata.R:269,276](R/phybreakdata.R#L269)

`class(last.negative) != class(sample.times)` and same for `admission.times` compare character vectors; breaks on subclasses.

**Fix:** use `inherits(last.negative, class(sample.times))`.

### I7. Positional-arg breakage in `phybreakdata()`
**File:** [R/phybreakdata.R](R/phybreakdata.R)

New `admission.times`/`last.negative` were inserted between existing args, silently breaking positional callers.

**Fix:** move new args to the end of the signature.

### I8. `is.nan(logLikseq)` checked before assignment
**File:** [R/mcmc-environment-functions.R:181-195](R/mcmc-environment-functions.R#L181-L195)

`logLikseq` is referenced before any visible assignment in the function. Either `.likseqenv` writes into the calling env (fragile), or this throws `object 'logLikseq' not found`.

**Fix:** make the binding explicit; `logLikseq <- .likseqenv(...)`.

### I9. `prepare_pbe` does `stop()` with no message
**File:** [R/mcmc-environment-functions.R:250-253](R/mcmc-environment-functions.R#L250-L253)

Hard-aborts the chain instead of letting the proposal be rejected.

**Fix:** set `pbe1$logLikseq <- -Inf`, return, and let the MH step reject. Or `stop("informative message")` if truly unrecoverable.

### I10. C++ negative-clamp masks upstream issue
**File:** [src/likseqenv.cpp:82-94](src/likseqenv.cpp#L82-L94)

`if (val > -1e-3) val = 0.0` silently changes the Felsenstein recursion. Negative partial likelihoods should not arise unless edge lengths are non-positive — almost certainly tied to the `eps = 1e-6` coalescent-time shift heuristic at [R/mcmc-environment-functions.R:115-135](R/mcmc-environment-functions.R#L115-L135).

**Fix:** investigate why coalescent times are not strictly older than tip times at build, instead of patching both the tree and the C++.

### I11. `(ns+1):N` blows up when `ns == 1`
**File:** [R/mcmc-environment-functions.R:111](R/mcmc-environment-functions.R#L111)

`(2):1` reverse-iterates.

**Fix:** guard with `if (N >= ns + 1)` or `seq.int(ns+1, N, length.out = max(0, N-ns))`.

### I12. Detailed-balance audit for moves touching `t_inf`
**Files:** all `update_path*` and `rewire_*` in [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R), [R/mcmc-rewire_paths_*](R/)

After C1 is fixed (likelihood includes truncation normalizer), every move that touches any `t_inf[i]` must (a) check the constraint and (b) include the change in the normalizer in its acceptance ratio. Audit all paths systematically.

### I13. Parallel-tempering swap correctness depends on where the normalizer lives
**File:** [R/swap_heats.R:38-43](R/swap_heats.R#L38-L43), [R/sample_phybreak.R:243-245](R/sample_phybreak.R#L243-L245)

PT `likelihood` vector is built by summing all `pbe$logLik*` slots. Currently both cold and tempered chains miss the truncation factor by the same amount, so swaps are internally consistent — but only by accident. Once C1 is fixed by adding the normalizer to `lik_sampletimes` (or any `logLik*` slot), swap_heats will pick it up correctly. But if the normalizer is instead implemented as a proposal-only correction (current `lastneg_logratio_correction` pattern), the swap will silently miss it.

**Fix:** add the truncation normalizer to a `logLik*` slot, not just as a proposal correction. Verify by inspecting `swap_heats(heats, likelihood)` after the fix to confirm `likelihood` includes the new term.

### I14. Missing NA-guard in sanity loop
**File:** [R/mcmc-environment-functions.R:117-133](R/mcmc-environment-functions.R#L117-L133)

The first sanity loop uses `while (node != 0 && node != i)` without NA-guard, while the second loop at line 141 explicitly tests `!is.na(anc)`. If any `v$nodeparents` entry is `NA_integer_` (possible during partially-initialized trees), the first loop throws `"missing value where TRUE/FALSE needed"`.

**Fix:** mirror the second loop: `while (!is.na(node) && node != 0 && node != i)`.

### I15. Silent abort in `update_host_keepphylo` without diagnostic counter
**File:** [R/mcmc-updatehost-paths.R:75-88](R/mcmc-updatehost-paths.R#L75-L88)

When `propose_tinf_lastneg` returns `NA`, the function exits mid-prepared with no accept/reject record. The rejection counter is not updated. Combined with up to 1000 rejection retries, this can dominate runtime if `last.negative` is tight, with no visible diagnostic.

**Fix:** register the abort as a rejected proposal in the existing diagnostic counter, or log a warning at most once per N rejections.

### I16. `plotPhyloTrans` `xhost1` semantic confusion
**File:** [R/plotPhyloTrans.R:566-568](R/plotPhyloTrans.R#L566-L568)

`xhost1[valid_lastneg] <- pmin(xhost1[valid_lastneg], xlastneg_for_hosts[valid_lastneg])` silently extends the host's shaded region back to the last-negative date. `xhost1` is used downstream as the host's *infection time* (lower bound of the host bar). Code that reads `xhost1` as an infection-time coordinate after this mutation will be wrong.

**Fix:** introduce a separate `xhost1_shade` for the shaded-region plotting, and leave `xhost1` with its infection-time semantic intact.

### I17. `as.Date(numeric, origin = NULL)` errors when `xaxis.breaks` is supplied
**File:** [R/plotPhyloTrans.R:741-759](R/plotPhyloTrans.R#L741-L759)

`reference_date` defaults to NULL but `as.Date(tick_positions, origin = reference_date)` requires a non-NULL origin. Any call to `plotPhyloTrans()` that supplies `xaxis.breaks` but not the new `reference_date` arg crashes.

**Fix:** default `reference_date` to `"1970-01-01"` or pull from `plotinput$d$reference.date` when NULL.

### I18. Numeric x-axis always formatted as date
**File:** [R/plotPhyloTrans.R:754](R/plotPhyloTrans.R#L754)

In the numeric branch, `format(as.Date(tick_positions, origin = reference_date), "%b %d")` always formats numeric tick positions as dates. If the user passes purely numeric times (no Date class), tick labels become nonsensical date offsets ("Jan 01" + offset).

**Fix:** in the numeric branch, use `labels = tick_positions` (or accept a `tick_format` argument).

---

## Minor issues (cleanup / style)

- **Leftover debug output:** `cat()` calls in [R/mcmc-environment-functions.R](R/mcmc-environment-functions.R) at lines 55, 134, 151, 182, 251, 318, 329; `std::cout` blocks in [src/likseqenv.cpp](src/likseqenv.cpp).
- **Stale roxygen in [R/infect_distribution.R](R/infect_distribution.R):** `@param lastneg.time`, `vars`, `p` are documented but not in the signature; `le` is in the signature but undocumented.
- **[tests/testthat/test-lastneg.R](tests/testthat/test-lastneg.R):** only re-implements `1 − pgamma(...)` against itself — zero coverage of package functions. Replace with tests that actually call `propose_tinf_lastneg` and `lastneg_logratio_correction`.
- **Typo in [R/phybreakdata.R:82](R/phybreakdata.R#L82):** `"shoud"` and missing close-quote on `"Date"`.
- **[R/mcmc-updatehost-paths.R:22](R/mcmc-updatehost-paths.R#L22):** `d$last.negative[hostID]` relies on the orderedhosts invariant; make explicit (`d$last.negative[d$hostnames[hostID]]`) or document.
- **Missing trailing newlines:** [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R), [R/treetransformations.R](R/treetransformations.R).
- **Unreachable branch in [R/modules.R:698-702](R/modules.R#L698-L702)** (pre-existing, but touched by diff).
- **Perfect-test assumption undocumented.** Note in [R/phybreakdata.R:22](R/phybreakdata.R#L22) that `last.negative` is treated as a perfect-sensitivity exact observation.
- **Package-global `tinf.lastneg.maxtries`.** [R/mcmc-updatehost-paths.R:62](R/mcmc-updatehost-paths.R#L62) — top-level mutable; a user `assign("tinf.lastneg.maxtries", 5, envir = .GlobalEnv)` would silently change MCMC behavior. Move into a constant inside `propose_tinf_lastneg` or expose via a function argument.
- **Degenerate Gamma proposal at low `sample.shape`.** [R/mcmc-updatehost-paths.R:62-95](R/mcmc-updatehost-paths.R#L62-L95) — if `sample.shape < 1.5`, then `tinf.prop.shape.mult * sample.shape < 1` and the Gamma proposal becomes degenerate near 0, inflating rejection. Add a check or document.
- **Snake_case vs dot.case naming.** [R/plotPhyloTrans.R:200-204](R/plotPhyloTrans.R#L200-L204) — new arg `reference_date` uses snake_case while the rest of the package uses `reference.date`. Pick one.
- **Roxygen formatting.** [R/phybreakdata.R:22](R/phybreakdata.R#L22) — `@param last.negative` is indented by 2 extra spaces, may be merged into `removal.times`'s docs by roxygen2.
- **Roxygen merge bug.** [R/plotPhyloTrans.R:65-67](R/plotPhyloTrans.R#L65-L67) — a `#'` line was joined into the previous line, breaking `@section Details`.
- **`update_pathI` proposal symmetry note.** Once C11 is fixed (`Z` normalizers for path I), document the asymmetry in the proposal kernel explicitly so future maintainers don't reintroduce the bug.

---

## Recommended order to address

1. **C1** — add the truncation normalizer to the target distribution. The rest only makes sense once the target is right.
2. **C8** — rewrite `update_mS`/`update_mG` as proper MH that include the normalizer's derivative w.r.t. `sample.mean`/`sample.shape`. Without this, the cure is incomplete: parameter moves still drift to biased values.
3. **C9** — hoist `update_move_sampleedges()` out of the new nested guards (regression in within-host mixing).
4. **C2, C10, C11** — fix path I/J `tinf.prop` handling and add `Z` normalizer for path I.
5. **C3, C4, C14** — gate paths C/F under `last.negative` and enforce the constraint in within-host rewires.
6. **C12, C13** — fix the C++ negative-clamp `else` branch and stop running `lik_func` on corrupt state.
7. **I2** — replace rejection sampling with one-shot inverse-CDF truncated Gamma draw.
8. **C5, I10** — stop fabricating likelihood components on NaN; investigate the upstream cause (the `eps=1e-6` coalescent-time shift heuristic).
9. **C7, I3** — fix the broken test and uncomment `admission.times` plumbing.
10. **I13** — verify that, after C1, the truncation normalizer flows through `swap_heats` correctly (it should, if added to a `logLik*` slot rather than as a pure proposal correction).
11. **I1, I5–I9, I11, I14–I18** — clean up guards, dead/unsafe code, plot regressions.
12. **Calibration check:** simulate with known generation/sampling parameters, set `last.negative[i] = t_inf[i] − k_i` for known `k_i > 0`, fit, and check that posteriors for `sample.mean`/`sample.shape` cover truth. Without C1 + C8, this calibration will fail; with the full fix, it should pass.
13. Remove all debug `cat()` / `std::cout` before merging.

---

## Method notes

This report consolidates findings from two rounds of specialized review (R-code bug review and statistical methods review, each run twice with the second pass instructed to find issues beyond the first) plus inspection of the diff `master...HEAD`. Both reviewers independently identified the missing truncation-normalizer (C1) as the headline issue, from different angles (Hastings-ratio math and likelihood specification). The second pass surfaced the parameter-update consequences of C1 (C8 — `update_mS`/`update_mG` Gibbs misspecification), the within-host-mixing regression (C9), the path I/J/K omissions (C10, C11), the C++ `else`-branch clamp gap (C12), the corrupt-state `lik_func` propagation (C13), and the kinetic boundary-leakage characterization (C14). All file:line references were verified against the current branch state.
