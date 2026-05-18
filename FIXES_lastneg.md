# Fix log: `add-lastneg-likelihood` branch

Companion to [REVIEW_lastneg.md](REVIEW_lastneg.md). Each section documents one issue from the review and the corresponding fix applied to the code.

---

## C1 — Truncation normalizer added to the target distribution

**Date:** 2026-05-17
**Issue:** [REVIEW_lastneg.md § C1](REVIEW_lastneg.md)
**Severity:** Critical
**Status:** Fixed

### What was wrong

`lik_sampletimes` was the unmodified Gamma density on `nodetime − inftime`. The last-negative-test constraint was implemented only at the proposal level (`propose_tinf_lastneg`, `lastneg_logratio_correction`), never in the joint log-posterior. MCMC therefore converged to neither the constrained nor the unconstrained target, breaking detailed balance and biasing `sample.mean` / `sample.shape` toward shorter intervals.

### Math

Let `D_i = nodetime_i − tinf_i` be host *i*'s sampling interval, modelled as `Gamma(shape, scale = mean/shape)`. A last-negative test at `lastneg_i` under perfect-test assumption forces `tinf_i > lastneg_i`, equivalently `D_i < M_i` where `M_i = nodetime_i − lastneg_i`.

The truncated Gamma density on `(0, M_i)` is

```
f(D_i) / F(M_i, shape)        for D_i ∈ (0, M_i)
0                              elsewhere
```

so the log-likelihood contribution per host with a last-negative date is

```
dgamma(D_i, shape, scale, log = TRUE) − pgamma(M_i, shape, scale, log.p = TRUE)   if D_i < M_i
-Inf                                                                              otherwise
```

The `−log F(M_i, shape)` normalizer is what makes inference of `sample.mean` / `sample.shape` unbiased. Without it, the constraint `D_i < M_i` is absorbed into the data without the corresponding parameter penalty for parameters that put little mass in `(0, M_i)`.

**Note:** the reviewer suggested `−log(1 − pgamma(lastneg − inftimes, …))`, which is algebraically inconsistent with both the perfect-test assumption and the existing test (which expects `-Inf` when `D > M`). The form implemented here `−log F(M, shape)` matches both.

### Files changed

- [R/logLik_phybreak.R](R/logLik_phybreak.R)
  - **`lik_sampletimes()`** (lines 101–129): added `last.negative = NULL` argument. For each host with a finite `last.negative`, returns the truncated Gamma log-density: `-Inf` when `D >= M`, otherwise subtracts `pgamma(M, shape, scale, log.p = TRUE)` from the standard `dgamma` contribution.
  - **`logLik.phybreak()`** (line 51): call site updated to pass `d$last.negative`.
- [R/mcmc-environment-functions.R](R/mcmc-environment-functions.R)
  - **`build_pbe()`** (line 192): call site updated to pass `d$last.negative`.
  - **`propose_pbe()`** (line 327): call site updated to pass `d$last.negative`.

### Side-effect fixes

- **C7 (broken test signature)** is now resolved: `tests/testthat/test_lik_sampletimes.R` calls `lik_sampletimes(..., lastneg)` and that signature now exists.

### Verification

- `test_lik_sampletimes.R` passes both cases:
  - Case 1: `D = 5, M = 0` → `-Inf` (constraint violated). ✓
  - Case 2: `D = [5, 6], M = [6, NA]` → finite. ✓
- `tests/testthat/test-lastneg.R` (2 tests) and `tests/testthat/test_phybreakdata.R` (30 tests) still pass.
- The one failure in `tests/testthat/test_phybreak.R` (`use.tree works correctly`) is pre-existing on this branch and unrelated to C1.

### Open follow-ups

C1 only fixes the *target distribution*. The proposal-side issues remain open and are still required for a correct sampler:

- **C8** — `update_mS` / `update_mG` Gibbs/MH ratios do not include the new normalizer's derivative w.r.t. `sample.mean` / `sample.shape`. Without this, parameter updates still drift to biased values.
- **C2, C10, C11** — path I/J `tinf.prop` handling and path I missing `Z` normalizer.
- **C3, C4** — paths C/F and within-host rewires bypass the constraint.
- **C14** — kinetic boundary leakage; auto-resolved once all constraint checks are in place.

These are the recommended next fixes (see [REVIEW_lastneg.md § Recommended order](REVIEW_lastneg.md)).

---

## C2 — Path I refuses moves under last-negative

**Date:** 2026-05-17
**Issue:** [REVIEW_lastneg.md § C2](REVIEW_lastneg.md)
**Severity:** Critical
**Status:** Fixed

### What was wrong

In `update_host_keepphylo` ([R/mcmc-updatehost-paths.R:125](R/mcmc-updatehost-paths.R#L125)), `tinf.prop <- v$inftimes[hostiorID]` overwrote the value drawn by `propose_tinf_lastneg` before dispatching to `update_pathI()`. The downstream `lastneg_logratio_correction(hostID, …, tinf.prop, p, d)` at line 979 then received the infector's existing infection time (a deterministic value), not a sample from the rejection-distorted proposal. Detailed balance was broken whenever `d$last.negative[hostID]` was set.

### Choice of fix

The review suggested two options: carry the original sampled value to the correction, or skip the correction. Both are imperfect because they miss the *gating-integral* distortion:

- In path I, hostID's new infection time is **deterministic** (`v$inftimes[hostiorID]`), not the discarded `tinf.prop_orig`.
- The discarded sample only *gates* entry into path I via `1{tinf.prop_orig < timemrca}`.
- Under last-negative, the gating probability `P(tinf.prop_orig < timemrca)` is distorted by the rejection sampler in `propose_tinf_lastneg`.
- The forward/reverse Hastings ratio therefore needs the ratio of these gating integrals, which is a numerical integral over the proposal density times the survival weight `s(t)`. Neither option presented in the review computes this.

The cleanest correctness-preserving fix is to **refuse path I when hostID has a finite last-negative date**. Path I is a specific NNYY tree structure (hostID not index, `tinf.prop < timemrca`, `tinf2.prop > timemrca`, hostiorID is index) — a relatively rare move even without last-negative — and the other paths (A, B, D, E, G, H, J) remain available for the affected hosts, so the impact on mixing is small.

### Files changed

- [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R)
  - **`update_pathI()`** (new early-return after local-variable setup): refuses the move when `length(d$last.negative) > 0 && !is.na(d$last.negative[hostID])`.
  - **`update_pathI()`** (proposal-ratio block): removed the trailing `+ lastneg_logratio_correction(hostID, pbe0$v$inftimes[hostID], tinf.prop, p, d)` term (now unreachable for hosts with last-negative, and the previous call passed the wrong `tinf.prop`). Updated the surrounding comment.
  - **`update_host_keepphylo()`** (comment near line 86): updated to note path I refuses under last-negative.

### Verification

- All `test_lik_sampletimes.R`, `test-lastneg.R`, and `test_phybreakdata.R` tests still pass.
- The pre-existing failure in `test_phybreak.R::"use.tree works correctly"` is unchanged and unrelated.
- No new code path is reachable when `d$last.negative` is NULL or all-NA, so previously-correct behavior is preserved bit-for-bit in that case.

### Open follow-ups

- **C10** (Path J missing `lastneg_logratio_correction` for `hostiorID`'s `tinf2.prop`) is structurally similar but distinct: in path J, `tinf2.prop` for `hostiorID` is drawn from a plain `rgamma` (no last-negative constraint) and used directly. The analogous fix is either to draw `tinf2.prop` via `propose_tinf_lastneg(hostiorID, …)` with the matching correction, or to refuse path J when `hostiorID` has a last-negative date. Recommend handling this when C10 is taken up.
- **C11** (path I/L `Z` normalizer): path I is now refused under last-negative so the `Z` issue no longer applies to it. Path L still needs handling.
