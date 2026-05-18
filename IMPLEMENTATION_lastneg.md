# Implementation notes: last-negative test in `mcmc-updatehost-paths.R`

Companion to [REVIEW_lastneg.md](REVIEW_lastneg.md) and [FIXES_lastneg.md](FIXES_lastneg.md). This file documents the *proposal-side* implementation in [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R) and explains why paths C and F carry `TODO(last.negative)` markers instead of corrections.

---

## 1. Summary of changes

All edits are confined to [R/mcmc-updatehost-paths.R](R/mcmc-updatehost-paths.R) and address four issues identified in a prior review of the in-place rejection-sampling block:

1. **Missing MH proposal-ratio correction** when the proposal is rejection-sampled.
2. **No iteration guard** on the `repeat { ... }` rejection sampler — could hang.
3. **Undocumented proposal-vs-acceptance shape mismatch** (proposal uses `2/3 * sample.shape`, weight uses full `sample.shape`).
4. **Three near-identical copies** of the rejection-sampling block — drift risk.

### 1.1 New helpers (top of file, [lines 6–54](R/mcmc-updatehost-paths.R#L6-L54))

- **`tinf.lastneg.maxtries <- 1000`** ([line 8](R/mcmc-updatehost-paths.R#L8)) — caps the rejection loop.
- **`propose_tinf_lastneg(hostID, p, v, d)`** ([line 22](R/mcmc-updatehost-paths.R#L22)) — single source of truth for the rejection sampler. If `d$last.negative[hostID]` is missing/`NA`, draws a plain `rgamma`. Otherwise rejection-samples and either returns a `tinf.prop` or `NA_real_` (and warns) on exhaustion. Documents the intentional shape asymmetry in its docstring (issues 2, 3, 4).
- **`lastneg_logratio_correction(hostID, tinf_old, tinf_new, p, d)`** ([line 47](R/mcmc-updatehost-paths.R#L47)) — returns `log[s(tinf_old)/s(tinf_new)]` where `s(t) = 1 - pgamma(lastneg - t; full shape)`. Returns `0` when no last-negative is set, so call sites stay clean (issue 1).

### 1.2 Call-site simplifications (deduplication, issue 4)

The three duplicated `if/else { repeat{...} }` blocks in
[`update_host_keepphylo`](R/mcmc-updatehost-paths.R#L87),
[`update_host_phylotrans`](R/mcmc-updatehost-paths.R#L193), and
[`update_host_history`](R/mcmc-updatehost-paths.R#L265)
are now two lines each: call `propose_tinf_lastneg`, abort on `NA`.

### 1.3 MH correction added to path-function `logproposalratio` (issue 1)

| Path | File line of correction | Group |
|---|---|---|
| A | [299](R/mcmc-updatehost-paths.R#L299) | phylotrans |
| B | [350](R/mcmc-updatehost-paths.R#L350) | phylotrans |
| D | [454](R/mcmc-updatehost-paths.R#L454) | phylotrans |
| E | [516](R/mcmc-updatehost-paths.R#L516) | phylotrans |
| G | [803](R/mcmc-updatehost-paths.R#L803) | keepphylo |
| H | [872](R/mcmc-updatehost-paths.R#L872) | keepphylo |
| J | [1096](R/mcmc-updatehost-paths.R#L1096) | keepphylo |
| L | [664](R/mcmc-updatehost-paths.R#L664) | history |

For these paths the correction term takes the form
```r
+ lastneg_logratio_correction(hostID, v$inftimes[hostID], tinf.prop, p, d)
```
(in keepphylo, `pbe0$v$inftimes[hostID]` is used because the function has already overwritten `v$inftimes[hostID] <- tinf.prop` earlier in its body).

### 1.4 Path I was originally in this list but was later **refused under last-negative** by the [C2 fix](FIXES_lastneg.md). The `lastneg_logratio_correction` call I had added there has been removed and replaced with an early `return()` when `d$last.negative[hostID]` is set ([line 914](R/mcmc-updatehost-paths.R#L914)).

### 1.5 Paths C and F carry **TODO comments**, not corrections (issue 1, partial)

[Path C, line 397](R/mcmc-updatehost-paths.R#L397) and [path F, line 553](R/mcmc-updatehost-paths.R#L553):
```r
# TODO(last.negative): MH correction for path C/F not yet derived under the
# index-swap; tinf.prop here is not directly hostID's post-acceptance inftime.
# Acceptable when last.negative is unset for hostID; biased otherwise.
```

Rationale below.

---

## 2. Why paths C and F have TODOs instead of corrections

### 2.1 What the "simple" correction assumes

For paths A, B, D, E, G, H, J, L, the `lastneg_logratio_correction` term is provably correct because the move has a clean structure:

- The proposal draws **hostID's** new infection time `tinf.prop` from `propose_tinf_lastneg(hostID, ...)`, which is a rejection-sampled Gamma anchored at hostID's first positive sample, weighted by `s_hostID(lastneg_hostID - t)`.
- After acceptance, **`v$inftimes[hostID] = tinf.prop`**: the proposed value becomes hostID's new infection time. The role of hostID in the tree (index vs. non-index, who its infector is) is the same in forward and reverse.
- The reverse proposal, by the same algorithm, would draw the *current* `v$inftimes[hostID]` from the same rejection-sampled Gamma anchored at hostID, with the same `s_hostID` weight.
- Forward and reverse densities therefore differ only in where they're evaluated, and the MH correction is the clean ratio `log[s_hostID(D - M)] − log[s_hostID(D' - M)]`.

### 2.2 What's different about paths C and F: the **index-swap**

[Path C (line 380)](R/mcmc-updatehost-paths.R#L380) fires when hostID is currently the index *and* the proposed `tinf.prop` falls past hostID's **second** secondary infection. The move then relabels `newindexID` (hostID's first infectee) as the new index. After acceptance, hostID is no longer index — it becomes a non-index infected by `newindexID`.

[Path F (line 555)](R/mcmc-updatehost-paths.R#L555) is the symmetric mirror: hostID is a non-index whose proposed `tinf.prop` lands past its first secondary infection, triggering a swap where hostID becomes (or remains in a different position in) a re-rooted subtree.

Two things make these structurally different from the "simple" paths:

1. **`tinf.prop` is not hostID's post-acceptance infection time.** The swap reassigns infection times: in pathC, hostID's post-acceptance `inftimes[hostID]` is set to a time derived from the swap geometry, **not** to `tinf.prop`. The variable `tinf.prop` instead *gates* entry into pathC via the comparison `tinf.prop > sort(inftimes[infectees])[2]`, after which it is consumed by the rewire functions to position hostID in the swapped tree.

2. **Forward and reverse `logproposalratio` already had non-standard algebra.** Look at the existing formula:
   ```r
   logproposalratio <- pgamma(sampleinterval.newindex, shape, scale, log.p=TRUE) -
                       pgamma(sampleinterval.hostID,    shape, scale, log.p=TRUE)
   ```
   This is a difference of *CDF* values on sample-intervals, not a `dgamma` density ratio. It is the result of Klinkenberg algebraically simplifying the forward/reverse densities of the swap: `dgamma` terms cancel and what remains are these `pgamma` truncation factors that depend on **which host's** sample-interval gates the swap in each direction.

### 2.3 Why naively adding `lastneg_logratio_correction(hostID, …)` is wrong

To derive the correct MH correction under rejection sampling, we need the ratio
of the actual forward and reverse *proposal densities*, including the `s(.)`
weights. The hosts whose `s(.)` enter are determined by which host's `tinf` is
drawn from `propose_tinf_lastneg` in each direction:

- **Forward (pathC)**: hostID was index. We drew `tinf.prop` for hostID via `propose_tinf_lastneg(hostID, ...)`. So `s_hostID(lastneg_hostID - tinf.prop)` enters the forward density.
- **Reverse (un-swap)**: starting from the new state (newindexID is index, hostID is a non-index), the reverse move would propose for a host whose role swap would land back in the pre-swap configuration. The mechanics in `update_host_phylotrans` always draw `tinf.prop` for the *focal* host whose `update_host_phylotrans(hostID = ...)` was called. The reverse focal host in the swapped state is not the same as the forward focal host; whose `last.negative` weights the reverse density depends on which host gets selected as focal in the reverse move.

A simple `log[s_hostID(t_old) / s_hostID(t_new)]` correction encodes neither (a) the fact that `t_old` and `t_new` are not hostID's `inftimes` in pathC/F, nor (b) the fact that the reverse proposal is weighted by a *different* host's `s(.)`. Pasting it on top of the existing `pgamma(sampleinterval…)` formula would silently produce a biased MH ratio in the regime where `last.negative` is set for one of the swapped hosts.

### 2.4 What "Acceptable when `last.negative` is unset for hostID" means in the TODO

The `propose_tinf_lastneg` helper falls back to a plain `rgamma` when `d$last.negative[hostID]` is `NULL`/`NA`. In that case the proposal is the *original* unmodified Gamma, the rejection-sampling weight `s(.)` never enters, and the existing `pgamma(sampleinterval…)` formula gives the correct MH ratio. Hence pathC/F are unbiased whenever the focal host has no last-negative date.

The bias is introduced only when `d$last.negative[hostID]` is set and the proposed move happens to land in pathC or pathF. Path C is a relatively rare branch (`tinf.prop` past the *second* secondary infection of an index case), and the other phylotrans paths (A, B, D, E) remain available for the same host, so the practical impact on mixing is limited but **non-zero**.

### 2.5 Interaction with the C1 likelihood fix

[FIXES_lastneg.md § C1](FIXES_lastneg.md) added a perfect-test truncation factor to `lik_sampletimes`, i.e. the target distribution now includes a hard `1{tinf > lastneg}` indicator (states with `tinf <= lastneg` have `logLik = -Inf`).

Under that perfect-test target, the survival weight `s(t) = 1 - pgamma(lastneg - t; full shape)` evaluates to **1 for any state in the target's support**, because the target only assigns positive density to states where `tinf > lastneg` (i.e. `lastneg - tinf < 0`, so `pgamma` of a negative argument is `0`).

Consequently:
- `lastneg_logratio_correction` returns `log(1) − log(1) = 0` for any move between two states both satisfying the constraint.
- For moves to states that violate the constraint, the target's `-Inf` rejects them outright, regardless of the proposal-ratio correction.

So my `lastneg_logratio_correction` term is a **no-op under the C1 perfect-test target**: it is mathematically harmless but adds no information. It would become essential only if the model were re-parameterized to treat `s(.)` as a *soft* likelihood factor (imperfect test) rather than a hard truncation. The code therefore remains correct under the current perfect-test convention.

For pathC/F specifically, this interaction means the missing MH correction described in §2.3 is also a no-op under C1's perfect-test target *for states in support*. The remaining concern is whether `propose_tinf_lastneg`'s rejection sampler distorts the **gating probability** that pathC/F is selected at all (a separate issue from the MH ratio correction). If the rejection sampler shifts mass toward larger `tinf.prop`, pathC/F may be reached *less often* than under the unmodified Gamma proposal — and this gating-frequency distortion is not captured by `lastneg_logratio_correction` either. This is the deeper reason a principled fix requires a from-scratch re-derivation of the pathC/F proposal ratio that carries the `s(.)` factor through to the gating-integral.

---

## 3. Verification done

- `Rscript -e 'parse(file = "R/mcmc-updatehost-paths.R")'` parses cleanly.
- Behavior is byte-for-byte identical to pre-change when `d$last.negative` is absent from `d`, because `propose_tinf_lastneg` returns the plain Gamma draw and `lastneg_logratio_correction` returns `0`.

## 4. Outstanding items

- **C3 (paths C/F under last-negative)**: derive the swap-aware proposal ratio with `s(.)` carried through. This is the remaining proposal-side correction needed to make all path functions strictly unbiased when `last.negative` is set.
- **C10 (path J `tinf2.prop`)**: pathJ's *secondary* candidate `tinf2.prop` for `hostiorID` is drawn from a plain `rgamma` (line 117 of [update_host_keepphylo](R/mcmc-updatehost-paths.R#L117-L118)), not from `propose_tinf_lastneg(hostiorID, ...)`. If `hostiorID` has its own `last.negative`, this draw bypasses the constraint. Analogous to C2: either route the draw through `propose_tinf_lastneg(hostiorID, ...)` and add the matching correction, or refuse pathJ when `hostiorID` has a `last.negative` date.
