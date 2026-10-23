library(phybreak)
context("Initial state has a finite log-likelihood")

# Regression test for the NJ-tree initialisation of the history host
# (R/treetransformations.R, R/mcmc-sample-minitrees.R). For some datasets
# phybreak() built an initial state whose sequence log-likelihood was NaN:
# sample_coaltimes_history(init = TRUE) shifted the sampled history-host
# coalescent times so that the latest one coincided with the earliest
# infection time, and whether its strict "children later than parent" test
# then passed depended on floating-point rounding. When it failed, the
# failure flag was ignored and the coalescent nodes kept the times assigned
# earlier under a different topology, giving negative branch lengths and
# NaN partial likelihoods. burnin_phybreak() and sample_phybreak() then ran to
# completion rejecting every move (one mh_accept() warning per move) and
# returned the initial state as if it were a posterior sample.
#
# The dataset below reproduces this: five introductions spread over 425 days
# with a cluster-wide mutation rate, generated with the L'Ecuyer RNG. The
# initial tree is random, so a range of initialisation seeds is tried.

old_rng <- RNGkind("L'Ecuyer-CMRG")
set.seed(3005)
sim <- sim_phybreak(obsize = 6, introductions = 5, intro.rate = 0.01, R0 = 2,
                    gen.shape = 5, gen.mean = 4, sample.shape = 3, sample.mean = 2,
                    wh.model = "linear", wh.history = 50, mu = 2.2e-6,
                    sequence.length = 29903)

make_fit <- function(seed) {
  set.seed(seed)
  phybreak(dataset = sim, gen.shape = 5, gen.mean = 4, sample.shape = 3,
           sample.mean = 2, est.gen.mean = FALSE, est.sample.mean = FALSE,
           wh.model = "linear", introductions = 1, multiple.introductions = TRUE)
}

test_that("every initialisation gives a finite log-likelihood", {
  logliks <- sapply(1:20, function(seed) as.numeric(logLik(make_fit(seed))))
  expect_true(all(is.finite(logliks)), 
              info = paste("non-finite for seeds:", paste(which(!is.finite(logliks)), collapse = " ")))
})

test_that("the initial tree has no negative branch lengths", {
  for (seed in 1:5) {
    v <- make_fit(seed)$v
    haspar <- v$nodeparents > 0
    expect_true(all(v$nodetimes[haspar] - v$nodetimes[v$nodeparents[haspar]] >= 0))
  }
})

test_that("a short two-chain burn-in from this start emits no mh_accept warnings", {
  f <- make_fit(1)
  w <- character(0)
  set.seed(42)
  y <- withCallingHandlers(
    suppressMessages(burnin_phybreak(f, ncycles = 30, nchains = 2)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("NA/NaN", w)), 0, info = paste(head(unique(w)), collapse = "\n"))
  expect_true(is.finite(logLik(y)))
})

test_that("the samplers refuse to start from a non-finite log-likelihood", {
  f <- make_fit(6)
  f$p$mu <- NaN   # makes logLikseq NaN
  expect_error(suppressMessages(burnin_phybreak(f, ncycles = 1)), "non-finite")
  expect_error(suppressMessages(sample_phybreak(f, nsample = 1)), "non-finite")
})

do.call(RNGkind, as.list(old_rng))
