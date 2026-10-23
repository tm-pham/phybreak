library(phybreak)
context("Module parameter updaters run during burn-in")

# Regression test for R/burnin_phybreak.R: the dispatch to the module updaters
# (update_wh_history, update_ir, update_R, ...) fetched the function from
# x$updaters without calling it, so intro.rate, wh.history and R stayed at
# their initial values for the whole burn-in and only started moving once
# sample_phybreak() was run.

set.seed(20241)
sim <- sim_phybreak(obsize = 6, sample.mean = 3, sample.shape = 3,
                    gen.mean = 3, gen.shape = 3)
dat <- phybreakdata(sequences = sim$sequences, sample.times = sim$sample.times)

x <- phybreak(dat, sample.mean = 3, gen.mean = 3,
              multiple.introductions = TRUE,
              intro.rate = 0.1, wh.history = 50,
              prior.intro.rate.mean = 0.05, prior.intro.rate.shape = 2,
              prior.wh.history.mean = 20, prior.wh.history.shape = 2)

test_that("introductions module is enabled with the chosen initial values", {
  expect_true(x$p$mult.intro)
  expect_true(x$h$est.ir)
  expect_true(x$h$est.wh.h)
  expect_equal(x$p$intro.rate, 0.1)
  expect_equal(x$p$wh.history, 50)
})

test_that("burn-in moves intro.rate and wh.history away from their initial values", {
  set.seed(20242)
  y <- burnin_phybreak(x, ncycles = 300, nchains = 1)
  expect_false(isTRUE(all.equal(y$p$intro.rate, x$p$intro.rate)))
  expect_false(isTRUE(all.equal(y$p$wh.history, x$p$wh.history)))
})
