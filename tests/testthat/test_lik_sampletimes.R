library(phybreak)
context("Test lik_sampletimes")

test_that("lik_sampletimes respects last.neg", {
  shape <- 2; mean <- 4
  
  samp <- as.numeric(as.Date(c("2020-01-10", "2020-01-11")))
  inf  <- as.numeric(as.Date(c("2020-01-05", "2020-01-05")))
  lastneg <- as.numeric(as.Date(c("2020-01-10", NA)))
  
  # Case‐1: for host 1, D1 = samp[1]−inf[1] = 5 days,
  # but M1 = samp[1]−lastneg[1] = 1 day,
  # so a Gamma‐density on [0,∞) truncated to [0,1] must assign zero likelihood
  ll1 <- lik_sampletimes(1, shape, mean, samp, inf, lastneg)
  expect_true(is.infinite(ll1) && ll1 < 0)
  
  # Case‐2: if we move lastneg back so that M1 = 6 days,
  # then both hosts have feasible intervals and ll should be finite
  lastneg2 <- as.numeric(as.Date(c("2020-01-04", NA)))
  ll2 <- lik_sampletimes(2, shape, mean, samp, inf, lastneg2)
  expect_true(is.finite(ll2))
})

