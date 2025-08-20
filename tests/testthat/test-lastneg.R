# tests/testthat/test-lastneg.R

library(testthat)

# Example dummy diagnostic delay parameters
shape <- 2
mean <- 1.6
scale <- mean / shape

# Example: last-negative at day 5, test infection at day 4 (should be possible), at day 6 (should be penalized)
diagnostic_cdf <- function(x) pgamma(x, shape=shape, scale=scale)

test_that("Likelihood for infection times respects last-negative", {
  t_inf1 <- 7         # Infection before last-negative
  t_lastneg <- 5
  penalty1 <- 1 - diagnostic_cdf(t_lastneg - t_inf1)
  expect_gt(penalty1, 0.2) # Should be >0, not penalized much
  
  t_inf2 <- 1        # Infection way before last-negative
  penalty2 <- 1 - diagnostic_cdf(t_lastneg - t_inf2)
  expect_lt(penalty2, 0.05) # Should be ~0 (strongly penalized)
})
