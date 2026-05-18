library(phybreak)
context("Test lik_sampletimes")

# Helpers --------------------------------------------------------------------
# Plain Gamma log-density for one host's sampling interval, no last-negative.
ll_gamma <- function(D, shape, mean) {
  dgamma(D, shape = shape, scale = mean / shape, log = TRUE)
}

# Closed-form truncated-Gamma log-density: dgamma(D) / pgamma(M, ...).
ll_trunc <- function(D, M, shape, mean) {
  dgamma(D, shape = shape, scale = mean / shape, log = TRUE) -
    pgamma(M, shape = shape, scale = mean / shape, log.p = TRUE)
}

# Hard-constraint cases (D >= M) ---------------------------------------------

test_that("D >= M (constraint violated) returns -Inf", {
  shape <- 2; mean <- 4

  samp <- as.numeric(as.Date(c("2020-01-10", "2020-01-11")))
  inf  <- as.numeric(as.Date(c("2020-01-05", "2020-01-05")))
  # M_1 = samp[1] - lastneg[1] = 0 days; D_1 = samp[1] - inf[1] = 5 days
  # D_1 >= M_1 so the constraint t_inf > lastneg is violated -> -Inf.
  lastneg <- as.numeric(as.Date(c("2020-01-10", NA)))

  ll1 <- lik_sampletimes(1, shape, mean, samp, inf, lastneg)
  expect_true(is.infinite(ll1) && ll1 < 0)
})

test_that("D < M (constraint satisfied) is finite and equals the closed form", {
  shape <- 2; mean <- 4

  samp <- as.numeric(as.Date(c("2020-01-10", "2020-01-11")))
  inf  <- as.numeric(as.Date(c("2020-01-05", "2020-01-05")))
  # M_1 = samp[1] - lastneg[1] = 6 days; D_1 = 5 days. Constraint OK.
  # Host 2 has NA last-negative -> plain dgamma contribution.
  lastneg <- as.numeric(as.Date(c("2020-01-04", NA)))

  ll2 <- lik_sampletimes(2, shape, mean, samp, inf, lastneg)
  expect_true(is.finite(ll2))

  # Expected: ll_trunc for host 1 (D=5, M=6) + ll_gamma for host 2 (D=6).
  D <- samp - inf
  M1 <- samp[1] - lastneg[1]
  expected <- ll_trunc(D[1], M1, shape, mean) + ll_gamma(D[2], shape, mean)
  expect_equal(ll2, expected, tolerance = 1e-12)
})

# Truncation normalizer is actually subtracted -------------------------------

test_that("normalizer -log F(M) is applied (truncated-Gamma, not just indicator)", {
  shape <- 3; mean <- 5
  samp <- 20
  inf  <- 13     # D = 7
  lastneg <- 15  # M = 5; D > M would violate. Move lastneg earlier:
  lastneg <- 10  # M = 10; D = 7 < M, constraint OK.

  got <- lik_sampletimes(1, shape, mean, samp, inf, lastneg)
  expected <- ll_trunc(samp - inf, samp - lastneg, shape, mean)
  expect_equal(got, expected, tolerance = 1e-12)

  # The normalizer must make the truncated-Gamma value STRICTLY GREATER than
  # the plain Gamma (since F(M) < 1, -log F(M) > 0). If the implementation
  # were the hard indicator alone, the two would be equal.
  plain <- ll_gamma(samp - inf, shape, mean)
  expect_true(got > plain + 1e-9)
  expect_equal(got - plain, -pgamma(samp - lastneg, shape = shape, scale = mean / shape, log.p = TRUE),
               tolerance = 1e-12)
})

# No / NA / NULL last-negative reduce to plain sum(dgamma) -------------------

test_that("no last-negative reduces to sum(dgamma)", {
  shape <- 2.5; mean <- 6
  samp <- c(10, 12, 14)
  inf  <- c(3, 4, 7)
  D <- samp - inf

  expected <- sum(ll_gamma(D, shape, mean))

  expect_equal(lik_sampletimes(3, shape, mean, samp, inf), expected, tolerance = 1e-12)
  expect_equal(lik_sampletimes(3, shape, mean, samp, inf, NULL), expected, tolerance = 1e-12)
  expect_equal(lik_sampletimes(3, shape, mean, samp, inf, rep(NA_real_, 3)), expected, tolerance = 1e-12)
})

# Mixed hosts: only constrained hosts contribute normalizer ------------------

test_that("only hosts with non-NA last-negative pick up the normalizer", {
  shape <- 2; mean <- 5
  samp <- c(10, 12, 14, 16)
  inf  <- c(5, 7, 9, 11)            # D = c(5, 5, 5, 5)
  lastneg <- c(NA, 4, NA, 6)        # M = c(NA, 8, NA, 10) for present entries

  got <- lik_sampletimes(4, shape, mean, samp, inf, lastneg)

  # Hosts 1 and 3 contribute plain Gamma; 2 and 4 contribute truncated Gamma.
  D <- samp - inf
  M2 <- samp[2] - lastneg[2]
  M4 <- samp[4] - lastneg[4]
  expected <- ll_gamma(D[1], shape, mean) +
              ll_trunc(D[2], M2, shape, mean) +
              ll_gamma(D[3], shape, mean) +
              ll_trunc(D[4], M4, shape, mean)

  expect_equal(got, expected, tolerance = 1e-12)
})

# Boundary D == M -------------------------------------------------------------

test_that("boundary D == M returns -Inf (matches C1 D >= M semantics)", {
  shape <- 2; mean <- 5
  samp <- 10
  inf <- 5         # D = 5
  lastneg <- 5     # M = 5, so D == M

  expect_identical(lik_sampletimes(1, shape, mean, samp, inf, lastneg), -Inf)
})
