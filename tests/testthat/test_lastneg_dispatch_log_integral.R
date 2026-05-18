library(phybreak)
context("Test lastneg_dispatch_log_integral")

# Helpers
make_pv <- function(sample.shape, sample.mean, nodetime = 10) {
  p <- list(sample.shape = sample.shape, sample.mean = sample.mean)
  v <- list(nodetimes = c(nodetime))
  list(p = p, v = v)
}

shape.prop <- function(sample.shape) phybreak:::tinf.prop.shape.mult * sample.shape

test_that("no last-negative: matches pgamma exactly", {
  pv <- make_pv(sample.shape = 3, sample.mean = 6)
  d_null <- list(last.negative = NULL)
  d_na   <- list(last.negative = NA_real_)

  for (interval in c(0.5, 1.5, 4.0, 12.0)) {
    sp <- shape.prop(pv$p$sample.shape)
    expected <- pgamma(interval, shape = sp, scale = pv$p$sample.mean / sp, log.p = TRUE)
    expect_equal(
      phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d_null),
      expected
    )
    expect_equal(
      phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d_na),
      expected
    )
  }
})

test_that("last-negative with interval <= M: matches pgamma", {
  pv <- make_pv(sample.shape = 3, sample.mean = 6, nodetime = 10)
  # M = nodetime - lastneg = 10 - 4 = 6
  d <- list(last.negative = 4)
  sp <- shape.prop(pv$p$sample.shape)

  for (interval in c(0.5, 2.0, 5.999)) {
    expected <- pgamma(interval, shape = sp, scale = pv$p$sample.mean / sp, log.p = TRUE)
    expect_equal(
      phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d),
      expected
    )
  }
})

test_that("interval > M: agrees with direct numerical integration", {
  pv <- make_pv(sample.shape = 4, sample.mean = 5, nodetime = 10)
  # M = 10 - 6 = 4
  d <- list(last.negative = 6)
  sp <- shape.prop(pv$p$sample.shape)
  sf <- pv$p$sample.shape
  sc_prop <- pv$p$sample.mean / sp
  sc_full <- pv$p$sample.mean / sf
  M <- 4

  for (interval in c(5.0, 8.0, 15.0)) {
    direct <- integrate(function(D) {
      dgamma(D, shape = sp, scale = sc_prop) *
        (1 - pgamma(D - M, shape = sf, scale = sc_full))
    }, lower = 0, upper = interval)$value
    expected <- log(direct)
    got <- phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d)
    expect_equal(got, expected, tolerance = 1e-6)
  }
})

test_that("M <= 0 (lastneg >= nodetime): purely numerical, no closed-form head", {
  pv <- make_pv(sample.shape = 2, sample.mean = 4, nodetime = 10)
  # lastneg = 10 means M = 0; lastneg = 12 means M = -2
  sp <- shape.prop(pv$p$sample.shape)
  sf <- pv$p$sample.shape
  sc_prop <- pv$p$sample.mean / sp
  sc_full <- pv$p$sample.mean / sf

  for (lastneg in c(10, 12)) {
    d <- list(last.negative = lastneg)
    M <- pv$v$nodetimes[1] - lastneg
    interval <- 6.0
    direct <- integrate(function(D) {
      dgamma(D, shape = sp, scale = sc_prop) *
        (1 - pgamma(D - M, shape = sf, scale = sc_full))
    }, lower = 0, upper = interval)$value
    expected <- log(direct)
    got <- phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d)
    expect_equal(got, expected, tolerance = 1e-6)
  }
})

test_that("result is bounded above by pgamma(interval) and below by pgamma(M)", {
  # When last-negative is active and interval > M:
  #   pgamma(M) <= P_distorted(D < interval) <= pgamma(interval)
  pv <- make_pv(sample.shape = 3, sample.mean = 7, nodetime = 10)
  d <- list(last.negative = 5)  # M = 5
  sp <- shape.prop(pv$p$sample.shape)
  sc_prop <- pv$p$sample.mean / sp
  M <- 5
  interval <- 9.0

  got <- phybreak:::lastneg_dispatch_log_integral(interval, 1, pv$p, pv$v, d)
  upper <- pgamma(interval, shape = sp, scale = sc_prop, log.p = TRUE)
  lower <- pgamma(M, shape = sp, scale = sc_prop, log.p = TRUE)

  expect_true(got <= upper + 1e-12)
  expect_true(got >= lower - 1e-12)
})

test_that("interval <= 0 returns -Inf", {
  pv <- make_pv(sample.shape = 3, sample.mean = 6)
  d <- list(last.negative = 4)
  expect_identical(phybreak:::lastneg_dispatch_log_integral(0, 1, pv$p, pv$v, d), -Inf)
  expect_identical(phybreak:::lastneg_dispatch_log_integral(-1, 1, pv$p, pv$v, d), -Inf)
})
