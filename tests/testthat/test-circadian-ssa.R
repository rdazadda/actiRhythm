# Singular spectrum analysis.

test_that("circadian.ssa separates trend, circadian, and noise", {
  set.seed(1)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 1800, length.out = 10 * 48)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  trend  <- 50 + 0.05 * seq_along(ts)
  circ   <- 60 * cos(2 * pi * (h - 14) / 24)
  counts <- trend + circ + stats::rnorm(length(ts), 0, 8)

  ssa <- circadian.ssa(counts, ts)
  expect_s3_class(ssa, "actiRhythm_ssa")
  expect_false(ssa$insufficient)
  expect_equal(sum(ssa$partial_variances), 1, tolerance = 0.01)
  expect_true(ssa$fundamental_period >= 22 && ssa$fundamental_period <= 26)
  expect_gt(stats::cor(ssa$circadian, circ), 0.75)
  expect_gt(stats::cor(ssa$trend, trend), 0.85)
})

test_that("circadian.ssa partial variances are sorted and the wcor diag is 1", {
  set.seed(3)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 1800, length.out = 8 * 48)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- 80 + 40 * cos(2 * pi * (h - 10) / 24) + stats::rnorm(length(ts), 0, 6)
  ssa <- circadian.ssa(counts, ts)
  expect_true(all(diff(ssa$partial_variances) <= 1e-9))
  expect_equal(diag(ssa$wcor), rep(1, ncol(ssa$wcor)))
})

test_that("circadian.ssa returns insufficient on short data and never errors", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 1800, length.out = 24)
  r <- circadian.ssa(stats::rnorm(24, 100, 10), ts)
  expect_true(r$insufficient)
  expect_no_error(print(r))
})
