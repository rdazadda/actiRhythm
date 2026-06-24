# Bootstrap confidence interval for the circadian period.

test_that("period.ci recovers a known 24h period inside the CI", {
  t_hours <- seq(0, 7 * 24 - 1/60, by = 1/60)
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
  set.seed(42)
  counts <- pmax(0, 100 + 80 * cos(2 * pi * (t_hours - 8) / 24) +
                   stats::rnorm(length(t_hours), 0, 20))

  r <- period.ci(counts, ts, n_boot = 80, seed = 1)
  expect_s3_class(r, "actiRhythm_period_ci")
  expect_lt(abs(r$tau - 24), 0.5)                 # point estimate near 24 h
  expect_lt(r$ci_lower, 24)                        # CI covers the truth
  expect_gt(r$ci_upper, 24)
  expect_true(r$ci_lower <= r$ci_upper)
})

test_that("parabolic refinement gives a non-degenerate SE (not grid-snapped)", {
  t_hours <- seq(0, 6 * 24, by = 1/30)
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
  set.seed(7)
  counts <- pmax(0, 100 + 60 * cos(2 * pi * t_hours / 24) +
                   stats::rnorm(length(t_hours), 0, 15))

  r <- period.ci(counts, ts, n_boot = 60, seed = 2)
  expect_gt(r$se, 0)                               # would be 0 if peaks snapped to grid
  expect_true(is.finite(r$ci_lower) && is.finite(r$ci_upper))
})

test_that("period.ci is reproducible with a seed", {
  t_hours <- seq(0, 4 * 24, by = 1/30)
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
  set.seed(99)
  counts <- pmax(0, 100 + 60 * cos(2 * pi * t_hours / 24) +
                   stats::rnorm(length(t_hours), 0, 15))

  a <- period.ci(counts, ts, n_boot = 40, seed = 5)
  b <- period.ci(counts, ts, n_boot = 40, seed = 5)
  expect_equal(a$ci_lower, b$ci_lower)
  expect_equal(a$ci_upper, b$ci_upper)
})

test_that("period.ci returns an NA structure on insufficient data", {
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + (0:100) * 60   # < 2 days
  r <- period.ci(stats::rnorm(101), ts, n_boot = 10)
  expect_s3_class(r, "actiRhythm_period_ci")
  expect_true(is.na(r$ci_lower))
})

test_that(".lsp_peak_parabola refines toward the true peak", {
  scanned <- seq(20, 28, by = 0.5)
  power   <- exp(-((scanned - 24.1)^2) / 0.5)     # smooth peak near 24.1, off-grid
  pk <- actiRhythm:::.lsp_peak_parabola(scanned, power)
  expect_lt(abs(pk$tau - 24.1), 0.1)              # sub-grid vertex, not the 24.0 node
  expect_false(pk$at_boundary)
})
