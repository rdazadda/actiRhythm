# Functional linear model of the 24-hour profile.

test_that("circadian.flm recovers a known multi-harmonic profile", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 7 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- 100 + 50 * cos(2 * pi * (h - 8) / 24) + 20 * cos(2 * pi * (h - 3) / 12)

  flm <- circadian.flm(counts, ts, n_harmonics = 4, weights = "none")
  expect_s3_class(flm, "actiRhythm_flm")
  expect_gt(flm$r_squared, 0.999)
  expect_equal(nrow(flm$smooth_curve), 1440L)
})

test_that("circadian.flm with one harmonic matches the cosinor", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- 100 + 40 * cos(2 * pi * (h - 16) / 24)

  flm <- circadian.flm(counts, ts, n_harmonics = 1)
  cos <- cosinor.analysis(counts, ts)
  expect_equal(flm$harmonics$amplitude[1], unname(cos$amplitude), tolerance = 0.03)
  expect_equal(flm$harmonics$acrophase_hours[1], unname(cos$acrophase), tolerance = 0.1)
})

test_that("circadian.flm R-squared is non-decreasing in n_harmonics", {
  set.seed(2)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 7 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- pmax(0, 100 + 60 * cos(2 * pi * (h - 14) / 24) +
                    30 * cos(2 * pi * (h - 2) / 8) + stats::rnorm(length(ts), 0, 10))
  r2 <- vapply(1:5, function(k) circadian.flm(counts, ts, n_harmonics = k)$r_squared, numeric(1))
  expect_true(all(diff(r2) >= -1e-8))
})

test_that("circadian.flm bspline runs and edge cases never error", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- 100 + 50 * cos(2 * pi * (h - 12) / 24)

  fb <- circadian.flm(counts, ts, basis = "bspline")
  expect_s3_class(fb, "actiRhythm_flm")
  expect_true(is.finite(fb$r_squared))

  ts2 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 600)
  r <- circadian.flm(stats::rnorm(600, 100, 10), ts2)
  expect_true(is.na(r$r_squared))
  expect_no_error(print(r))
})
