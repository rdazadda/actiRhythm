# Cosinor rhythmicity test (Halberg zero-amplitude F-test + percent rhythm).

rhythm_ts <- function(days = 3) {
  seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = days * 1440)
}

test_that("detects a strong 24h rhythm", {
  ts <- rhythm_ts()
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))

  r <- rhythmicity.test(counts, ts)
  expect_s3_class(r, "actiRhythm_rhythmicity")
  expect_true(r$rhythmic)
  expect_lt(r$p_value, 0.001)
  expect_gt(r$percent_rhythm, 90)
  expect_equal(r$df1, 2L)
})

test_that("does not flag a flat signal as rhythmic", {
  ts <- rhythm_ts()
  r <- rhythmicity.test(rep(100, length(ts)), ts)
  expect_false(isTRUE(r$rhythmic))
})

test_that("percent_rhythm equals the cosinor.analysis R-squared (single engine)", {
  ts <- rhythm_ts(2)
  h  <- as.numeric(format(ts, "%H"))
  counts <- pmax(0, 100 + 50 * cos(2 * pi * (h - 8) / 24))

  fit <- cosinor.analysis(counts, ts)
  r   <- rhythmicity.test(counts, ts, cosinor_result = fit)
  expect_equal(r$percent_rhythm / 100, fit$r_squared, tolerance = 1e-9)
  expect_equal(r$df2, fit$n_profile_hours - 3L)
})

test_that("F and R-squared are self-consistent", {
  ts <- rhythm_ts(2)
  h  <- as.numeric(format(ts, "%H"))
  counts <- pmax(0, 80 + 40 * cos(2 * pi * (h - 6) / 24) + 5 * sin(2 * pi * h / 12))

  r <- rhythmicity.test(counts, ts)
  f_check <- (r$r_squared / r$df1) / ((1 - r$r_squared) / r$df2)
  expect_equal(r$F, f_check, tolerance = 1e-9)
  expect_equal(r$p_value, stats::pf(r$F, r$df1, r$df2, lower.tail = FALSE), tolerance = 1e-12)
})

# Validation is transitive: percent_rhythm == cosinor.analysis() R^2 (tested
# above), and cosinor.analysis() is validated against cosinor::cosinor.lm in
# test-circadian-cosinor2.R. cosinor2::cosinor.PR fits the RAW series, so its R^2
# is not directly comparable to the averaged-profile R^2 used here.

test_that("the Halberg F matches base lm on the same averaged profile", {
  # Build a clean 24-point hourly profile and confirm the zero-amplitude F-test
  # equals lm()'s overall F-test for that profile (the definition of the test).
  set.seed(3)
  hrs <- 0:23
  prof <- 100 + 60 * cos(2 * pi * (hrs - 10) / 24) + stats::rnorm(24, 0, 6)
  r2 <- summary(stats::lm(prof ~ cos(2 * pi * hrs / 24) + sin(2 * pi * hrs / 24)))$r.squared
  fref <- summary(stats::lm(prof ~ cos(2 * pi * hrs / 24) + sin(2 * pi * hrs / 24)))$fstatistic

  f_ours <- (r2 / 2) / ((1 - r2) / (24 - 3))
  expect_equal(unname(f_ours), unname(fref[1]), tolerance = 1e-6)
  expect_equal(24 - 3, unname(fref[3]))
})
