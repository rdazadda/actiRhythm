# Curve registration and residual circadian spectrum.

test_that("curve.registration aligns days and reports a tight phase when timing is stable", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
  h <- as.numeric(format(ts, "%H"))
  cr <- curve.registration(ifelse(h >= 23 | h < 7, 5, 300), ts)
  expect_s3_class(cr, "actiRhythm_registration")
  expect_false(cr$insufficient)
  expect_gte(cr$n_days, 4L)
  expect_lt(cr$phase_sd, 1.5)                         # consistent day-to-day timing
  expect_true(cr$mean_M10 >= 10 && cr$mean_M10 <= 18) # active midpoint near midday
})

test_that("residual.spectrum finds ultradian power left after removing the 24h cosine", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  y <- 100 + 80 * cos(2 * pi * th / 24) + 30 * sin(2 * pi * th / 4)
  rcs <- residual.spectrum(y, ts, period = 24)
  expect_s3_class(rcs, "actiRhythm_rcs")
  expect_false(rcs$insufficient)
  ultra <- rcs$bands$fraction[rcs$bands$band == "ultradian"]
  hf <- rcs$bands$fraction[rcs$bands$band == "high_freq"]
  expect_gt(ultra, hf)                                # the 4h component lands in the ultradian band
})

test_that("registration / RCS validate inputs and never error on degenerate data", {
  expect_error(curve.registration(1:10, 1:9), "same length")
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 200)
  expect_true(curve.registration(rep(0, 200), ts)$insufficient)
  expect_no_error(residual.spectrum(rep(0, 200), ts))
})
