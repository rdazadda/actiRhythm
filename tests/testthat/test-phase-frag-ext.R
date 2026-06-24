# Phase / shape extensions: multicomponent cosinor, AonT/AoffT, phase tests,
# rest-activity bout fragmentation.

test_that("cosinor.multicomponent recovers a known two-harmonic shape", {
  set.seed(1)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  y <- 100 + 40 * cos(2 * pi * th / 24) + 20 * cos(2 * pi * 2 * th / 24) +
       stats::rnorm(length(th), 0, 4)
  m <- cosinor.multicomponent(y, ts, criterion = "BIC")
  expect_s3_class(m, "actiRhythm_multicosinor")
  expect_equal(m$n_harmonics, 2L)                     # BIC selects the true 2-harmonic model
  expect_gt(m$r_squared, 0.95)
  expect_true(abs(m$harmonics$amplitude[1] - 40) < 6) # about 40 fundamental
  expect_true(abs(m$harmonics$amplitude[2] - 20) < 6) # about 20 second harmonic
})

test_that("activity.onset.offset finds onset at the rise and offset at the fall", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h <- as.numeric(format(ts, "%H"))
  ao <- activity.onset.offset(ifelse(h >= 23 | h < 7, 5, 300), ts)
  expect_s3_class(ao, "actiRhythm_aont")
  expect_lt(abs(ao$onset_h - 7), 2)
  expect_lt(min(abs(ao$offset_h - c(23, -1))), 2)     # offset near 23:00
})

test_that("phase.concentration flags clustered onsets and accepts scattered ones", {
  set.seed(2)
  clustered <- (23 + stats::rnorm(12, 0, 0.4)) %% 24
  pc <- phase.concentration(clustered, n_perm = 500)
  expect_s3_class(pc, "actiRhythm_phasetest")
  expect_gt(pc$R, 0.8)
  expect_lt(pc$rayleigh_p, 0.05)
  expect_lt(pc$hr_p, 0.05)                            # Hermans-Rasson also rejects clustering
  scattered <- seq(0, 24, length.out = 13)[1:12]      # evenly spread
  expect_gt(phase.concentration(scattered, n_perm = 500)$rayleigh_p, 0.05)
})

test_that("rest.activity.fragmentation summarises bout lengths and transitions", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  h <- as.numeric(format(ts, "%H"))
  f <- rest.activity.fragmentation(h >= 7 & h < 23, ts)
  expect_s3_class(f, "actiRhythm_rafrag")
  expect_true(f$mean_active_bout > 800 && f$mean_active_bout < 1000)   # about 16 h
  expect_true(f$median_rest_bout > 400 && f$median_rest_bout < 520)    # about 8 h (median robust to edge partials)
  expect_true(f$transitions_per_day > 1.5 && f$transitions_per_day < 2.5)
})

test_that("phase/frag extensions validate inputs and never error", {
  expect_error(cosinor.multicomponent(1:10, 1:9), "same length")
  expect_error(rest.activity.fragmentation(rep(TRUE, 5), 1:4), "same length")
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 50)
  expect_no_error(cosinor.multicomponent(rep(1, 50), ts))
  expect_no_error(activity.onset.offset(rep(0, 50), ts))
  expect_no_error(phase.concentration(c(1, 2)))
})
