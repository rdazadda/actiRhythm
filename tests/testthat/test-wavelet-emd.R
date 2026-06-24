# Wavelet (CWT, ultradian bands) and EMD / Hilbert-Huang.

test_that("circadian.wavelet recovers a 24-hour period", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  w <- circadian.wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
  expect_s3_class(w, "actiRhythm_wavelet")
  expect_false(w$insufficient)
  expect_lt(abs(w$peak_period_h - 24), 4)
  expect_lt(abs(stats::median(w$dominant_period, na.rm = TRUE) - 24), 5)  # COI masks edge cols
})

test_that("ultradian.bandpower localises a 4-hour rhythm to the 4h band", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  bp <- ultradian.bandpower(100 + 50 * sin(2 * pi * th / 4), ts)
  expect_s3_class(bp, "actiRhythm_bandpower")
  expect_equal(as.character(bp$bands$band[which.max(bp$bands$fraction)]), "4h")
})

test_that("circadian.emd extracts a circadian IMF and reconstructs the signal", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  e <- circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
  expect_s3_class(e, "actiRhythm_emd")
  expect_false(e$insufficient)
  expect_false(is.na(e$circadian_imf))
  expect_lt(abs(e$circadian_period - 24), 6)
  expect_lt(e$recon_error, 1e-6)               # completeness
})

test_that("hilbert.huang gives an instantaneous period near 24 h", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  e <- circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
  hh <- hilbert.huang(e)
  expect_s3_class(hh, "actiRhythm_hht")
  expect_false(hh$insufficient)
  expect_lt(abs(hh$mean_period - 24), 6)
  expect_gt(hh$frac_in_band, 0.4)
})

test_that("wavelet/EMD validate inputs and never error on degenerate data", {
  expect_error(circadian.wavelet(1:10, 1:9), "same length")
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 40)
  expect_no_error(circadian.wavelet(rep(0, 40), ts))
  expect_no_error(circadian.emd(rep(0, 40), ts))
  expect_no_error(ultradian.bandpower(rep(0, 40), ts))
  expect_true(circadian.emd(rep(5, 40), ts)$insufficient)
})
