# Spectrogram, consensus rhythmicity, and batch analysis.

test_that("consensus.rhythmicity flags a strong rhythm and combines methods", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))

  r <- consensus.rhythmicity(counts, ts)
  expect_s3_class(r, "actiRhythm_consensus")
  expect_equal(nrow(r$tests), 4L)
  expect_true(r$consensus_rhythmic)
  expect_gte(r$votes, 3)
  expect_lt(r$consensus_p, 0.05)
})

test_that("consensus.rhythmicity does not flag a flat signal", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  r <- consensus.rhythmicity(rep(100, length(ts)), ts)
  expect_false(isTRUE(r$consensus_rhythmic))
})

test_that("circadian.spectrogram produces a period-by-time map", {
  set.seed(1)
  t_hours <- seq(0, 6 * 24 - 1 / 60, by = 1 / 60)
  ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
  counts <- 100 + 80 * cos(2 * pi * t_hours / 24) + stats::rnorm(length(t_hours), 0, 20)

  sp <- circadian.spectrogram(counts, ts, window_hours = 72, step_hours = 24)
  expect_s3_class(sp, "actiRhythm_spectrogram")
  expect_gt(sp$n_windows, 1)
  expect_s3_class(sp$plot, "ggplot")
  expect_true(all(c("center_time", "period_h", "power") %in% names(sp$data)))
})

test_that("circadian.spectrogram returns insufficient when shorter than the window", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440)  # 1 day
  sg <- circadian.spectrogram(stats::rnorm(1440), ts, window_hours = 72)
  expect_true(isTRUE(sg$insufficient))
  expect_s3_class(sg, "actiRhythm_spectrogram")
  expect_s3_class(sg$plot, "ggplot")
})

test_that("circadian.batch analyses every .agd in a directory", {
  dir <- system.file("extdata", package = "actiRhythm")
  b <- circadian.batch(dir, verbose = FALSE)

  expect_true(is.data.frame(b))
  expect_gte(nrow(b), 2)                       # both bundled .agd files
  expect_true(all(c("file", "error", "IS", "IV", "RA", "rhythm_p_value", "kRA") %in% names(b)))
  expect_true(all(is.finite(b$IS)))            # core metric computed per file
  expect_true(all(is.na(b$error)))             # no read/analysis errors
})

test_that("circadian.batch writes a combined workbook", {
  skip_if_not_installed("openxlsx")
  dir <- system.file("extdata", package = "actiRhythm")
  f <- tempfile(fileext = ".xlsx"); on.exit(unlink(f), add = TRUE)
  circadian.batch(dir, file = f, verbose = FALSE)
  expect_true(file.exists(f))
  expect_true("Summary" %in% openxlsx::getSheetNames(f))
})
