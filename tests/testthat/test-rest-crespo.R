# Crespo rank-order + morphology rest/activity detection.

test_that("rest.crespo detects multiple consolidated rest bouts", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  counts <- ifelse(h >= 23 | h < 7, 5, 300)
  rp <- rest.crespo(counts, ts)
  expect_s3_class(rp, "actiRhythm_crespo")
  expect_gte(rp$n_rest_periods, 2)                       # several nights
  expect_gt(stats::median(rp$rest_periods$duration_min), 300)   # about 8 h nightly rest
})

test_that("rest.crespo consolidates a short daytime rest into the main period", {
  ts <- seq(as.POSIXct("2024-01-01 08:00", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  doy <- as.integer(format(ts, "%j"))
  counts <- ifelse(h >= 23 | h < 7, 5, 300)
  counts[h >= 13 & h < 16 & doy == sort(unique(doy))[2]] <- 5   # a 3 h daytime rest
  rp <- rest.crespo(counts, ts)
  # Crespo detects MAIN rest periods (Eq 4 ~8 h smoothing), so the short daytime
  # rest is suppressed rather than reported as its own bout.
  expect_lte(rp$mean_bouts_per_day, 1.25)
  expect_true(all(rp$rest_periods$duration_min > 180))
})

test_that("rest.crespo endpoints are active and the table is well formed", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  rp <- rest.crespo(ifelse(h >= 23 | h < 7, 5, 300), ts)
  expect_equal(rp$rest_state[1], "A")
  expect_equal(rp$rest_state[length(rp$rest_state)], "A")
  expect_true(all(c("bout", "onset", "offset", "onset_index", "offset_index",
                    "duration_min", "date") %in% names(rp$rest_periods)))
  expect_true(all(rp$rest_periods$offset_index >= rp$rest_periods$onset_index))
})

test_that("rest.crespo never errors on flat or short data", {
  ts1 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 300)
  r1 <- rest.crespo(stats::rnorm(300, 100, 10), ts1)
  expect_equal(r1$n_rest_periods, 0L)

  ts2 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  r2 <- rest.crespo(rep(100, 3 * 1440), ts2)
  expect_equal(r2$n_rest_periods, 0L)

  expect_no_error(print(r1))
  expect_no_error(print(r2))
  expect_error(rest.crespo(1:10, seq_len(9)), "same length")
})

test_that("rest.crespo runs on the bundled recording", {
  agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
  cr <- rest.crespo(agd$axis1, agd$timestamp)
  expect_gt(cr$n_rest_periods, 0)
  expect_true(cr$total_rest_min > 0)
})
