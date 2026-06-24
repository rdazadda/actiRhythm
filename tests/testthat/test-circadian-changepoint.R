# Change-point sleep/wake detection.

test_that("sleep.changepoints recovers known sleep and wake onsets", {
  set.seed(1)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 6 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- ifelse(h >= 8 & h < 23, 300, 5) + pmax(0, stats::rnorm(length(ts), 0, 5))

  cp <- sleep.changepoints(counts, ts)
  expect_s3_class(cp, "actiRhythm_changepoints")
  expect_false(cp$insufficient)
  expect_gte(cp$n_episodes, 3)

  sot <- cp$changepoints$time[cp$changepoints$type == "sleep onset"]
  wot <- cp$changepoints$time[cp$changepoints$type == "wake onset"]
  sot_h <- as.numeric(format(sot, "%H")) + as.numeric(format(sot, "%M")) / 60
  wot_h <- as.numeric(format(wot, "%H")) + as.numeric(format(wot, "%M")) / 60
  expect_true(all(abs(sot_h - 23) < 1.5))   # sleep onset near 23:00
  expect_true(all(abs(wot_h - 8) < 1.5))    # wake onset near 08:00
})

test_that("sleep.changepoints sleep durations are reasonable", {
  set.seed(2)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 7 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  counts <- ifelse(h >= 7 & h < 22, 250, 8) + pmax(0, stats::rnorm(length(ts), 0, 6))
  cp <- sleep.changepoints(counts, ts)
  expect_true(abs(cp$mean_sleep_duration - 9) < 2)   # about 9 h rest (22:00-07:00)
})

test_that("sleep.changepoints never errors on short or flat data", {
  ts1 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 600)
  r1 <- sleep.changepoints(stats::rnorm(600, 100, 10), ts1)
  expect_true(r1$insufficient)

  ts2 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  r2 <- sleep.changepoints(rep(100, 4 * 1440), ts2)
  expect_true(r2$insufficient)

  expect_no_error(print(r1))
  expect_no_error(print(r2))
})

test_that("sleep.changepoints validates its inputs", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 100)
  expect_error(sleep.changepoints(1:99, ts), "same length")
  expect_error(sleep.changepoints(1:100, seq_len(100)), "POSIXct")
})
