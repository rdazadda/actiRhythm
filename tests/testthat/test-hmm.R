# State-space (Gaussian HMM) rest-activity model.

test_that("rest.hmm recovers rest and active states", {
  set.seed(1)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  h <- as.numeric(format(ts, "%H"))
  true_rest <- h >= 23 | h < 7
  counts <- ifelse(true_rest, 2, 250) + pmax(0, stats::rnorm(length(ts), 0, 10))
  m <- rest.hmm(counts, ts, n_starts = 3L)
  expect_s3_class(m, "actiRhythm_hmm")
  expect_false(m$insufficient)
  expect_equal(nrow(m$emission), 2L)
  expect_lt(m$emission$mean_transformed[1], m$emission$mean_transformed[2])  # rest < active
  agree <- mean((m$state_path == 1L) == true_rest)
  expect_gt(max(agree, 1 - agree), 0.85)                                     # decoded path matches truth
  tod <- m$tod_profile
  expect_gt(mean(tod$p_rest[tod$hour %in% 0:5]), mean(tod$p_rest[tod$hour %in% 10:14]))
})

test_that("rest.hmm yields a usable sleep_state and finite AIC", {
  set.seed(2)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  h <- as.numeric(format(ts, "%H"))
  m <- rest.hmm(ifelse(h >= 23 | h < 7, 2, 250) + pmax(0, stats::rnorm(length(ts), 0, 10)),
                ts, n_starts = 3L)
  expect_true(all(m$sleep_state %in% c("S", "W")))
  expect_true(is.finite(m$AIC) && is.finite(m$BIC))
})

test_that("rest.hmm never errors on flat or short data", {
  ts12 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 12)
  expect_true(rest.hmm(stats::rnorm(12, 100, 5), ts12)$insufficient)   # below the fit threshold
  ts2 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 300)
  expect_true(rest.hmm(rep(50, 300), ts2)$insufficient)                # flat -> sd 0
  ts3 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 40)
  expect_no_error(rest.hmm(stats::rnorm(40, 100, 5), ts3))             # short noise: no crash
  expect_error(rest.hmm(1:10, 1:9), "same length")
  expect_no_error(print(rest.hmm(rep(50, 300), ts2)))
})
