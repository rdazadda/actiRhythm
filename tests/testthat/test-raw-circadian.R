# Raw-native circadian: ENMO/MAD/anglez, van Hees calibration + z-angle sleep,
# circadian.raw bridge, ABI, transition probability.

test_that("auto.calibrate recovers a known gain/offset miscalibration", {
  set.seed(1); fs <- 30; win <- fs * 10; nd <- 60
  dirs <- matrix(stats::rnorm(nd * 3), nd, 3); dirs <- dirs / sqrt(rowSums(dirs^2))
  st <- c(1.03, 0.97, 1.01); of <- c(0.04, -0.03, 0.02)
  raw <- do.call(rbind, lapply(seq_len(nd), function(i) {
    base <- dirs[i, ] / st + of
    matrix(rep(base, each = win), win, 3) + matrix(stats::rnorm(win * 3, 0, 0.004), win, 3)
  }))
  colnames(raw) <- c("x", "y", "z")
  cal <- auto.calibrate(as.data.frame(raw), fs)
  expect_true(cal$calibrated)
  expect_lt(cal$cal_error_end, cal$cal_error_start)
  expect_equal(cal$scale, st, tolerance = 0.01)
  expect_equal(cal$offset, of, tolerance = 0.01)
})

test_that("auto.calibrate returns identity on too-little data", {
  cal <- auto.calibrate(data.frame(x = rep(0, 50), y = rep(0, 50), z = rep(1, 50)), 30)
  expect_false(cal$calibrated)
  expect_equal(cal$scale, c(1, 1, 1)); expect_equal(cal$offset, c(0, 0, 0))
})

test_that("auto.calibrate is robust to idle-sleep all-zero windows", {
  # real .gt3x files carry imputed all-zero (idle-sleep) rows; those windows are
  # off the gravity sphere and must not produce NaN in the fit
  set.seed(1); fs <- 30; win <- fs * 10; nd <- 60
  dirs <- matrix(stats::rnorm(nd * 3), nd, 3); dirs <- dirs / sqrt(rowSums(dirs^2))
  good <- do.call(rbind, lapply(seq_len(nd), function(i)
    matrix(rep(dirs[i, ], each = win), win, 3) + stats::rnorm(win * 3, 0, 0.004)))
  raw <- rbind(good, matrix(0, win * 20, 3))           # + 20 all-zero windows
  colnames(raw) <- c("x", "y", "z")
  cal <- auto.calibrate(as.data.frame(raw), fs)
  expect_true(cal$calibrated)
  expect_true(all(is.finite(cal$scale)) && all(is.finite(cal$offset)))
})

test_that(".raw_metrics_compute gives ENMO near 0 and anglez near 90 for flat 1g on z", {
  flat <- cbind(rep(0, 1500), rep(0, 1500), rep(1, 1500))
  m <- actiRhythm:::.raw_metrics_compute(flat, 30, 10, c("ENMO", "MAD", "anglez"),
                                         as.POSIXct("2024-01-01", tz = "UTC"), "UTC")
  expect_lt(mean(m$ENMO), 1e-6)
  expect_equal(mean(m$anglez), 90, tolerance = 1e-6)
})

test_that("the z-angle sleep pipeline finds nights and scores them", {
  ep <- 5; pd <- 24 * 3600 / ep
  ts <- seq(as.POSIXct("2024-01-01 12:00", tz = "UTC"), by = ep, length.out = 2 * pd)
  hh <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  night <- hh >= 23 | hh < 7
  set.seed(2); az <- numeric(length(ts))
  az[!night] <- -30 + 40 * sin(2 * pi * seq_len(sum(!night)) / 200) + stats::rnorm(sum(!night), 0, 5)
  az[night]  <- -60 + stats::rnorm(sum(night), 0, 0.02)
  spt <- rest.spt(az, ts, epoch_length = 5)
  sib <- sib.vanhees(az, epoch_length = 5)
  slp <- sleep.from.spt(spt, sib, ts, epoch_length = 5)
  expect_s3_class(spt, "actiRhythm_spt")
  expect_equal(nrow(spt), 2L)
  expect_true(all(abs(spt$duration - 8) < 1))          # about 8h nights
  expect_gt(mean(sib[night] == "S"), 0.95)             # still nights scored sleep
  expect_equal(nrow(slp), 2L)
  expect_true(all(abs(slp$tst - 8) < 1))
  # HorAngle variant runs
  expect_s3_class(rest.spt(az, ts, epoch_length = 5, algo = "HorAngle"), "actiRhythm_spt")
})

test_that("rest.spt and sleep.from.spt are empty-safe on short input", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 5, length.out = 100)
  spt <- rest.spt(rep(0, 100), ts, epoch_length = 5)
  expect_equal(nrow(spt), 0L)
  expect_equal(nrow(sleep.from.spt(spt, rep("S", 100), ts)), 0L)
})

test_that("activity.balance.index peaks at alpha = 1", {
  expect_equal(activity.balance.index(1.0), 1)
  expect_lt(activity.balance.index(0.5), 1)
  expect_lt(activity.balance.index(1.5), 1)
  abi <- activity.balance.index(list(alpha = 1, alpha1 = 0.8, alpha2 = 1.2))
  expect_equal(abi$ABI_overall, 1)
  expect_true(is.finite(abi$ABI_short) && is.finite(abi$ABI_long))
})

test_that("transition.probability counts bouts correctly", {
  cnt <- c(rep(0, 50), rep(100, 20), rep(0, 40), rep(80, 30), rep(0, 60))
  tp <- transition.probability(cnt)
  expect_equal(tp$n_active_bouts, 2L)
  expect_equal(tp$mean_active_bout, 25)
  expect_equal(tp$tp_ar_mle, 2 / 50)                   # active at-risk = 50 (trailing bout is rest)
  expect_equal(tp$tp_ra_mle, 2 / 149)                  # rest at-risk = 150 - 1 (trailing rest bout)
  expect_true(tp$tp_ar_bayes > 0 && tp$tp_ar_bayes < 1)
})

test_that("example_raw + data-frame input run end to end", {
  raw <- example_raw(days = 1)
  expect_s3_class(raw, "data.frame")
  expect_true(all(c("time", "x", "y", "z") %in% names(raw)))
  expect_equal(attr(raw, "fs"), 30)
  expect_equal(nrow(raw), 1L * 86400L * 30L)
  m <- raw.metrics(raw, epoch = 60)                     # data-frame (non-file) input
  expect_true(all(c("time", "ENMO", "MAD", "anglez") %in% names(m)))
  expect_equal(nrow(m), 1440L)                          # 24h at 60s
  hh <- as.numeric(format(m$time, "%H"))
  expect_gt(mean(m$ENMO[hh >= 10 & hh < 18]), mean(m$ENMO[hh < 6]))   # day/night rhythm
  bad <- data.frame(a = 1:10); attr(bad, "fs") <- 30
  expect_error(raw.metrics(bad), "column")                           # missing x/y/z
})

test_that("detect.nonwear.raw + the rest.spt wear gate exclude a device-off period", {
  raw  <- example_raw(days = 2, device_off = 1, fs = 15)   # 2 worn days + 1 off day
  wear <- detect.nonwear.raw(raw, epoch = 5)
  expect_gt(mean(wear), 0.6); expect_lt(mean(wear), 0.72)  # about 2/3 worn
  m5 <- raw.metrics(raw, epoch = 5, metrics = "anglez")
  spt0 <- rest.spt(m5$anglez, m5$time, epoch_length = 5)               # ungated
  sptW <- rest.spt(m5$anglez, m5$time, epoch_length = 5, wear = wear)  # gated
  expect_gt(max(spt0$duration), 20)                        # device-off -> 24h "sleep"
  expect_equal(nrow(sptW), 2L)                             # gate removes it
  expect_true(all(sptW$duration < 12))                     # only the real nights remain
  # a mismatched mask is ignored, not an error
  expect_s3_class(rest.spt(m5$anglez, m5$time, wear = c(TRUE, FALSE)), "actiRhythm_spt")
})

test_that("raw.metrics and circadian.raw work on a real .gt3x", {
  skip_if_not_installed("read.gt3x"); skip_if_not_installed("agcounts")
  g <- system.file("extdata", "TAS1H30182785_2019-09-17.gt3x", package = "read.gt3x")
  skip_if(!nzchar(g))
  rm <- raw.metrics(g, epoch = 60)
  expect_true(all(c("time", "ENMO", "MAD", "anglez") %in% names(rm)))
  expect_true(all(rm$ENMO >= 0))
  expect_true(all(rm$anglez >= -90 & rm$anglez <= 90))
  cr <- circadian.raw(g, metric = "ENMO", epoch = 60, min_valid_hours = 0)
  expect_s3_class(cr, "actiRhythm_circadian")
})
