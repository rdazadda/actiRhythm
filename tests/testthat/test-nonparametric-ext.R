# Nonparametric extensions: IVm, MX/LX, Dichotomy Index, per-day metrics.

.ra <- function(days = 5) {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = days * 1440)
  h <- as.numeric(format(ts, "%H"))
  list(ts = ts, counts = ifelse(h >= 23 | h < 7, 5, 300), h = h)
}

test_that("activity.extrema puts the least-active window at night, most-active by day", {
  d <- .ra()
  ae <- activity.extrema(d$counts, d$ts, windows = c(5, 10))
  expect_s3_class(ae, "actiRhythm_mxlx")
  r5 <- ae$table[ae$table$window_h == 5, ]
  expect_lt(r5$L_mean, r5$M_mean)
  expect_true(r5$L_onset_h >= 22 || r5$L_onset_h <= 3)        # L5 onset in the night
  m10 <- ae$table[ae$table$window_h == 10, ]
  expect_true(m10$M_onset_h >= 7 && m10$M_onset_h <= 14)      # M10 onset in the day
})

test_that("dichotomy.index is high when rest is quiet relative to the active day", {
  d <- .ra(2)
  di <- dichotomy.index(d$counts, rest = d$h >= 23 | d$h < 7)
  expect_s3_class(di, "actiRhythm_dichotomy")
  expect_gt(di$IO, 90)
})

test_that("intradaily.variability.multiscale returns a finite averaged IVm", {
  set.seed(1)
  d <- .ra(3)
  iv <- intradaily.variability.multiscale(d$counts + stats::rnorm(length(d$counts), 0, 2), d$ts)
  expect_s3_class(iv, "actiRhythm_ivm")
  expect_true(is.finite(iv$IVm))
  expect_equal(nrow(iv$table), 60L)
})

test_that("circadian.daily returns one row per day with the expected columns", {
  d <- .ra(5)
  cd <- circadian.daily(d$counts, d$ts)
  expect_s3_class(cd, "actiRhythm_daily")
  expect_gte(cd$n_days, 4L)
  expect_true(all(c("L5", "M10", "RA", "IV") %in% names(cd$daily)))
  expect_true(all(cd$daily$M10 > cd$daily$L5))
})

test_that("nonparametric extensions validate inputs and never error on degenerate data", {
  expect_error(activity.extrema(1:10, 1:9), "same length")
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 200)
  expect_no_error(intradaily.variability.multiscale(rep(0, 200), ts))
  expect_no_error(circadian.daily(rep(0, 200), ts))
  expect_no_error(activity.extrema(rep(0, 200), ts))
})
