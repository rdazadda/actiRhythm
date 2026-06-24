# Consolidated multi-episode rest-period detection (Roenneberg/MASDA).

test_that("rest.periods detects multiple consolidated bouts including a nap", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 6 * 1440)
  h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
  doy <- as.integer(format(ts, "%j"))
  counts <- ifelse(h >= 23 | h < 7, 5, 300)
  counts[h >= 14 & h < 15.5 & doy == (min(doy) + 2L)] <- 5   # a nap on day 3

  rp <- rest.periods(counts, ts)
  expect_s3_class(rp, "actiRhythm_roenneberg")
  expect_gte(rp$n_bouts, 3)                          # multiple bouts (the defining property)
  expect_true(any(rp$rest_periods$type == "nap"))    # the daytime nap is detected
  nap <- rp$rest_periods[rp$rest_periods$type == "nap", ][1, ]
  nap_h <- as.numeric(format(nap$onset, "%H"))
  expect_true(nap_h >= 12 && nap_h <= 16)
})

test_that("the closed-form correlation matches a brute-force cor()", {
  set.seed(1)
  s <- as.numeric(stats::rbinom(40, 1, 0.5))
  cf <- actiRhythm:::.rp_corr(s)
  bf <- vapply(seq_len(length(s) - 1L), function(m) {
    suppressWarnings(stats::cor(s, c(rep(1, m), rep(0, length(s) - m))))
  }, numeric(1))
  ok <- is.finite(cf) & is.finite(bf)
  expect_lt(max(abs(cf[ok] - bf[ok])), 1e-9)
})

test_that("rest.periods returns more than one bout per day on polyphasic input", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  counts <- ifelse(h < 6 | (h >= 12 & h < 14), 5, 300)   # night rest + midday rest
  rp <- rest.periods(counts, ts)
  expect_gt(rp$mean_bouts_per_day, 1)
})

test_that("rest.periods never errors on flat or short data", {
  ts1 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 300)
  r1 <- rest.periods(stats::rnorm(300, 100, 10), ts1)
  expect_equal(r1$n_bouts, 0L)

  ts2 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  r2 <- rest.periods(rep(100, 3 * 1440), ts2)
  expect_equal(r2$n_bouts, 0L)

  expect_no_error(print(r1))
  expect_no_error(print(r2))
  expect_error(rest.periods(1:10, seq_len(9)), "same length")
})

test_that("rest.periods table has the documented columns", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  rp <- rest.periods(ifelse(h >= 23 | h < 7, 5, 300), ts)
  expect_true(all(c("bout", "onset", "offset", "onset_index", "offset_index",
                    "duration_min", "date", "is_main", "type") %in%
                  names(rp$rest_periods)))
})
