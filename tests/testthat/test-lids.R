# LIDS ultradian sleep oscillation (Winnebeck 2018).

test_that("recovers a known 90-min ultradian period", {
  set.seed(1)
  base <- as.POSIXct("2024-01-01 23:00:00", tz = "UTC")
  n  <- 480L                                   # 8 h of 1-min epochs
  ts <- base + (0:(n - 1L)) * 60
  tt <- 0:(n - 1L)
  activity <- pmax(0, 5 + 4 * cos(2 * pi * tt / 90) + stats::rnorm(n, 0, 1))
  sp <- data.frame(in_bed_time = base, out_bed_time = base + n * 60)

  r <- lids(activity, ts, sp)
  expect_s3_class(r, "actiRhythm_lids")
  expect_equal(r$mean_period_min, 90, tolerance = 10)
  expect_gt(r$mean_MRI, 0)
  expect_equal(r$n_periods, 1L)
})

test_that("a strong ultradian rhythm yields a higher MRI than noise", {
  set.seed(2)
  base <- as.POSIXct("2024-01-01 23:00:00", tz = "UTC")
  n  <- 480L
  ts <- base + (0:(n - 1L)) * 60
  tt <- 0:(n - 1L)
  strong <- pmax(0, 5 + 5 * cos(2 * pi * tt / 90))
  noisy  <- pmax(0, 5 + stats::rnorm(n, 0, 5))
  sp <- data.frame(in_bed_time = base, out_bed_time = base + n * 60)

  rs <- lids(strong, ts, sp)
  rn <- lids(noisy,  ts, sp)
  expect_gt(rs$mean_MRI, rn$mean_MRI)
})

test_that("handles multiple sleep periods", {
  set.seed(3)
  base <- as.POSIXct("2024-01-01 23:00:00", tz = "UTC")
  mk <- function(start_h) {
    st <- base + start_h * 3600
    tt <- 0:479
    list(ts = st + tt * 60,
         act = pmax(0, 5 + 4 * cos(2 * pi * tt / 90) + stats::rnorm(480, 0, 1)),
         in_bed = st, out_bed = st + 480 * 60)
  }
  a <- mk(0); b <- mk(24)
  ts  <- c(a$ts, b$ts)
  act <- c(a$act, b$act)
  sp  <- data.frame(in_bed_time = c(a$in_bed, b$in_bed),
                    out_bed_time = c(a$out_bed, b$out_bed))
  r <- lids(act, ts, sp)
  expect_equal(r$n_periods, 2L)
  expect_equal(nrow(r$periods), 2L)
})
