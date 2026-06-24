# Tier A/B additions: Bingham joint test, wavelet COI, non-wear, CSV ingest, scorer NA.

.two_group <- function(acroA = 8, acroB = 11, n = 5, seed = 1) {
  set.seed(seed); hrs <- 0:23; act <- ts <- subj <- grp <- NULL
  for (g in c("A", "B")) for (i in seq_len(n)) {
    acro <- if (g == "A") acroA else acroB
    y <- 100 + 40 * cos(2 * pi * (hrs - acro) / 24) + stats::rnorm(24, 0, 4)
    act <- c(act, y); subj <- c(subj, rep(paste0(g, i), 24)); grp <- c(grp, rep(g, 24))
    ts <- c(ts, as.POSIXct("2024-01-01", tz = "UTC") + hrs * 3600)
  }
  list(act = act, ts = as.POSIXct(ts, tz = "UTC", origin = "1970-01-01"),
       subj = subj, grp = grp)
}

test_that("cosinor.compare reports a valid Bingham joint test that rejects a real difference", {
  d <- .two_group(8, 11)
  cc <- cosinor.compare(d$act, d$ts, d$subj, d$grp)
  expect_true(cc$joint$valid)
  expect_equal(cc$joint$df1, 3L); expect_equal(cc$joint$df2, 6L)   # K=10 -> K-4
  expect_lt(cc$joint$p_value, 0.05)                                # groups differ
  # identical groups -> joint not significant
  s <- .two_group(8, 8, seed = 7)
  expect_gt(cosinor.compare(s$act, s$ts, s$subj, s$grp)$joint$p_value, 0.05)
})

test_that("circadian.wavelet adds a Torrence-Compo cone of influence", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  w <- circadian.wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
  expect_equal(w$coi_multiplier, 0.7305, tolerance = 1e-3)         # ff/sqrt(2) at omega0=6
  expect_equal(length(w$coi_period_h), length(ts))
  expect_equal(which.max(w$coi_period_h), length(ts) %/% 2L, tolerance = 2)  # peaks mid-series
  expect_lt(w$coi_period_h[1], w$coi_period_h[length(ts) %/% 2L])  # rises from the edge
  expect_lt(abs(w$peak_period_h - 24), 4)                          # still recovers 24h
})

test_that("wavelet.coi matches the standalone formula", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 200)
  coi <- wavelet.coi(ts, epoch_length = 600)
  expect_equal(length(coi), 200L)
  expect_true(coi[1] < 0.01 && coi[200] < 0.01)                   # sliver at both ends
})

test_that("non-wear detectors recover an injected gap and differ as designed", {
  counts <- c(rep(200, 200), rep(0, 200), rep(200, 200))
  choi <- detect.nonwear.choi(counts); tro <- detect.nonwear.troiano(counts)
  expect_equal(sum(!choi), 200L); expect_equal(sum(!tro), 200L)   # the 200-min zero block
  expect_true(all(choi[1:200]) && all(choi[401:600]))             # worn stretches kept
  # Troiano's 100-count stop level: a >100 spike splits the gap below the frame
  hi <- c(rep(0, 30), 500, rep(0, 30))                            # 500 spike -> two 30-min halves
  lo <- c(rep(0, 30), 50, rep(0, 30))                             # 50 spike (<= stoplevel) absorbed
  expect_equal(sum(!detect.nonwear.troiano(hi)), 0L)             # high spike not tolerated
  expect_gt(sum(!detect.nonwear.troiano(lo)), 60)               # low spike absorbed -> non-wear
  # Choi tolerates a small flanked spike inside a long gap
  gc <- c(rep(0, 60), 50, rep(0, 60))                            # 1-min spike flanked by zeros
  expect_gt(sum(!detect.nonwear.choi(gc, frame = 100)), 100)     # still one non-wear bout
})

test_that("read.actigraph.csv parses an ActiLife epoch CSV", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "------------ Data File Created By ActiGraph GT3X+ ActiLife v6.13.3 date format M/d/yyyy at 30 Hz Filter Normal -----------",
    "Serial Number: NEO1B41100262", "Start Time 10:00:00", "Start Date 1/1/2024",
    "Epoch Period (hh:mm:ss) 00:01:00", "Download Time 09:00:00", "Download Date 1/8/2024",
    "Current Memory Address: 0", "Current Battery Voltage: 4.19 Mode = 12",
    "------------ Data Table File Created By ActiGraph -----------",
    "Axis1,Axis2,Axis3", "100,50,20", "200,80,40", "0,0,0"), f)
  d <- read.actigraph.csv(f)
  expect_equal(nrow(d), 3L)
  expect_true(all(c("timestamp", "axis1", "axis2", "axis3") %in% names(d)))
  expect_s3_class(d$timestamp, "POSIXct")
  expect_equal(format(d$timestamp[1], "%Y-%m-%d %H:%M:%S"), "2024-01-01 10:00:00")
  expect_equal(as.numeric(difftime(d$timestamp[2], d$timestamp[1], units = "secs")), 60)
  expect_equal(d$axis1, c(100, 200, 0))
})

test_that("counts.from.data.frame returns timestamp + counts", {
  df <- data.frame(a = c(1, 2, 3),
                   t = c("2024-01-01 00:00:00", "2024-01-01 00:01:00", "2024-01-01 00:02:00"))
  out <- counts.from.data.frame(df, "a", "t")
  expect_s3_class(out$timestamp, "POSIXct")
  expect_equal(out$counts, c(1, 2, 3))
  expect_error(counts.from.data.frame(df, "missing", "t"), "not found")
  expect_warning(counts.from.data.frame(df, "a"), "synthesized")
})

test_that("sleep scorers honour na_action", {
  x <- c(0, NA, 400, 0, NA)
  expect_equal(sleep.cole.kripke(x), c("S", NA, "S", "S", NA))           # default na
  expect_equal(sleep.sadeh(x, na_action = "wake")[c(2, 5)], c("W", "W")) # NA -> wake
  expect_false(anyNA(sleep.cole.kripke(x, na_action = "zero")))          # NA -> scored
  # the NA states pass through SRI as missing without breaking it
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
  h <- as.numeric(format(ts, "%H"))
  cnt <- ifelse(h >= 23 | h < 7, 2, 250); cnt[1000:1100] <- NA
  st <- sleep.cole.kripke(cnt)
  expect_true(anyNA(st) && is.finite(sleep.regularity.index(st, ts)))
})
