# Tier-1/2 batch: wavelet red-noise significance, error-contract harmonization,
# MSFsc/SJLsc, multicomponent-cosinor WLS profile fit.

test_that("circadian.wavelet flags AR(1)-significant power and rejects red noise", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 8 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  w  <- circadian.wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
  expect_true(is.matrix(w$significant))
  expect_equal(length(w$sig_power), length(w$period_hours))
  expect_true(is.finite(w$phi))
  fsig <- function(x) sum(x$significant) / max(1, sum(!x$in_coi))
  set.seed(1)
  wn <- circadian.wavelet(cumsum(stats::rnorm(length(ts))), ts, epoch_length = 600)
  expect_gt(fsig(w), fsig(wn))                       # rhythm has more significant power
  expect_lt(fsig(wn), 0.05)                          # red noise at/below the 5% rate
  sigrows <- which(rowSums(w$significant) > 0)       # significant cells sit at about 24h
  expect_true(any(abs(w$period_hours[sigrows] - 24) < 6))
})

test_that("spectrogram / lids / state.transitions degrade instead of erroring", {
  ts1 <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440)
  sg <- circadian.spectrogram(stats::rnorm(1440), ts1, window_hours = 72)
  expect_true(isTRUE(sg$insufficient)); expect_s3_class(sg$plot, "ggplot")
  expect_true(isTRUE(state.transitions(c(0, 1))$insufficient))
  sp <- data.frame(in_bed_time = as.POSIXct("2024-06-01 23:00", tz = "UTC"),
                   out_bed_time = as.POSIXct("2024-06-02 06:00", tz = "UTC"))
  expect_true(isTRUE(lids(rep(0, 5), ts1[1:5], sp)$insufficient))   # no in-window epochs
  expect_error(lids(1, ts1, data.frame(a = 1)), "needs")           # input-contract stop kept
})

test_that("social.jet.lag adds the MCTQ MSFsc / SJLsc correction", {
  ib <- as.POSIXct(c(paste0("2024-01-0", 1:5, " 23:00"),
                     "2024-01-06 01:00", "2024-01-07 01:00"), tz = "UTC")
  ob <- as.POSIXct(c(paste0("2024-01-0", 2:6, " 06:00"),
                     "2024-01-06 10:00", "2024-01-07 10:00"), tz = "UTC")
  sjl <- social.jet.lag(data.frame(in_bed_time = ib, out_bed_time = ob))
  expect_equal(sjl$MSF, 5.5, tolerance = 0.1)
  expect_lt(sjl$MSFsc, sjl$MSF)                      # free-day sleep longer -> corrected down
  expect_equal(sjl$social_jet_lag_sc_hours, round(sjl$MSFsc - sjl$MSW, 2), tolerance = 0.05)
})

test_that("cosinor.multicomponent fits the weighted averaged profile", {
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  set.seed(1)
  mc <- cosinor.multicomponent(100 + 50 * cos(2 * pi * th / 24) +
          20 * cos(2 * pi * 2 * th / 24) + stats::rnorm(length(ts), 0, 5), ts)
  expect_false(mc$insufficient)
  expect_gte(mc$n_harmonics, 1L)
  expect_gt(mc$r_squared, 0.9)
  expect_true(all(c("harmonic", "amplitude", "acrophase_h") %in% names(mc$harmonics)))
})
