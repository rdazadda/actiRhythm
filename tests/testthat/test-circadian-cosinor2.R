# Cross-validation of the single-component cosinor against cosinor::cosinor.lm
# (the fitting engine behind the cosinor2 package). Both estimate
# Y = MESOR + amplitude * cos(2*pi*t/period - acrophase).

test_that("single cosinor MESOR / amplitude / acrophase match cosinor.lm and ground truth", {
  skip_if_not_installed("cosinor")

  epl <- 60
  n <- 7 * 1440
  t0 <- as.POSIXct("2024-01-06 00:00:00", tz = "UTC")
  ts <- t0 + (seq_len(n) - 1) * epl
  th <- as.numeric(difftime(ts, t0, units = "hours"))
  set.seed(1)
  # Ground truth: MESOR 100, amplitude 80, acrophase (peak) at 14:00.
  act <- 100 + 80 * cos(2 * pi * (th - 14) / 24) + rnorm(n, 0, 3)

  ca <- cosinor.analysis(act, ts)

  clm <- cosinor::cosinor.lm(Y ~ time(time), period = 24,
                             data = data.frame(time = th %% 24, Y = act))
  tt <- summary(clm)$transformed.table

  # Agreement with cosinor.lm (actiRhythm fits the averaged 24 h profile, so the
  # amplitude differs by a fraction of a percent from the raw-series fit).
  expect_equal(as.numeric(ca$mesor),     as.numeric(tt["(Intercept)", "estimate"]), tolerance = 0.02)
  expect_equal(as.numeric(ca$amplitude), as.numeric(tt["amp", "estimate"]),         tolerance = 0.02)

  # Recovery of the known parameters.
  expect_equal(as.numeric(ca$mesor),     100, tolerance = 0.01)
  expect_equal(as.numeric(ca$amplitude), 80,  tolerance = 0.01)
  expect_lt(abs(as.numeric(ca$acrophase) - 14), 0.2)
})

test_that("cosinor recovers a shifted acrophase", {
  skip_if_not_installed("cosinor")

  epl <- 60
  n <- 7 * 1440
  t0 <- as.POSIXct("2024-01-06 00:00:00", tz = "UTC")
  ts <- t0 + (seq_len(n) - 1) * epl
  th <- as.numeric(difftime(ts, t0, units = "hours"))
  set.seed(7)
  act <- 50 + 40 * cos(2 * pi * (th - 6) / 24) + rnorm(n, 0, 3)

  ca <- cosinor.analysis(act, ts)
  expect_equal(as.numeric(ca$mesor),     50, tolerance = 0.02)
  expect_equal(as.numeric(ca$amplitude), 40, tolerance = 0.02)
  expect_lt(abs(as.numeric(ca$acrophase) - 6), 0.2)
})

# Cross-validate the Bingham (1982) population-mean cosinor and its confidence
# intervals against cosinor2::population.cosinor.lm. cosinor2 takes one row per
# subject; its time reference is set to our hourly bin centres (h + 0.5), and its
# acrophase (negative radians) is converted to clock hours.
test_that("population.cosinor matches cosinor2 population mean and CIs", {
  skip_if_not_installed("cosinor2")

  set.seed(11); K <- 8L; hrs <- 0:23
  mat <- t(vapply(seq_len(K), function(i) {
    Mi <- 100 + stats::rnorm(1, 0, 5); Ai <- 40 + stats::rnorm(1, 0, 5)
    pii <- 8 + stats::rnorm(1, 0, 1)
    Mi + Ai * cos(2 * pi * (hrs - pii) / 24) + stats::rnorm(24, 0, 3)
  }, numeric(24)))

  invisible(utils::capture.output(
    pc <- cosinor2::population.cosinor.lm(data = as.data.frame(mat), time = hrs + 0.5,
                                          period = 24, plot = FALSE)))
  rad2h <- function(r) ((-r) * 24 / (2 * pi)) %% 24

  base <- as.POSIXct("2024-01-01", tz = "UTC")
  ts   <- rep(base + hrs * 3600, times = K)
  ours <- population.cosinor(as.vector(t(mat)), ts, rep(paste0("S", seq_len(K)), each = 24))

  expect_equal(ours$mesor,     pc$coefficients[["MESOR"]],            tolerance = 1e-4)
  expect_equal(ours$amplitude, pc$coefficients[["Amplitude"]],        tolerance = 1e-4)
  expect_equal(ours$acrophase, rad2h(pc$coefficients[["Acrophase"]]), tolerance = 1e-3)

  expect_equal(sort(ours$ci_mesor),     sort(as.numeric(pc$conf.ints[, "MESOR"])),     tolerance = 1e-3)
  expect_equal(sort(ours$ci_amplitude), sort(as.numeric(pc$conf.ints[, "Amplitude"])), tolerance = 1e-3)
  expect_equal(sort(ours$ci_acrophase), sort(rad2h(as.numeric(pc$conf.ints[, "Acrophase"]))), tolerance = 1e-3)
})
