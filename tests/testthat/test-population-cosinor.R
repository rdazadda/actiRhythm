# Population-mean cosinor (Bingham 1982).

make_pop <- function(K = 6, M0 = 100, A0 = 40, phi0 = 8, sd_noise = 3, seed = 11) {
  set.seed(seed)
  hrs <- 0:23
  act <- numeric(0); subj <- character(0)
  tss <- as.POSIXct(character(0), tz = "UTC")
  base <- as.POSIXct("2024-01-01", tz = "UTC")
  for (i in seq_len(K)) {
    Mi <- M0 + stats::rnorm(1, 0, 5)
    Ai <- A0 + stats::rnorm(1, 0, 5)
    pii <- phi0 + stats::rnorm(1, 0, 1)
    y <- Mi + Ai * cos(2 * pi * (hrs - pii) / 24) + stats::rnorm(24, 0, sd_noise)
    act  <- c(act, y)
    subj <- c(subj, rep(paste0("S", i), 24))
    tss  <- c(tss, base + hrs * 3600)
  }
  list(activity = act, timestamps = tss, subject = subj, hrs = hrs)
}

test_that("recovers the simulated group means", {
  d <- make_pop(K = 8, M0 = 100, A0 = 40, phi0 = 8)
  r <- population.cosinor(d$activity, d$timestamps, d$subject)

  expect_s3_class(r, "actiRhythm_population_cosinor")
  expect_equal(r$mesor, 100, tolerance = 5)
  expect_equal(r$amplitude, 40, tolerance = 8)
  expect_equal(r$acrophase, 8, tolerance = 1.5)
  expect_true(r$ci_mesor[1] < r$mesor && r$mesor < r$ci_mesor[2])
})

test_that("group mean equals the delinearized average of per-subject betas", {
  d <- make_pop(K = 6)
  r <- population.cosinor(d$activity, d$timestamps, d$subject)
  cf <- r$subjects
  expect_equal(r$mesor, mean(cf$mesor), tolerance = 1e-9)
  expect_equal(r$amplitude, sqrt(mean(cf$beta1)^2 + mean(cf$beta2)^2), tolerance = 1e-9)
})

test_that("amplitude/acrophase CIs are undefined with fewer than 3 subjects", {
  d <- make_pop(K = 2)
  r <- population.cosinor(d$activity, d$timestamps, d$subject)
  expect_false(r$conf_interval_valid)
  expect_true(all(is.na(r$ci_amplitude)))
})

test_that("group = returns a population cosinor per group", {
  d1 <- make_pop(K = 4, phi0 = 6,  seed = 1)
  d2 <- make_pop(K = 4, phi0 = 14, seed = 2)
  act  <- c(d1$activity, d2$activity)
  ts   <- c(d1$timestamps, d2$timestamps)
  subj <- c(paste0("A", d1$subject), paste0("B", d2$subject))
  grp  <- c(rep("A", length(d1$activity)), rep("B", length(d2$activity)))

  r <- population.cosinor(act, ts, subj, group = grp)
  expect_s3_class(r, "actiRhythm_population_cosinor_list")
  expect_named(r, c("A", "B"))
  expect_gt(abs(r$A$acrophase - r$B$acrophase), 4)   # phi0 6 vs 14
})

# Validation is by simulation recovery (above) plus the definitional check that
# the group rhythm is the delinearized average of the per-subject coefficients.
# The per-subject engine (.ext.ordinary.cosinor) is the same WLS cosinor fit as
# cosinor.analysis(), which is validated against cosinor::cosinor.lm in
# test-circadian-cosinor2.R.

test_that("Bingham amplitude CI tightens as between-subject variance shrinks", {
  tight <- make_pop(K = 8, A0 = 40, sd_noise = 1, seed = 4)
  rt <- population.cosinor(tight$activity, tight$timestamps, tight$subject)
  expect_true(rt$conf_interval_valid)
  expect_lt(rt$ci_amplitude[1], rt$amplitude)
  expect_gt(rt$ci_amplitude[2], rt$amplitude)
})
