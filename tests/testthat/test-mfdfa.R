# Multifractal DFA, validated against fractal.dfa and known fractal exponents.

test_that("white noise is roughly monofractal with h(2) near 0.5", {
  set.seed(1)
  m <- mfdfa(stats::rnorm(4096))

  expect_s3_class(m, "actiRhythm_mfdfa")
  expect_equal(m$alpha_dfa, 0.5, tolerance = 0.12)
  hq <- m$h_q[is.finite(m$h_q)]
  expect_lt(max(hq) - min(hq), 0.3)          # roughly flat h(q) = monofractal
})

test_that("h(2) matches fractal.dfa alpha (start-only convention)", {
  set.seed(2)
  x <- cumsum(stats::rnorm(4096))            # Brownian motion, alpha near 1.5
  m <- mfdfa(x, both_ends = FALSE)
  d <- fractal.dfa(x)
  expect_equal(m$alpha_dfa, d$alpha, tolerance = 0.15)
  expect_gt(m$alpha_dfa, 1.2)                # > white noise, consistent with fBm
})

test_that("returns an NA structure on insufficient data", {
  m <- mfdfa(stats::rnorm(20))
  expect_s3_class(m, "actiRhythm_mfdfa")
  expect_true(is.na(m$alpha_dfa))
  expect_no_error(print(m))                  # print must not throw on the NA result
})

test_that("tau(q) and the Holder spectrum are consistent", {
  set.seed(3)
  m <- mfdfa(cumsum(stats::rnorm(4096)))
  # tau(q) = q*h(q) - 1 by construction
  expect_equal(m$tau_q, m$q_values * m$h_q - 1, tolerance = 1e-9)
  # f(alpha) = q*alpha - tau(q)
  expect_equal(m$f_alpha, m$q_values * m$alpha - m$tau_q, tolerance = 1e-9)
})
