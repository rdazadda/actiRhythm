# Rest-activity state transitions, validated against a known Markov process.

markov_states <- function(n, p, q, seed = 1) {
  # p = P(rest -> active), q = P(active -> rest); geometric bouts -> constant hazard.
  set.seed(seed)
  u <- stats::runif(n)
  s <- integer(n)
  for (i in 2:n) {
    s[i] <- if (s[i - 1L] == 0L) (if (u[i] < p) 1L else 0L)
            else                 (if (u[i] < q) 0L else 1L)
  }
  s
}

test_that("recovers known Markov transition rates", {
  s <- markov_states(120000, p = 0.05, q = 0.10)
  r <- state.transitions(s * 100, threshold = 1)

  expect_s3_class(r, "actiRhythm_transitions")
  expect_equal(r$pRA, 0.05, tolerance = 0.05)   # 1 / mean rest bout = p
  expect_equal(r$pAR, 0.10, tolerance = 0.05)
  expect_equal(r$kRA, 0.05, tolerance = 0.10)   # sustained hazard near p (geometric)
  expect_equal(r$kAR, 0.10, tolerance = 0.10)
})

test_that("more fragmented rest gives a higher kRA", {
  calm <- state.transitions(markov_states(60000, p = 0.02, q = 0.10, seed = 2) * 100)
  frag <- state.transitions(markov_states(60000, p = 0.12, q = 0.10, seed = 3) * 100)
  expect_gt(frag$kRA, calm$kRA)
})

test_that("threshold controls the rest/active split", {
  counts <- rep(c(0, 5, 50), length.out = 6000)
  lo <- state.transitions(counts, threshold = 1)    # 5 and 50 are active (longer active bouts)
  hi <- state.transitions(counts, threshold = 10)   # only 50 is active (shorter active bouts)
  expect_gt(hi$pAR, lo$pAR)                          # shorter active bouts -> higher pAR
})

test_that("degrades gracefully on too-short input", {
  r <- state.transitions(c(0, 1))
  expect_true(isTRUE(r$insufficient))
  expect_s3_class(r, "actiRhythm_transitions")
  expect_true(is.na(r$kRA))
})
