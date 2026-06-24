# Epoch-level sleep/wake scoring (Cole-Kripke, Sadeh).

# Three days of a clear pattern: 9 h rest (zero) then 15 h active.
.rest_active <- function() {
  rep(c(rep(0, 540), rep(400, 900)), 3)
}
.rest_mask <- function() rep(c(rep(TRUE, 540), rep(FALSE, 900)), 3)

test_that("sleep.cole.kripke scores a clear rest/active pattern", {
  counts <- .rest_active()
  st <- sleep.cole.kripke(counts)
  expect_equal(length(st), length(counts))
  expect_true(all(st %in% c("S", "W")))
  rest <- .rest_mask()
  expect_gt(mean(st[rest] == "S"), 0.8)        # zero blocks -> sleep
  expect_gt(mean(st[!rest] == "W"), 0.8)       # active blocks -> wake
})

test_that("sleep.sadeh scores a clear rest/active pattern", {
  counts <- .rest_active()
  st <- sleep.sadeh(counts)
  expect_equal(length(st), length(counts))
  rest <- .rest_mask()
  expect_gt(mean(st[rest] == "S"), 0.8)
  expect_gt(mean(st[!rest] == "W"), 0.8)
})

test_that("Webster rescoring only converts sleep to wake", {
  counts <- .rest_active()
  raw <- sleep.cole.kripke(counts, apply_rescoring = FALSE)
  res <- sleep.cole.kripke(counts, apply_rescoring = TRUE)
  expect_lte(sum(res == "S"), sum(raw == "S"))
  # a wake epoch is never turned into sleep by rescoring
  expect_true(all(res[raw == "W"] == "W"))
})

test_that("the scored state feeds sleep.regularity.index", {
  agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
  st <- sleep.cole.kripke(agd$axis1)
  expect_equal(length(st), nrow(agd))
  sri <- sleep.regularity.index(st, agd$timestamp)
  expect_true(is.finite(sri))
})

test_that("scorers handle NA and empty input", {
  s <- sleep.cole.kripke(c(0, NA, 400, 0))   # na_action = "na": NA propagates, no warning
  expect_equal(length(s), 4L)
  expect_true(is.na(s[2]))
  expect_error(sleep.cole.kripke(numeric(0)), "empty")
  expect_error(sleep.sadeh(numeric(0)), "empty")
})
