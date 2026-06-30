# Danilevicz et al. (2024) ML/Bayesian transition probabilities, eps = 0.5.
# Series ends in rest: rest at-risk = 8 - 1 = 7, active at-risk = 3.
test_that("transition.probability applies final-bout censoring with eps = 0.5", {
  counts <- c(0, 0, 0, 5, 5, 0, 0, 5, 0, 0, 0)
  tp <- transition.probability(counts, threshold = 1)
  expect_equal(tp$tp_ra_mle, 2 / 7)
  expect_equal(tp$tp_ra_bayes, 2.5 / 7.5)
  expect_equal(tp$tp_ar_mle, 2 / 3)
  expect_equal(tp$tp_ar_bayes, 2.5 / 3.5)
})
