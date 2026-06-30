# Pin the Activity Balance Index to Danilevicz et al. (2024) eq.5:
# ABI(alpha) = exp(-|alpha - 1| / exp(-2)) = exp(-e^2 |alpha - 1|).
test_that("activity.balance.index matches Danilevicz 2024 eq.5", {
  expect_equal(activity.balance.index(1.0), 1)
  expect_equal(activity.balance.index(0.7), exp(-exp(2) * 0.3), tolerance = 1e-9)
  expect_equal(round(activity.balance.index(0.7), 5), 0.10897)
  expect_equal(round(activity.balance.index(0.0), 6), 0.000618)
  expect_equal(round(activity.balance.index(2.0), 6), 0.000618)
})
