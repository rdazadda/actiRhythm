# Bingham et al. (1982) Eq 64 population acrophase CI. The joint Fieller-t interval
# must be symmetric about the acrophase; the c23 sign error broke that symmetry.
test_that("population acrophase CI follows Bingham Eq 64", {
  M  <- c(200, 205, 198, 210, 195, 202, 199, 207, 201, 196, 203, 208)
  b1 <- c(-70, -60, -75, -55, -68, -72, -58, -65, -71, -62, -69, -57)
  b2 <- c(38, 45, 30, 50, 40, 35, 48, 42, 33, 47, 39, 44)
  res <- actiRhythm:::.bingham.population(M, b1, b2)
  expect_true(res$conf_interval_valid)
  expect_gte(res$acrophase, res$ci_acrophase[1])
  expect_lte(res$acrophase, res$ci_acrophase[2])
  lo <- res$acrophase - res$ci_acrophase[1]
  hi <- res$ci_acrophase[2] - res$acrophase
  expect_lt(abs(lo - hi), 0.05)
  expect_equal(round(res$ci_acrophase, 3), c(9.580, 10.127))
})

test_that("population acrophase CI stays contiguous near midnight", {
  M  <- c(200, 205, 198, 210, 195, 202, 199, 207, 201, 196, 203, 208)
  b1 <- c(75, 72, 78, 70, 76, 74, 71, 77, 73, 79, 72, 75)
  b2 <- c(18, 22, 15, 25, 19, 16, 23, 20, 17, 21, 18, 24)
  res <- actiRhythm:::.bingham.population(M, b1, b2)
  expect_true(res$conf_interval_valid)
  expect_lt(res$acrophase %% 24, 2)                 # acrophase just after midnight
  expect_gte(res$acrophase, res$ci_acrophase[1])    # contiguous band brackets it
  expect_lte(res$acrophase, res$ci_acrophase[2])
  expect_lt(diff(res$ci_acrophase), 6)              # not the complementary ~22 h arc
})
