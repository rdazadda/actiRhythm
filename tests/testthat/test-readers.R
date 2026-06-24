# Device-file readers exercised against the bundled sample .agd.

sample_agd <- function() {
  system.file("extdata", "MOS2E3923063660sec.agd", package = "actiRhythm")
}

test_that("a sample .agd ships with the package", {
  expect_true(nzchar(sample_agd()))
  expect_true(file.exists(sample_agd()))
})

test_that("read.agd + agd.counts ingest a real .agd into counts", {
  agd <- read.agd(sample_agd(), verbose = FALSE)
  expect_false(is.null(agd))

  cnt <- agd.counts(agd)
  expect_true(all(c("axis1", "timestamp") %in% names(cnt)))
  expect_gt(length(cnt$axis1), 0)
  expect_true(is.numeric(cnt$axis1))
  expect_gt(sum(cnt$axis1, na.rm = TRUE), 0)
  expect_s3_class(cnt$timestamp, "POSIXct")
})

test_that("circadian.rhythm runs on counts read from a .agd file", {
  cnt <- agd.counts(read.agd(sample_agd(), verbose = FALSE))
  r <- circadian.rhythm(cnt$axis1, cnt$timestamp)

  expect_true(all(c("IS", "IV", "RA") %in% names(r)))
  expect_true(r$IS >= 0 && r$IS <= 1)
  expect_true(r$IV >= 0)
  expect_true(r$RA >= 0 && r$RA <= 1)
})

test_that("gt3x.counts is available and guards its optional dependencies", {
  expect_true(is.function(gt3x.counts))
  # Without agcounts it must fail with an actionable message, not a cryptic one.
  if (!requireNamespace("agcounts", quietly = TRUE)) {
    expect_error(gt3x.counts("missing.gt3x"), "agcounts")
  }
})
