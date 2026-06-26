# The bibliography lives in two places by necessity: inst/REFERENCES.bib is what
# Rdpack (\insertRef in the man pages) reads, and vignettes/REFERENCES.bib is what
# the vignettes and pkgdown articles cite (a relative path into inst/ breaks under
# pkgdown's temp-dir rendering). This guards the two copies against drift; it runs
# from the source tree and skips when only the installed package is available.
test_that("vignettes/REFERENCES.bib stays in sync with inst/REFERENCES.bib", {
  skip_on_cran()
  inst_bib <- testthat::test_path("..", "..", "inst", "REFERENCES.bib")
  vig_bib  <- testthat::test_path("..", "..", "vignettes", "REFERENCES.bib")
  skip_if_not(file.exists(inst_bib) && file.exists(vig_bib),
              "source bib files not available (installed package)")
  expect_identical(readLines(inst_bib, warn = FALSE),
                   readLines(vig_bib, warn = FALSE))
})
