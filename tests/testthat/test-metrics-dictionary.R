# The metrics dictionary (inst/extdata/metrics_dictionary.csv) is the single source
# of truth for the codebook vignette and the workbook's "Data Dictionary" sheet.
# These tests keep it valid and in step with the metrics the code actually emits.

test_that("metrics dictionary is structurally valid", {
  dict <- read.csv(
    system.file("extdata", "metrics_dictionary.csv", package = "actiRhythm"),
    stringsAsFactors = FALSE
  )

  expect_setequal(
    names(dict),
    c("metric", "producing_function", "family", "definition", "formula_units",
      "range_interpretation", "reference", "output_object", "controlling_arg")
  )

  # required fields are non-empty; metric names are unique
  for (col in c("metric", "definition", "output_object", "producing_function", "family")) {
    expect_false(any(is.na(dict[[col]]) | dict[[col]] == ""),
                 info = paste("empty value in column", col))
  }
  expect_identical(anyDuplicated(dict$metric), 0L)

  # families come from the allowed set
  allowed <- c("io", "nonparametric", "cosinor", "period", "fractal", "phase", "sleep")
  expect_true(all(dict$family %in% allowed),
              info = paste("unknown family:",
                           paste(setdiff(dict$family, allowed), collapse = ", ")))

  # every producing_function is an exported function of the package
  exports <- getNamespaceExports("actiRhythm")
  expect_true(all(unique(dict$producing_function) %in% exports),
              info = paste("not exported:",
                           paste(setdiff(dict$producing_function, exports), collapse = ", ")))

  # every non-empty reference key resolves in the shared bibliography
  keys <- unique(dict$reference[nzchar(dict$reference)])
  bibkeys <- names(Rdpack::get_bibentries(package = "actiRhythm"))
  expect_true(all(keys %in% bibkeys),
              info = paste("unresolved reference key:",
                           paste(setdiff(keys, bibkeys), collapse = ", ")))
})

test_that("dictionary documents exactly the circadian.batch / workbook summary metrics", {
  skip_on_cran()
  dict <- read.csv(
    system.file("extdata", "metrics_dictionary.csv", package = "actiRhythm"),
    stringsAsFactors = FALSE
  )
  ts  <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
  h   <- as.numeric(format(ts, "%H"))
  act <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
  summ <- actiRhythm:::.circadian.analyze(act, ts, include_period_ci = TRUE, n_boot = 10)$summary
  expect_setequal(dict$metric, names(summ))
})
