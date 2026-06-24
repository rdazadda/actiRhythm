# Circadian workbook export + plot save helper.

test_that("circadian.workbook writes a multi-sheet xlsx including the new methods", {
  skip_if_not_installed("openxlsx")
  set.seed(1)
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
  h  <- as.numeric(format(ts, "%H"))
  act <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24) + stats::rnorm(length(ts), 0, 20))

  f <- tempfile(fileext = ".xlsx"); on.exit(unlink(f), add = TRUE)
  wb <- circadian.workbook(act, ts, file = f, n_boot = 40)

  expect_true(file.exists(f))
  sheets <- openxlsx::getSheetNames(f)
  expect_true(all(c("Summary", "Hourly Profile", "MF-DFA", "Data Dictionary") %in% sheets))
  expect_true(any(grepl("Transitions", sheets)))
  expect_true(any(grepl("Periodogram", sheets)))

  smry <- openxlsx::read.xlsx(f, sheet = "Summary")
  # core + every new-method column present and populated
  expect_true(all(c("IS", "rhythm_F", "period_ci_lower", "mfdfa_h2", "kRA", "mse_area") %in% names(smry)))
  expect_true(is.finite(smry$IS))         # nonparametric engine ran
  expect_true(is.finite(smry$rhythm_F))   # rhythmicity test
  expect_true(is.finite(smry$mfdfa_h2))   # MF-DFA
  expect_true(is.finite(smry$kRA))        # state transitions
})

test_that("circadian.workbook errors clearly without openxlsx", {
  skip_if(requireNamespace("openxlsx", quietly = TRUE), "openxlsx installed")
  ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2880)
  expect_error(circadian.workbook(stats::rnorm(2880), ts), "openxlsx")
})

test_that("save.circadian.plot writes an image file", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)
  save.circadian.plot(p, f, width = 4, height = 3)
  expect_true(file.exists(f))
})
