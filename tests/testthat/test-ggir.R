# Cross-check of the raw pipeline against GGIR on a real .gt3x recording kept in
# data-raw/ (which is excluded from the CRAN tarball via .Rbuildignore). This is
# slow (GGIR re-reads and processes the whole file, ~5 min) and needs GGIR plus the
# recording, so it is opt-in: set RUN_GGIR_VALIDATION=true to run it. The same
# comparison renders live in the "From raw acceleration" article whenever the source
# repository and GGIR are present.
test_that("raw metrics agree with GGIR on the validation recording", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("RUN_GGIR_VALIDATION")),
              "set RUN_GGIR_VALIDATION=true to run the (slow) GGIR cross-check")
  skip_if_not_installed("GGIR")
  skip_if_not_installed("read.gt3x")
  gt3x <- testthat::test_path("..", "..", "data-raw", "MOS2E39230594.gt3x")
  skip_if_not(file.exists(gt3x), "validation .gt3x not present (data-raw/ is not shipped)")

  P <- GGIR::load_params()
  P$params_general$windowsizes <- c(5, 900, 3600)
  I <- GGIR::g.inspectfile(gt3x, params_rawdata = P$params_rawdata,
                           params_general = P$params_general)
  ggir <- GGIR::g.getmeta(gt3x, params_rawdata = P$params_rawdata,
                          params_general = P$params_general,
                          params_cleaning = P$params_cleaning,
                          inspectfileobject = I)$metashort
  acti <- raw.metrics(gt3x, epoch = 5)

  gg <- data.frame(k = gsub("T", " ", substr(as.character(ggir$timestamp), 1, 19)),
                   enmo = ggir$ENMO * 1000, anglez = ggir$anglez)
  aa <- data.frame(k = format(acti$time, "%Y-%m-%d %H:%M:%S"),
                   enmo = acti$ENMO, anglez = acti$anglez)
  m <- merge(gg, aa, by = "k", suffixes = c("_ggir", "_acti"))

  expect_gt(nrow(m), 1000)
  expect_gt(cor(m$enmo_ggir,   m$enmo_acti),   0.95)
  expect_gt(cor(m$anglez_ggir, m$anglez_acti), 0.95)
  expect_lt(mean(abs(m$enmo_ggir   - m$enmo_acti)), 10)   # mg
  expect_lt(mean(abs(m$anglez_ggir - m$anglez_acti)), 5)  # degrees

  # Sleep-period window: same z-angle into rest.spt() and GGIR's HASPT (HDCZA).
  ps   <- P$params_sleep
  spt  <- rest.spt(acti$anglez, acti$time, epoch_length = 5)
  day  <- as.Date(acti$time - 12 * 3600)
  gspt <- do.call(rbind, lapply(unique(day), function(dd) {
    idx <- which(day == dd)
    if (length(idx) < 720) return(NULL)
    h <- GGIR:::HASPT(acti$anglez[idx], params_sleep = ps, ws3 = 5,
                      HASPT.algo = "HDCZA", invalid = rep(0L, length(idx)))
    if (length(h$SPTE_start) == 0 || is.na(h$SPTE_start)) return(NULL)
    data.frame(date = as.Date(dd, origin = "1970-01-01"),
               g_on = acti$time[idx[round(h$SPTE_start)]],
               g_off = acti$time[idx[round(h$SPTE_end)]])
  }))
  cmp <- merge(spt, gspt, by = "date")
  expect_gt(nrow(cmp), 3)
  expect_lt(max(abs(as.numeric(difftime(cmp$onset,  cmp$g_on,  units = "secs")))), 30)
  expect_lt(max(abs(as.numeric(difftime(cmp$offset, cmp$g_off, units = "secs")))), 30)
})
