# Raw multi-brand ingest: gt3x (hardened public-API path), Axivity, GENEActiv.
# Uses the sample files bundled with read.gt3x and GGIRread; skips if absent.

test_that("gt3x.counts reproduces the agcounts baseline via the public API", {
  skip_if_not_installed("agcounts"); skip_if_not_installed("read.gt3x")
  g <- system.file("extdata", "TAS1H30182785_2019-09-17.gt3x", package = "read.gt3x")
  skip_if(!nzchar(g))
  cc <- gt3x.counts(g)
  expect_true(all(c("time", "axis1", "axis2", "axis3", "vm") %in% names(cc)))
  expect_s3_class(cc$time, "POSIXct")
  expect_equal(sum(cc$axis1), 27065)   # byte-parity with the old private-internals path
})

test_that("axivity.counts reads a .cwa and returns epoch counts", {
  skip_if_not_installed("GGIRread"); skip_if_not_installed("agcounts")
  cwa <- system.file("testfiles", "ax3_testfile.cwa", package = "GGIRread")
  skip_if(!nzchar(cwa))
  ax <- axivity.counts(cwa, epoch = 15)
  expect_true(all(c("time", "axis1", "axis2", "axis3", "vm") %in% names(ax)))
  expect_gt(nrow(ax), 0L)
  expect_true(all(ax$vm >= 0))
})

test_that("geneactiv.counts reads a .bin and resamples 85.7Hz to 30Hz", {
  skip_if_not_installed("GGIRread"); skip_if_not_installed("agcounts")
  bin <- system.file("testfiles", "GENEActiv_testfile.bin", package = "GGIRread")
  skip_if(!nzchar(bin))
  ge <- geneactiv.counts(bin, epoch = 15)
  expect_true(all(c("time", "axis1", "axis2", "axis3", "vm") %in% names(ge)))
  expect_gt(nrow(ge), 0L)
  # a recording shorter than one epoch errors clearly, not cryptically
  expect_error(geneactiv.counts(bin, epoch = 60), "shorter than one")
})

test_that("read.raw dispatches by file extension", {
  skip_if_not_installed("GGIRread"); skip_if_not_installed("agcounts")
  skip_if_not_installed("read.gt3x")
  g   <- system.file("extdata", "TAS1H30182785_2019-09-17.gt3x", package = "read.gt3x")
  cwa <- system.file("testfiles", "ax3_testfile.cwa", package = "GGIRread")
  skip_if(!nzchar(g) || !nzchar(cwa))
  expect_identical(read.raw(g), gt3x.counts(g))
  expect_equal(nrow(read.raw(cwa, epoch = 15)), nrow(axivity.counts(cwa, epoch = 15)))
  expect_error(read.raw("recording.foo"), "cannot infer device")
})

test_that(".resample_to_30hz lands on a 30Hz grid", {
  fs <- 85.7; n <- round(fs * 10)                       # 10 seconds of data
  xyz <- data.frame(X = sin(seq_len(n) / fs), Y = cos(seq_len(n) / fs), Z = rep(1, n))
  out <- actiRhythm:::.resample_to_30hz(xyz, fs)
  expect_gt(nrow(out), 295L); expect_lt(nrow(out), 305L)  # about 10s * 30Hz
  expect_true(all(is.finite(out$X)) && all(is.finite(out$Z)))
})

test_that(".locf_idle_sleep carries gravity through idle-sleep zero runs", {
  # validated against GGIR: read.gt3x imputes idle-sleep as zeros, which otherwise
  # collapse the z-angle to the calibration-offset angle; LOCF restores it.
  g <- actiRhythm:::.locf_idle_sleep(c(0.1, 0.2, 0, 0, 0.3),
                                     c(0.9, 0.9, 0, 0, 0.9),
                                     c(0.4, 0.4, 0, 0, 0.4))
  expect_equal(g$x, c(0.1, 0.2, 0.2, 0.2, 0.3))         # carried forward
  expect_equal(g$z, c(0.4, 0.4, 0.4, 0.4, 0.4))
  g2 <- actiRhythm:::.locf_idle_sleep(c(0, 0, 0.5), c(0, 0, 0.5), c(0, 0, 0.5))
  expect_equal(g2$x, c(0, 0, 0.5))                      # leading zeros stay zero
})
