# Between-group cosinor comparison (Bingham summary-statistics approach).

build_two_groups <- function(acroA = 8, acroB = 12, mesorA = 100, mesorB = 100,
                             nper = 6, sd_noise = 3, seed = 1) {
  set.seed(seed)
  hrs <- 0:23
  base <- as.POSIXct("2024-01-01", tz = "UTC")
  act <- numeric(0); subj <- character(0); grp <- character(0)
  ts <- as.POSIXct(character(0), tz = "UTC")
  for (g in c("A", "B")) {
    acro <- if (g == "A") acroA else acroB
    M    <- if (g == "A") mesorA else mesorB
    for (i in seq_len(nper)) {
      y <- M + 40 * cos(2 * pi * (hrs - acro) / 24) + stats::rnorm(24, 0, sd_noise)
      act  <- c(act, y)
      subj <- c(subj, rep(paste0(g, i), 24))
      grp  <- c(grp, rep(g, 24))
      ts   <- c(ts, base + hrs * 3600)
    }
  }
  list(activity = act, timestamps = ts, subject = subj, group = grp)
}

test_that("detects a known acrophase difference but no MESOR difference", {
  d <- build_two_groups(acroA = 8, acroB = 12, mesorA = 100, mesorB = 100)
  r <- cosinor.compare(d$activity, d$timestamps, d$subject, d$group)

  expect_s3_class(r, "actiRhythm_cosinor_compare")
  acro <- r$tests[r$tests$parameter == "acrophase", ]
  mes  <- r$tests[r$tests$parameter == "mesor", ]
  expect_lt(acro$p_value, 0.05)
  expect_equal(abs(acro$difference), 4, tolerance = 1)
  expect_gt(mes$p_value, 0.05)                 # MESOR equal in both groups
})

test_that("the parameter tests match base t.test exactly", {
  d <- build_two_groups(acroA = 9, acroB = 9, mesorA = 100, mesorB = 120, seed = 2)
  r <- cosinor.compare(d$activity, d$timestamps, d$subject, d$group)

  cf <- r$subjects
  tt <- stats::t.test(cf$mesor[cf$group == "A"], cf$mesor[cf$group == "B"])
  row <- r$tests[r$tests$parameter == "mesor", ]
  expect_equal(row$statistic, unname(tt$statistic), tolerance = 1e-9)
  expect_equal(row$p_value, tt$p.value, tolerance = 1e-9)
  expect_lt(row$p_value, 0.05)                 # MESOR differs (100 vs 120)
})

test_that("errors when group does not have exactly two levels", {
  d <- build_two_groups()
  g3 <- d$group; g3[1] <- "C"                  # introduce a third level
  expect_error(cosinor.compare(d$activity, d$timestamps, d$subject, g3), "two levels")
})

test_that("acrophase comparison handles the midnight wrap", {
  d <- build_two_groups(acroA = 23.5, acroB = 0.5, sd_noise = 2, seed = 4)
  r <- cosinor.compare(d$activity, d$timestamps, d$subject, d$group)
  acro <- r$tests[r$tests$parameter == "acrophase", ]
  expect_lt(abs(acro$difference), 2)           # about 1 h apart, not about 23 h
})
