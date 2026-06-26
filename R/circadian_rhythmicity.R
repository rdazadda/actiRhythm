#' Cosinor Rhythmicity Test
#'
#' Tests whether a rhythm of the given period is present, using the Halberg
#' zero-amplitude F-test (H0: amplitude = 0) on the single-component cosinor fit,
#' and reports the percent rhythm (the proportion of variance the cosinor explains,
#' R^2). It reuses \code{\link{cosinor.analysis}} as the single fitting engine, so
#' the MESOR/amplitude/acrophase and valid-day gating are identical.
#'
#' @param counts Numeric vector of activity counts (or any activity measure).
#' @param timestamps POSIXct timestamps, one per count.
#' @param period Rhythm period in hours (default 24).
#' @param alpha Significance level for the \code{rhythmic} flag (default 0.05).
#' @param wear_time Optional logical wear-time mask passed to the cosinor fit.
#' @param min_valid_hours Minimum valid hours per day for the cosinor fit.
#' @param cosinor_result Optional existing \code{\link{cosinor.analysis}} result to
#'   test without refitting; when supplied, \code{counts}/\code{timestamps} are
#'   ignored.
#'
#' @return An object of class \code{actiRhythm_rhythmicity}: a list with the F
#'   statistic, numerator/denominator degrees of freedom (\code{df1}, \code{df2}),
#'   \code{p_value}, \code{percent_rhythm} (and \code{r_squared}), a logical
#'   \code{rhythmic} flag, and the cosinor parameters.
#'
#' @references
#' \insertRef{nelson1979}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440)
#' h <- as.numeric(format(ts, "%H"))
#' counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
#' rhythmicity.test(counts, ts)
#'
#' @export
rhythmicity.test <- function(counts, timestamps, period = 24, alpha = 0.05,
                             wear_time = NULL, min_valid_hours = 10,
                             cosinor_result = NULL) {
  fit <- cosinor_result
  if (is.null(fit)) {
    fit <- cosinor.analysis(counts, timestamps, period = period,
                            wear_time = wear_time, min_valid_hours = min_valid_hours)
  }

  r2  <- fit$r_squared
  n   <- fit$n_profile_hours
  df1 <- 2L
  df2 <- if (is.null(n)) NA_integer_ else as.integer(n) - 3L

  # Recompute F/p from R^2 and df so the test is exact (the cosinor fit rounds
  # its own f_statistic/p_value for display). F = (R^2/df1) / ((1-R^2)/df2) is
  # algebraically identical to the weighted SS-based F.
  if (is.na(r2) || is.na(df2) || df2 <= 0) {
    f_stat <- NA_real_; p_value <- NA_real_
  } else if (r2 >= 1) {
    f_stat <- Inf; p_value <- 0
  } else {
    f_stat  <- (r2 / df1) / ((1 - r2) / df2)
    p_value <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)
  }

  result <- list(
    F              = f_stat,
    df1            = df1,
    df2            = df2,
    p_value        = p_value,
    percent_rhythm = if (is.na(r2)) NA_real_ else 100 * r2,
    r_squared      = r2,
    rhythmic       = isTRUE(is.finite(p_value) && p_value < alpha),
    alpha          = alpha,
    amplitude      = fit$amplitude,
    mesor          = fit$mesor,
    acrophase      = fit$acrophase,
    period         = period,
    n              = n,
    method         = "Halberg zero-amplitude cosinor F-test"
  )
  class(result) <- "actiRhythm_rhythmicity"
  result
}


#' @export
print.actiRhythm_rhythmicity <- function(x, ...) {
  cat("Cosinor Rhythmicity Test (Halberg zero-amplitude F-test)\n\n")
  cat("H0: amplitude = 0 (no rhythm)\n")
  cat(sprintf("  Period:         %g h\n", x$period))
  cat(sprintf("  F(%d, %d):        %s\n", x$df1, x$df2,
              formatC(x$F, format = "f", digits = 2)))
  cat(sprintf("  P-value:        %s\n", format.pval(x$p_value, digits = 3, eps = 1e-300)))
  cat(sprintf("  Percent rhythm: %s%% (R-squared = %s)\n",
              formatC(x$percent_rhythm, format = "f", digits = 1),
              formatC(x$r_squared, format = "f", digits = 4)))
  cat(sprintf("  Rhythmic:       %s (alpha = %g)\n",
              if (isTRUE(x$rhythmic)) "YES" else "no", x$alpha))
  invisible(x)
}
