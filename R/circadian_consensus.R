#' Consensus Rhythmicity Across Methods
#'
#' Combines the rhythmicity verdicts of several independent tests into one
#' call: the cosinor zero-amplitude F-test (\code{\link{rhythmicity.test}}), the
#' Bingham confidence ellipse (\code{\link{cosinor.confidence.ellipse}}), the
#' Lomb-Scargle Baluev false-alarm probability (\code{\link{circadian.period}}),
#' and the chi-square (Sokolove-Bushell) periodogram (\code{\link{chi.sq.periodogram}}).
#' The available p-values are pooled by Fisher's method, and a majority vote across
#' all methods is also reported.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Rhythm period in hours for the cosinor tests (default 24).
#' @param alpha Significance level (default 0.05).
#' @param wear_time Optional logical wear-time mask for the cosinor fit.
#'
#' @return An object of class \code{actiRhythm_consensus}: the Fisher-combined
#'   p-value and consensus call, the vote count, a per-method \code{tests} data
#'   frame, and the agreement fraction.
#'
#' @references
#' \insertRef{fisher1925}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h <- as.numeric(format(ts, "%H"))
#' counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
#' consensus.rhythmicity(counts, ts)
#'
#' @export
consensus.rhythmicity <- function(counts, timestamps, period = 24, alpha = 0.05,
                                  wear_time = NULL) {
  saf <- function(expr) tryCatch(expr, error = function(e) NULL)

  cos <- saf(cosinor.analysis(counts, timestamps, period = period, wear_time = wear_time))
  rhy <- saf(rhythmicity.test(counts, timestamps, cosinor_result = cos, alpha = alpha))
  ell <- saf(cosinor.confidence.ellipse(cos))
  per <- saf(circadian.period(counts, timestamps))
  chi <- saf(chi.sq.periodogram(counts, timestamps))

  p_cos <- if (is.null(rhy)) NA_real_ else rhy$p_value
  p_lsp <- if (is.null(per)) NA_real_ else per$p_value
  p_chi <- if (is.null(chi)) NA_real_ else chi$p_value

  tests <- data.frame(
    method = c("cosinor F-test", "Bingham ellipse", "Lomb-Scargle FAP", "chi-square periodogram"),
    p_value = c(p_cos, NA_real_, p_lsp, p_chi),
    rhythmic = c(
      if (is.na(p_cos)) NA else p_cos < alpha,
      if (is.null(ell)) NA else isTRUE(ell$rhythm_detected),
      if (is.na(p_lsp)) NA else p_lsp < alpha,
      if (is.null(chi)) NA else isTRUE(chi$significant)
    ),
    stringsAsFactors = FALSE
  )

  # Fisher combination of the available, valid p-values.
  ps <- c(p_cos, p_lsp, p_chi)
  ps <- pmin(pmax(ps[is.finite(ps)], 1e-300), 1)
  fisher_p <- if (length(ps) >= 1) {
    stats::pchisq(-2 * sum(log(ps)), df = 2 * length(ps), lower.tail = FALSE)
  } else NA_real_

  votes   <- sum(tests$rhythmic, na.rm = TRUE)
  n_tests <- sum(!is.na(tests$rhythmic))

  structure(list(
    consensus_p = fisher_p,
    consensus_rhythmic = isTRUE(is.finite(fisher_p) && fisher_p < alpha),
    votes = votes, n_tests = n_tests,
    agreement = if (n_tests > 0) votes / n_tests else NA_real_,
    tests = tests, alpha = alpha, period = period
  ), class = "actiRhythm_consensus")
}


#' @export
print.actiRhythm_consensus <- function(x, ...) {
  cat("Consensus Rhythmicity (multi-method)\n\n")
  t <- x$tests
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %-22s %-9s %s\n", t$method[i],
                if (is.na(t$p_value[i])) "" else paste0("p=", format.pval(t$p_value[i], digits = 2)),
                if (is.na(t$rhythmic[i])) "-" else if (t$rhythmic[i]) "rhythmic" else "no"))
  }
  cat(sprintf("\n  Votes:        %d / %d methods\n", x$votes, x$n_tests))
  cat(sprintf("  Fisher p:     %s\n", format.pval(x$consensus_p, digits = 3)))
  cat(sprintf("  Consensus:    %s (alpha = %g)\n",
              if (isTRUE(x$consensus_rhythmic)) "RHYTHMIC" else "not rhythmic", x$alpha))
  invisible(x)
}
