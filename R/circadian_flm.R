#' Functional Linear Model of the 24-Hour Activity Profile
#'
#' Fits the averaged 24-hour activity profile with a periodic basis expansion
#' (Fourier by default, B-spline alternative) by weighted least squares, giving a
#' smooth functional form of the daily activity pattern. The single-component
#' cosinor is the one-harmonic special case; adding harmonics fits the
#' non-sinusoidal shape of a real rest-activity profile. This follows the
#' functional form of Wang et al. (2011), as implemented for actigraphy by pyActigraphy.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param basis Basis type: \code{"fourier"} (default) or \code{"bspline"}.
#' @param n_harmonics Number of Fourier harmonics (default 4, giving a 9-term
#'   expansion as in Wang et al. (2011)).
#' @param nbasis Number of B-spline basis functions (default 9).
#' @param spline_order B-spline order (default 4, cubic).
#' @param period Profile period in hours (default 24).
#' @param wear_time Optional logical wear-time mask.
#' @param min_valid_hours Minimum valid hours per day (default 10).
#' @param weights Profile weighting: \code{"n"} (square root of the per-hour
#'   observation count, default, matching \code{\link{cosinor.analysis}}) or
#'   \code{"none"} (plain least squares, matching pyActigraphy).
#' @param n_eval Length of the dense within-day evaluation grid (default 1440).
#'
#' @return An object of class \code{actiRhythm_flm}: the fitted
#'   \code{coefficients}, the smooth daily curve (\code{smooth_curve}), the fitted
#'   profile, per-harmonic amplitudes and acrophases (\code{harmonics}, Fourier
#'   only), the peak and trough, and the fit statistics (\code{r_squared},
#'   \code{aic}, \code{f_statistic}, \code{p_value}). The function never errors;
#'   on insufficient data it returns the same structure with \code{r_squared} NA.
#'
#' @references
#' \insertRef{wang2011flm}{actiRhythm}
#'
#' \insertRef{ramsay2005}{actiRhythm}
#'
#' \insertRef{hammad2021}{actiRhythm}
#'
#' @seealso \code{\link{cosinor.analysis}}, \code{\link{circadian.ssa}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
#' counts <- 100 + 60 * cos(2 * pi * (h - 14) / 24) + 25 * cos(2 * pi * (h - 6) / 12)
#' circadian.flm(counts, ts)
#'
#' @export
circadian.flm <- function(counts, timestamps, basis = c("fourier", "bspline"),
                          n_harmonics = 4, nbasis = 9, spline_order = 4,
                          period = 24, wear_time = NULL, min_valid_hours = 10,
                          weights = c("n", "none"), n_eval = 1440) {
  basis   <- match.arg(basis)
  weights <- match.arg(weights)
  nbasis  <- max(nbasis, spline_order)

  na_result <- function() structure(list(
    basis = basis, order = if (basis == "fourier") n_harmonics else nbasis,
    period = period, coefficients = numeric(0),
    fitted_profile = data.frame(t = numeric(0), observed = numeric(0),
                                n = numeric(0), fitted = numeric(0)),
    smooth_curve = data.frame(t = numeric(0), activity = numeric(0)),
    r_squared = NA_real_, aic = NA_real_, f_statistic = NA_real_, p_value = NA_real_,
    harmonics = NULL, peak_time = NA_real_, peak_value = NA_real_,
    trough_time = NA_real_, trough_value = NA_real_, n_days = 0L, n_used = 0L,
    n_profile_points = 0L, weights = weights,
    method = paste0("flm_", basis)), class = c("actiRhythm_flm", "list"))

  if (length(counts) != length(timestamps) || !length(counts)) return(na_result())

  build_basis <- function(t) {
    if (basis == "fourier") {
      omega <- 2 * pi / period
      cbind(1, do.call(cbind, lapply(seq_len(n_harmonics),
        function(k) cbind(cos(k * omega * t), sin(k * omega * t)))))
    } else {
      n_interior <- nbasis - spline_order
      interior <- if (n_interior > 0)
        seq(0, period, length.out = n_interior + 2)[-c(1, n_interior + 2)] else numeric(0)
      knots <- c(rep(0, spline_order), interior, rep(period, spline_order))
      splines::splineDesign(knots, x = t, ord = spline_order, outer.ok = TRUE)
    }
  }

  counts <- .gate_invalid_days(counts, timestamps, wear_time, min_valid_hours)
  if (!is.null(wear_time)) counts[!wear_time] <- NA

  prof <- tryCatch(.ext.hourly.profile(counts, timestamps), error = function(e) NULL)
  if (is.null(prof) || length(prof$t) < 12L || stats::sd(prof$y) == 0) return(na_result())

  t <- prof$t; y <- prof$y; nobs <- prof$n
  X <- build_basis(t)
  w <- if (weights == "n") sqrt(nobs) else rep(1, length(nobs))

  fit <- tryCatch(stats::lm.fit(X * w, y * w), error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(fit$coefficients))) return(na_result())
  coef <- fit$coefficients

  t_dense <- seq(0, period, length.out = n_eval + 1)[-(n_eval + 1)]
  smooth  <- as.numeric(build_basis(t_dense) %*% coef)
  fitted  <- as.numeric(X %*% coef)

  ybar_w   <- sum(w^2 * y) / sum(w^2)
  ss_total <- sum(w^2 * (y - ybar_w)^2)
  ss_resid <- sum(w^2 * (y - fitted)^2)
  r_squared <- if (ss_total > 0) 1 - ss_resid / ss_total else NA_real_
  p <- ncol(X); n_grid <- length(y)
  df_model <- p - 1L; df_resid <- n_grid - p
  f_statistic <- if (df_resid > 0 && ss_resid > 0)
    ((ss_total - ss_resid) / df_model) / (ss_resid / df_resid) else NA_real_
  p_value <- if (is.finite(f_statistic))
    stats::pf(f_statistic, df_model, df_resid, lower.tail = FALSE) else NA_real_
  aic <- if (ss_resid > 0) n_grid * log(ss_resid / n_grid) + 2 * p else NA_real_

  harmonics <- NULL
  if (basis == "fourier") {
    k <- seq_len(n_harmonics)
    cos_c <- coef[2 * k]; sin_c <- coef[2 * k + 1]
    harmonics <- data.frame(
      harmonic = k,
      frequency_cycles_per_day = k * (24 / period),
      amplitude = sqrt(cos_c^2 + sin_c^2),
      acrophase_hours = (atan2(sin_c, cos_c) * period / (2 * pi)) %% period)
  }

  pk <- which.max(smooth); tr <- which.min(smooth)

  structure(list(
    basis = basis, order = if (basis == "fourier") n_harmonics else nbasis,
    period = period, coefficients = coef,
    fitted_profile = data.frame(t = t, observed = y, n = nobs, fitted = fitted),
    smooth_curve = data.frame(t = t_dense, activity = smooth),
    r_squared = r_squared, aic = aic, f_statistic = f_statistic, p_value = p_value,
    harmonics = harmonics, peak_time = t_dense[pk], peak_value = smooth[pk],
    trough_time = t_dense[tr], trough_value = smooth[tr],
    n_days = length(unique(as.Date(timestamps))), n_used = sum(!is.na(counts)),
    n_profile_points = n_grid, weights = weights, method = paste0("flm_", basis)),
    class = c("actiRhythm_flm", "list"))
}


#' @export
print.actiRhythm_flm <- function(x, ...) {
  cat("Functional Linear Model (24-hour activity profile)\n\n")
  if (is.na(x$r_squared)) { cat("  Insufficient data for FLM fit\n\n"); return(invisible(x)) }
  cat(sprintf("  Basis:        %s (order %d, %d terms)\n",
              x$basis, x$order, length(x$coefficients)))
  cat(sprintf("  Period:       %g hours\n", x$period))
  cat(sprintf("  Days / profile points: %d / %d\n", x$n_days, x$n_profile_points))
  cat("\n  Model fit:\n")
  cat(sprintf("    R-squared:   %.4f\n", x$r_squared))
  cat(sprintf("    AIC:         %.2f\n", x$aic))
  cat(sprintf("    F-statistic: %.2f (p = %s)\n",
              x$f_statistic, format.pval(x$p_value, digits = 2)))
  if (!is.null(x$harmonics)) {
    h <- x$harmonics[which.max(x$harmonics$amplitude), ]
    cat(sprintf("\n  Dominant harmonic: H%d, amplitude %.2f, acrophase %.2f h\n",
                h$harmonic, h$amplitude, h$acrophase_hours))
  }
  cat(sprintf("\n  Peak %.1f at %.1f h; trough %.1f at %.1f h\n",
              x$peak_value, x$peak_time, x$trough_value, x$trough_time))
  cat("\n  Reference: Wang et al. (2011)\n\n")
  invisible(x)
}
