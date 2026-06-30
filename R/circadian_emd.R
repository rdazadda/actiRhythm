# Indices of interior local maxima and minima (sign change of the slope).
.emd_extrema <- function(x) {
  s <- diff(sign(diff(x)))
  list(maxi = which(s < 0) + 1L, mini = which(s > 0) + 1L)
}

# Sift one intrinsic mode function out of the residual h. Stops on the Wu & Huang
# (2009) Cauchy-type energy ratio sum((h - h1)^2) / sum(h^2); the ~1e-3 threshold
# gives proper convergence (the 0.2-0.3 range belongs to Huang's per-point SD, a
# different statistic).
.emd_sift_one <- function(x, max_sift = 50L, cauchy_tol = 1e-3) {
  h <- x; n <- length(h)
  for (it in seq_len(max_sift)) {
    e <- .emd_extrema(h)
    if (length(e$maxi) < 2L || length(e$mini) < 2L) return(list(imf = h, residue = TRUE))
    # clamp endpoints so the spline envelopes do not overshoot the edges
    mx <- unique(c(1L, e$maxi, n)); mn <- unique(c(1L, e$mini, n))
    up <- stats::spline(mx, h[mx], xout = seq_len(n), method = "natural")$y
    lo <- stats::spline(mn, h[mn], xout = seq_len(n), method = "natural")$y
    h1 <- h - (up + lo) / 2
    cauchy <- sum((h - h1)^2) / (sum(h^2) + 1e-10)
    h <- h1
    if (cauchy < cauchy_tol) break
  }
  list(imf = h, residue = FALSE)
}

# Mean period (hours) of an IMF from its zero-crossing count.
.imf_period <- function(imf, epoch_seconds) {
  zc <- sum(abs(diff(sign(imf))) > 0)
  if (zc < 2) return(NA_real_)
  (length(imf) / (zc / 2)) * epoch_seconds / 3600
}

#' Empirical Mode Decomposition of the Activity Rhythm
#'
#' Decomposes the activity series into intrinsic mode functions by empirical mode
#' decomposition, a data-adaptive, nonlinear alternative to the linear SSA
#' decomposition (Huang et al. 1998). The intrinsic mode function whose period is
#' nearest 24 hours is taken as the circadian component. Optional ensemble
#' EMD adds noise to reduce mode mixing (Wu and Huang 2009).
#'
#' @param counts Numeric activity vector (a coarse epoch is recommended for speed).
#' @param timestamps POSIXct timestamps, one per value.
#' @param max_imf Maximum number of modes to extract (default 10).
#' @param ensemble Ensemble size for EEMD (default 1 = plain EMD).
#' @param noise_sd Added-noise SD as a fraction of the series SD (default 0.2).
#' @param period_range Period window (hours) for the circadian mode (default 20 to 28).
#' @param epoch_length Epoch length in seconds (default 60).
#' @param seed Optional seed for the EEMD noise.
#'
#' @return An object of class \code{actiRhythm_emd}: the IMF matrix, the residual
#'   trend, per-IMF period and variance share, the circadian IMF index, and the
#'   reconstruction error. The first and last epochs are unreliable (the spline
#'   envelopes pin the IMFs near zero at the edges), so read the instantaneous
#'   series away from the boundary. Never errors.
#'
#' @references
#' \insertRef{huang1998}{actiRhythm}
#'
#' \insertRef{wuhuang2009}{actiRhythm}
#'
#' @seealso \code{\link{hilbert.huang}}, \code{\link{circadian.ssa}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#'
#' @export
circadian.emd <- function(counts, timestamps, max_imf = 10L, ensemble = 1L,
                          noise_sd = 0.2, period_range = c(20, 28),
                          epoch_length = 60, seed = NULL) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  o <- order(ts); x <- x[o]; ts <- ts[o]
  n <- length(x)
  na_out <- structure(list(imfs = matrix(nrow = 0, ncol = 0), residual = numeric(0),
    periods = numeric(0), variance_share = numeric(0), circadian_imf = NA_integer_,
    circadian_period = NA_real_, recon_error = NA_real_, times = ts,
    epoch_seconds = epoch_length, insufficient = TRUE), class = c("actiRhythm_emd", "list"))
  if (n < 64L || stats::sd(x) == 0) return(na_out)

  one_emd <- function(sig) {
    r <- sig; out <- list()
    for (k in seq_len(max_imf)) {
      s <- .emd_sift_one(r)
      out[[k]] <- s$imf; r <- r - s$imf
      if (s$residue || length(.emd_extrema(r)$maxi) < 2L) break
    }
    list(imfs = out, residual = r)
  }

  xc <- x - mean(x)
  if (ensemble > 1L) {
    if (!is.null(seed)) set.seed(seed)
    sd_x <- stats::sd(xc)
    accum <- NULL; ncols <- NULL; col_counts <- NULL
    for (e in seq_len(ensemble)) {
      em <- one_emd(xc + stats::rnorm(n, 0, noise_sd * sd_x))
      m <- do.call(cbind, em$imfs)
      if (is.null(accum)) { accum <- m; ncols <- ncol(m); col_counts <- rep(1L, ncols) }
      else {
        k <- min(ncols, ncol(m))
        accum[, seq_len(k)] <- accum[, seq_len(k)] + m[, seq_len(k)]
        col_counts[seq_len(k)] <- col_counts[seq_len(k)] + 1L
      }
    }
    imfs_m <- sweep(accum, 2L, col_counts, "/")
    residual <- xc - rowSums(imfs_m)
  } else {
    em <- one_emd(xc)
    imfs_m <- do.call(cbind, em$imfs)
    residual <- em$residual
  }

  periods <- apply(imfs_m, 2, .imf_period, epoch_seconds = epoch_length)
  vshare <- apply(imfs_m, 2, stats::var) / sum(apply(imfs_m, 2, stats::var))
  in_band <- which(periods >= period_range[1] & periods <= period_range[2])
  circ <- if (length(in_band)) in_band[which.min(abs(periods[in_band] - 24))] else NA_integer_
  recon <- max(abs(rowSums(imfs_m) + residual - xc))

  structure(list(imfs = imfs_m, residual = residual, periods = periods,
    variance_share = vshare, circadian_imf = circ,
    circadian_period = if (is.na(circ)) NA_real_ else periods[circ],
    recon_error = recon, times = ts, epoch_seconds = epoch_length,
    insufficient = FALSE), class = c("actiRhythm_emd", "list"))
}

#' @export
print.actiRhythm_emd <- function(x, ...) {
  cat("Empirical Mode Decomposition\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  IMFs: %d   reconstruction error: %.2e\n", ncol(x$imfs), x$recon_error))
  d <- data.frame(IMF = seq_along(x$periods), period_h = round(x$periods, 2),
                  var_share = round(x$variance_share, 3))
  print(d, row.names = FALSE)
  if (!is.na(x$circadian_imf))
    cat(sprintf("\n  Circadian IMF: %d (period %.2f h)\n\n", x$circadian_imf, x$circadian_period))
  else cat("\n  No IMF in the circadian band\n\n")
  invisible(x)
}


# Analytic signal via FFT (one-sided multiplier); Re = y, Im = Hilbert transform.
.analytic_signal <- function(y) {
  N <- length(y); X <- stats::fft(y)
  h <- numeric(N); h[1] <- 1
  if (N %% 2 == 0) { h[N / 2 + 1] <- 1; if (N > 2) h[2:(N / 2)] <- 2 }
  else h[2:((N + 1) / 2)] <- 2
  stats::fft(X * h, inverse = TRUE) / N
}

#' Hilbert-Huang Instantaneous Phase and Frequency
#'
#' Computes the instantaneous amplitude, phase, and period of a circadian
#' intrinsic mode function via its analytic signal (Huang et al. 1998). Where the
#' cosinor gives one acrophase for the whole recording, the instantaneous phase
#' tracks the rhythm cycle by cycle, and the spread of the instantaneous period
#' measures how stationary the circadian band is.
#'
#' @param x A \code{actiRhythm_emd} object, or a numeric activity vector (then
#'   \code{timestamps} is required).
#' @param timestamps POSIXct timestamps (required when \code{x} is numeric).
#' @param imf Which IMF to analyse; default is the circadian IMF of the EMD object.
#' @param epoch_length Epoch length in seconds (default 60; used when \code{x} is numeric).
#'
#' @return An object of class \code{actiRhythm_hht}: the instantaneous amplitude,
#'   phase, and period series, plus the mean instantaneous period and its SD, mean
#'   amplitude and its CV, and the fraction of time the period stays in 20-28 h.
#'   Never errors.
#'
#' @references
#' \insertRef{huang1998}{actiRhythm}
#'
#' @seealso \code{\link{circadian.emd}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' e <- circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#' hilbert.huang(e)
#'
#' @export
hilbert.huang <- function(x, timestamps = NULL, imf = NULL, epoch_length = 60) {
  if (inherits(x, "actiRhythm_emd")) {
    e <- x
    ep <- if (!is.null(e$epoch_seconds)) e$epoch_seconds else epoch_length
    idx <- if (is.null(imf)) e$circadian_imf else imf
    sig <- if (!is.na(idx) && ncol(e$imfs) >= idx) e$imfs[, idx] else NULL
    ts <- e$times
  } else {
    if (is.null(timestamps)) stop("timestamps required when x is numeric")
    e <- circadian.emd(x, timestamps, epoch_length = epoch_length)
    ep <- epoch_length
    idx <- if (is.null(imf)) e$circadian_imf else imf
    sig <- if (!is.na(idx) && ncol(e$imfs) >= idx) e$imfs[, idx] else NULL
    ts <- e$times
  }
  na_out <- structure(list(amplitude = numeric(0), phase = numeric(0),
    period = numeric(0), mean_period = NA_real_, sd_period = NA_real_,
    mean_amplitude = NA_real_, amplitude_cv = NA_real_, frac_in_band = NA_real_,
    insufficient = TRUE), class = c("actiRhythm_hht", "list"))
  if (is.null(sig) || length(sig) < 4L) return(na_out)

  z <- .analytic_signal(sig)
  amp <- Mod(z); phi <- Arg(z)
  d <- diff(phi); d[d > pi] <- d[d > pi] - 2 * pi; d[d < -pi] <- d[d < -pi] + 2 * pi
  phi_u <- c(phi[1], phi[1] + cumsum(d))
  dt_h <- ep / 3600
  inst_freq <- c(NA, diff(phi_u)) / (2 * pi * dt_h)        # cycles per hour
  inst_period <- 1 / inst_freq
  inst_period[!is.finite(inst_period) | inst_period <= 0] <- NA
  structure(list(amplitude = amp, phase = phi_u, period = inst_period,
    mean_period = mean(inst_period, na.rm = TRUE),
    sd_period = stats::sd(inst_period, na.rm = TRUE),
    mean_amplitude = mean(amp), amplitude_cv = stats::sd(amp) / mean(amp),
    frac_in_band = mean(inst_period >= 20 & inst_period <= 28, na.rm = TRUE),
    insufficient = FALSE), class = c("actiRhythm_hht", "list"))
}

#' @export
print.actiRhythm_hht <- function(x, ...) {
  cat("Hilbert-Huang Instantaneous Dynamics\n\n")
  if (isTRUE(x$insufficient)) { cat("  No circadian IMF to analyse\n\n"); return(invisible(x)) }
  cat(sprintf("  Instantaneous period: %.2f h (SD %.2f)\n", x$mean_period, x$sd_period))
  cat(sprintf("  Instantaneous amplitude: %.1f (CV %.2f)\n", x$mean_amplitude, x$amplitude_cv))
  cat(sprintf("  Time in 20-28 h band: %.0f%%\n\n", 100 * x$frac_in_band))
  invisible(x)
}
