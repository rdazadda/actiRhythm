#' @title Endogenous Circadian Period Estimation (Lomb-Scargle)
#'
#' @description
#' Estimates the dominant (endogenous) circadian PERIOD (tau) of an activity
#' time series using the Lomb-Scargle periodogram. Unlike the classical Fast
#' Fourier Transform (FFT), the Lomb-Scargle method does not require evenly
#' sampled data, so it correctly handles the irregular and gappy sampling that
#' results from non-wear periods, dropped epochs, or mixed epoch lengths in
#' accelerometer recordings.
#'
#' @name circadian-period
NULL


# Standard-normalized Lomb-Scargle periodogram over a PERIOD window, with the
# Baluev (2008) analytic false-alarm probability. Replicates the relevant path
# of lomb::lsp(type = "period", normalize = "standard") so actiRhythm does not
# depend on lomb (and its plotly/data.table closure). Inputs are assumed finite
# and length-matched; returns NULL on a degenerate grid/series.
.lomb_scargle <- function(x, times, from, to, ofac = 4) {
  ofac <- max(1L, as.integer(floor(ofac)))
  o <- order(times)
  t <- as.numeric(times)[o]
  y <- as.numeric(x)[o]
  n <- length(y)
  tspan <- t[n] - t[1]
  if (!is.finite(tspan) || tspan <= 0) return(NULL)

  fr.d <- 1 / tspan
  step <- 1 / (tspan * ofac)
  # type = "period": the period window [from, to] maps to frequencies [1/to, 1/from].
  freq <- seq(fr.d, 1 / from, by = step)
  freq <- freq[freq >= 1 / to]
  n.out <- length(freq)
  if (n.out < 2L) return(NULL)

  y <- y - mean(y)
  ss <- sum(y^2)
  if (!is.finite(ss) || ss == 0) return(NULL)
  norm <- 1 / ss

  w <- 2 * pi * freq
  PN <- numeric(n.out)
  for (i in seq_len(n.out)) {
    wi <- w[i]
    tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t))) / wi
    arg <- wi * (t - tau)
    cs <- cos(arg)
    sn <- sin(arg)
    PN[i] <- (sum(y * cs))^2 / sum(cs * cs) + (sum(y * sn))^2 / sum(sn * sn)
  }
  PN <- norm * PN
  PN.max <- max(PN)
  peak.freq <- freq[which.max(PN)]

  # Baluev (2008) false-alarm probability (lomb's pbaluev/ggamma).
  ggamma <- function(N) sqrt(2 / N) * exp(lgamma(N / 2) - lgamma((N - 1) / 2))
  Dt <- mean(t^2) - mean(t)^2
  NH <- n - 1
  NK <- n - 3
  fsingle <- (1 - PN.max)^(0.5 * NK)
  W <- max(freq) * sqrt(4 * pi * Dt)
  tau_b <- ggamma(NH) * W * (1 - PN.max)^(0.5 * (NK - 1)) * sqrt(0.5 * NH * PN.max)
  p.value <- -(exp(-tau_b) - 1) + fsingle * exp(-tau_b)

  # Report on the ascending-period axis (lomb reverses the frequency order here).
  list(
    scanned = (1 / freq)[n.out:1],
    power   = PN[n.out:1],
    peak    = PN.max,
    peak.at = c(1 / peak.freq, peak.freq),
    p.value = p.value,
    n.out   = n.out,
    n       = n
  )
}


#' Estimate Endogenous Circadian Period via the Lomb-Scargle Periodogram
#'
#' Computes the Lomb-Scargle periodogram of an activity \code{counts} series
#' sampled at the supplied \code{timestamps} and returns the period (tau, in
#' hours) of the strongest spectral peak within the search window
#' \code{[from, to]}. The Lomb-Scargle method (Lomb 1976; Scargle 1982) is the
#' least-squares-equivalent spectral estimator for unevenly sampled time series
#' and is therefore appropriate for actigraphy data containing gaps, which the
#' FFT cannot accommodate.
#'
#' @param counts Numeric vector of activity counts (minute-level recommended).
#'   \code{NA} values (e.g. non-wear epochs) are dropped together with their
#'   timestamps before estimation.
#' @param timestamps A \code{POSIXct} vector (or anything coercible by
#'   \code{as.numeric}) of epoch timestamps, the same length as \code{counts}.
#'   Internally converted to hours elapsed since the first timestamp.
#' @param from Numeric. Lower bound of the period search window, in hours
#'   (default \code{18}).
#' @param to Numeric. Upper bound of the period search window, in hours
#'   (default \code{30}).
#' @param ofac Integer oversampling factor controlling the period-grid
#'   resolution. Higher values give a finer period grid and a more precise peak
#'   location at the cost of computation (default \code{4}).
#'
#' @return A named \code{list} with elements:
#'   \describe{
#'     \item{tau}{Numeric. Period (hours) of the strongest Lomb-Scargle peak in
#'       \code{[from, to]}, i.e. the estimated endogenous circadian period.
#'       \code{NA_real_} when the data are insufficient.}
#'     \item{peak_power}{Numeric. Normalized power of that peak (the
#'       Lomb-Scargle peak statistic). \code{NA_real_} when insufficient.}
#'     \item{p_value}{Numeric. P-value of the peak under the null hypothesis of
#'       Gaussian noise (Baluev 2008 analytic false-alarm probability).
#'       \code{NA_real_} when insufficient.}
#'     \item{oversampling}{The \code{ofac} oversampling factor used.}
#'     \item{n_used}{Integer. Number of non-\code{NA} observations actually
#'       passed to the periodogram (\code{NA_integer_} when not run).}
#'     \item{span_days}{Numeric. Total recording span in days (max minus min
#'       timestamp), used for the >= 2-day guard.}
#'     \item{scanned}{Numeric vector of trial periods (hours) of the full
#'       Lomb-Scargle spectrum (\code{numeric(0)} when not run).}
#'     \item{power}{Numeric vector of standard-normalized Lomb-Scargle power,
#'       aligned to \code{scanned} (\code{numeric(0)} when not run).}
#'   }
#'   On any edge case (too few points, too short a span, degenerate input, or an
#'   internal numerical failure) the function returns this same structure
#'   with \code{tau}, \code{peak_power} and \code{p_value} set to \code{NA};
#'   it never throws.
#'
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Timestamps are converted to hours since the first sample
#'     (\code{t_hours = (as.numeric(timestamps) - min) / 3600}).
#'   \item Pairs with a missing \code{count} or a missing/non-finite time are
#'     dropped.
#'   \item Two guards are applied so the estimate is never based on
#'     insufficient data: the recording must span at least \strong{2 days}
#'     (otherwise a 18-30 h period cannot be resolved) and at least
#'     \strong{10 non-\code{NA}} observations must remain.
#'   \item The standard-normalized Lomb-Scargle periodogram is evaluated over
#'     the period window. Its strongest peak gives the period (hours), and the
#'     Baluev (2008) analytic false-alarm probability gives the p-value.
#' }
#'
#' The Lomb-Scargle periodogram is chosen specifically because actigraphy series
#' are rarely gap-free: the FFT assumes uniform sampling, whereas Lomb-Scargle
#' fits sinusoids by least squares at each trial frequency and is unbiased for
#' irregular sampling.
#'
#' @references
#' \insertRef{lomb1976}{actiRhythm}
#'
#' \insertRef{scargle1982}{actiRhythm}
#'
#' \insertRef{ruf1999}{actiRhythm}
#'
#' \insertRef{refinetti2007}{actiRhythm}
#'
#' @seealso \code{\link{cosinor.analysis}} for parametric (fixed-period) rhythm
#'   estimation, \code{\link{circadian.rhythm}} for non-parametric L5/M10/IS/IV
#'   metrics.
#'
#' @examples
#' \donttest{
#' # Seven days of minute-level data with a 24 h rhythm
#' t_hours <- seq(0, 7 * 24 - 1/60, by = 1/60)
#' ts <- as.POSIXct("2024-01-01 00:00:00") + t_hours * 3600
#' counts <- 100 + 80 * cos(2 * pi * (t_hours - 8) / 24) + rnorm(length(t_hours), 0, 5)
#' circadian.period(counts, ts)$tau   # about 24.0
#' }
#'
#' @export
circadian.period <- function(counts, timestamps, from = 18, to = 30, ofac = 4) {

  # Structured "insufficient data" return so callers get a stable shape and the
  # function never errors on an edge case.
  na_result <- function(n_used = NA_integer_, span_days = NA_real_) {
    list(
      tau          = NA_real_,
      peak_power   = NA_real_,
      p_value      = NA_real_,
      oversampling = ofac,
      n_used       = n_used,
      span_days    = span_days,
      scanned      = numeric(0),
      power        = numeric(0)
    )
  }

  # Basic structural validation (never error; return NA structure)
  if (missing(counts) || missing(timestamps) ||
      length(counts) == 0L || length(timestamps) == 0L ||
      length(counts) != length(timestamps)) {
    return(na_result())
  }

  # Numeric time axis in hours since the first sample. as.numeric() works for
  # POSIXct (seconds since epoch) and for plain numeric seconds.
  t_sec <- suppressWarnings(as.numeric(timestamps))
  counts <- suppressWarnings(as.numeric(counts))

  # Drop NA / non-finite pairs in either counts or times
  keep <- is.finite(counts) & is.finite(t_sec)
  counts <- counts[keep]
  t_sec  <- t_sec[keep]

  n_used <- length(counts)
  if (n_used == 0L) {
    return(na_result(n_used = 0L))
  }

  # Hours elapsed since the first retained sample.
  t_hours <- (t_sec - min(t_sec)) / 3600

  span_days <- (max(t_hours) - min(t_hours)) / 24

  # Guards: require >= 2 days of span and >= 10 non-NA points
  # A period in the 18-30 h band cannot be resolved from < 2 days of data, and
  # too few samples make the periodogram meaningless.
  if (!is.finite(span_days) || span_days < 2) {
    return(na_result(n_used = n_used, span_days = span_days))
  }
  if (n_used < 10L) {
    return(na_result(n_used = n_used, span_days = span_days))
  }

  # A degenerate (constant) series has no spectral peak; guard so lsp does not
  # divide by a zero variance.
  if (stats::sd(counts) == 0 || !is.finite(stats::sd(counts))) {
    return(na_result(n_used = n_used, span_days = span_days))
  }

  # Lomb-Scargle periodogram over the requested period window
  # Native standard-normalized Lomb-Scargle; degrades to the NA structure on any
  # failure rather than throwing.
  lsp <- tryCatch(
    .lomb_scargle(x = counts, times = t_hours, from = from, to = to, ofac = ofac),
    error = function(e) NULL
  )

  if (is.null(lsp)) {
    return(na_result(n_used = n_used, span_days = span_days))
  }

  # peak.at[1] is the peak PERIOD (hours); its second element is the matching
  # frequency, which we ignore.
  tau        <- suppressWarnings(as.numeric(lsp$peak.at[1]))
  peak_power <- suppressWarnings(as.numeric(lsp$peak))
  p_value    <- suppressWarnings(as.numeric(lsp$p.value))

  if (length(tau) == 0L || !is.finite(tau)) {
    return(na_result(n_used = n_used, span_days = span_days))
  }

  list(
    tau          = tau,
    peak_power   = if (length(peak_power) == 0L) NA_real_ else peak_power,
    p_value      = if (length(p_value) == 0L) NA_real_ else p_value,
    oversampling = ofac,
    n_used       = n_used,
    span_days    = span_days,
    scanned      = as.numeric(lsp$scanned),
    power        = as.numeric(lsp$power)
  )
}


# Sub-grid peak location: fit a parabola through the argmax and its two
# neighbours on the ascending period axis and return the vertex. This is
# load-bearing for period.ci(): snapping a bootstrap replicate's peak to the
# nearest grid node quantizes the period to the grid step and collapses the
# bootstrap SE to zero for strong rhythms. Falls back to the grid peak at a
# boundary or a non-concave bracket.
.lsp_peak_parabola <- function(scanned, power) {
  if (length(power) < 3L) {
    ip <- which.max(power)
    return(list(tau = scanned[ip], power = power[ip], at_boundary = TRUE))
  }
  ip <- which.max(power)
  if (ip == 1L || ip == length(power)) {
    return(list(tau = scanned[ip], power = power[ip], at_boundary = TRUE))
  }
  x1 <- scanned[ip - 1L]; x2 <- scanned[ip]; x3 <- scanned[ip + 1L]
  y1 <- power[ip - 1L];   y2 <- power[ip];   y3 <- power[ip + 1L]
  d <- (x1 - x2) * (x1 - x3) * (x2 - x3)
  if (!is.finite(d) || d == 0) return(list(tau = x2, power = y2, at_boundary = FALSE))
  a <- (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / d
  b <- (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / d
  if (!is.finite(a) || a >= 0) return(list(tau = x2, power = y2, at_boundary = FALSE))
  cc <- (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 +
         x1 * x2 * (x1 - x2) * y3) / d
  list(tau = -b / (2 * a), power = cc - b^2 / (4 * a), at_boundary = FALSE)
}


#' Confidence Interval for the Endogenous Circadian Period
#'
#' Attaches a bootstrap confidence interval (and standard error) to the
#' Lomb-Scargle period estimate from \code{\link{circadian.period}}. Because
#' activity is strongly autocorrelated, an i.i.d. resample would destroy the
#' rhythm and give an invalid (too-narrow) interval, so a \emph{circular
#' moving-block} bootstrap of the cosinor residuals is used: the rhythm fitted at
#' the point estimate is held fixed while blocks of residuals are resampled, the
#' periodogram peak is re-located on each replicate (with sub-grid parabolic
#' refinement), and the interval is taken from the quantiles of the replicate
#' periods.
#'
#' @param counts Numeric activity vector; \code{NA} are dropped with their times.
#' @param timestamps POSIXct (or numeric-coercible) timestamps, same length.
#' @param from,to Period search window in hours (default 18, 30).
#' @param ofac Integer oversampling factor for the period grid (default 4).
#' @param n_boot Number of bootstrap replicates (default 200).
#' @param block_hours Moving-block length in hours (default 24).
#' @param level Confidence level (default 0.95).
#' @param seed Optional integer seed passed to \code{\link[base]{set.seed}} for
#'   reproducible bootstrap draws.
#'
#' @return An object of class \code{actiRhythm_period_ci}: a list with \code{tau}
#'   (the point estimate, hours), \code{ci_lower}/\code{ci_upper}, \code{se}, the
#'   \code{level}, and the number of valid replicates.
#'
#' @references
#' \insertRef{kunsch1989}{actiRhythm}
#'
#' \insertRef{politis1992}{actiRhythm}
#'
#' @seealso \code{\link{circadian.period}}
#'
#' @examples
#' \donttest{
#' t_hours <- seq(0, 5 * 24 - 1/60, by = 1/60)
#' ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
#' counts <- 100 + 80 * cos(2 * pi * (t_hours - 8) / 24) + rnorm(length(t_hours), 0, 20)
#' period.ci(counts, ts, n_boot = 50, seed = 1)
#' }
#'
#' @export
period.ci <- function(counts, timestamps, from = 18, to = 30, ofac = 4,
                      n_boot = 200, block_hours = 24, level = 0.95, seed = NULL) {

  na_out <- function(tau = NA_real_) structure(
    list(tau = tau, ci_lower = NA_real_, ci_upper = NA_real_, se = NA_real_,
         level = level, n_boot = n_boot, n_valid = 0L, block_hours = block_hours,
         method = "circular block residual bootstrap"),
    class = "actiRhythm_period_ci")

  t_sec <- suppressWarnings(as.numeric(timestamps))
  cnt   <- suppressWarnings(as.numeric(counts))
  keep  <- is.finite(cnt) & is.finite(t_sec)
  cnt <- cnt[keep]; t_sec <- t_sec[keep]
  n <- length(cnt)
  if (n < 10L) return(na_out())

  o <- order(t_sec); cnt <- cnt[o]; t_sec <- t_sec[o]
  t_hours <- (t_sec - t_sec[1]) / 3600
  span_h <- t_hours[n] - t_hours[1]
  if (!is.finite(span_h) || span_h / 24 < 2 || stats::sd(cnt) == 0) return(na_out())

  lsp0 <- tryCatch(.lomb_scargle(cnt, t_hours, from, to, ofac), error = function(e) NULL)
  if (is.null(lsp0)) return(na_out())
  tau_hat <- .lsp_peak_parabola(lsp0$scanned, lsp0$power)$tau
  if (!is.finite(tau_hat)) return(na_out())

  if (!is.null(seed)) set.seed(seed)

  # Hold the rhythm fitted at tau_hat fixed; bootstrap blocks of its residuals.
  omega  <- 2 * pi / tau_hat
  X      <- cbind(1, cos(omega * t_hours), sin(omega * t_hours))
  fit    <- stats::lm.fit(X, cnt)
  fitted <- as.numeric(X %*% fit$coefficients)
  resid  <- cnt - fitted

  L <- max(2L, min(as.integer(round(block_hours * n / span_h)), as.integer(n / 2)))
  nblocks <- ceiling(n / L)

  tau_boot <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    starts <- sample.int(n, nblocks, replace = TRUE)
    idx <- unlist(lapply(starts, function(s) ((s - 1L + 0:(L - 1L)) %% n) + 1L))[seq_len(n)]
    lspb <- tryCatch(.lomb_scargle(fitted + resid[idx], t_hours, from, to, ofac),
                     error = function(e) NULL)
    if (!is.null(lspb)) tau_boot[b] <- .lsp_peak_parabola(lspb$scanned, lspb$power)$tau
  }
  tau_boot <- tau_boot[is.finite(tau_boot)]
  if (length(tau_boot) < 2L) return(na_out(tau_hat))

  # Basic (pivotal) bootstrap interval: reflects the replicate quantiles around
  # the point estimate, correcting the bias of a plain percentile interval.
  a  <- (1 - level) / 2
  q  <- stats::quantile(tau_boot, c(a, 1 - a), names = FALSE)
  ci <- sort(c(2 * tau_hat - q[2], 2 * tau_hat - q[1]))
  structure(list(
    tau = tau_hat, ci_lower = ci[1], ci_upper = ci[2], se = stats::sd(tau_boot),
    level = level, n_boot = n_boot, n_valid = length(tau_boot),
    block_hours = block_hours, method = "circular block residual bootstrap"
  ), class = "actiRhythm_period_ci")
}


#' @export
print.actiRhythm_period_ci <- function(x, ...) {
  cat("Circadian Period with Bootstrap Confidence Interval\n\n")
  cat(sprintf("  tau:      %s h\n", formatC(x$tau, format = "f", digits = 3)))
  cat(sprintf("  %g%% CI:   [%s, %s] h\n", 100 * x$level,
              formatC(x$ci_lower, format = "f", digits = 3),
              formatC(x$ci_upper, format = "f", digits = 3)))
  cat(sprintf("  SE:       %s h\n", formatC(x$se, format = "f", digits = 3)))
  cat(sprintf("  Method:   %s (%d/%d valid reps)\n", x$method, x$n_valid, x$n_boot))
  invisible(x)
}
