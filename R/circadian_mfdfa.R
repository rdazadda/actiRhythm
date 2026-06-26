#' Multifractal Detrended Fluctuation Analysis (MF-DFA)
#'
#' Generalizes detrended fluctuation analysis (\code{\link{fractal.dfa}}) to a
#' spectrum of moment orders q (Kantelhardt et al. 2002). A flat generalized
#' Hurst exponent h(q) indicates a monofractal signal; a decreasing h(q) and a
#' wide multifractal spectrum indicate multifractality. h(2) equals the standard
#' DFA scaling exponent.
#'
#' @param x Numeric series (e.g. activity counts); the longest gap-free run is
#'   analysed.
#' @param scale_min Smallest window size in samples (default 8).
#' @param scale_max Largest window size (default \code{floor(N/4)}).
#' @param q_values Moment orders to evaluate (default \code{seq(-5, 5, 0.5)}).
#' @param both_ends If \code{TRUE} (default) windows are taken from both ends so
#'   the tail of the profile is not discarded; \code{FALSE} reproduces the
#'   start-only convention of \code{\link{fractal.dfa}}.
#'
#' @return An object of class \code{actiRhythm_mfdfa}: a list with \code{q_values},
#'   \code{h_q} (generalized Hurst exponent), \code{tau_q} (mass exponent),
#'   \code{alpha}/\code{f_alpha} (the multifractal spectrum), \code{alpha_dfa}
#'   (\code{= h(2)}) and the spectrum \code{width}. Returns an NA structure on
#'   insufficient data; never errors.
#'
#' @references
#' \insertRef{kantelhardt2002}{actiRhythm}
#'
#' @seealso \code{\link{fractal.dfa}}
#'
#' @examples
#' set.seed(1)
#' mfdfa(stats::rnorm(4096))$alpha_dfa   # near 0.5 for white noise
#'
#' @export
mfdfa <- function(x, scale_min = 8L, scale_max = NULL,
                  q_values = seq(-5, 5, by = 0.5), both_ends = TRUE) {

  na_out <- function() structure(list(
    q_values = q_values, h_q = rep(NA_real_, length(q_values)),
    tau_q = rep(NA_real_, length(q_values)), alpha = rep(NA_real_, length(q_values)),
    f_alpha = rep(NA_real_, length(q_values)), alpha_dfa = NA_real_,
    width = NA_real_, scales = integer(0), both_ends = both_ends),
    class = "actiRhythm_mfdfa")

  x <- .longest_non_na_run(as.numeric(x))
  N <- length(x)
  if (N < 4L * scale_min) return(na_out())
  if (is.null(scale_max)) scale_max <- floor(N / 4)
  if (scale_max <= scale_min) return(na_out())

  scales <- unique(round(exp(seq(log(scale_min), log(scale_max), length.out = 20))))
  scales <- scales[scales >= scale_min & scales <= scale_max]
  if (length(scales) < 4L) return(na_out())

  Y <- cumsum(x - mean(x))

  Fqs <- matrix(NA_real_, nrow = length(scales), ncol = length(q_values))
  for (si in seq_along(scales)) {
    f2 <- .mfdfa_variances(Y, scales[si], both_ends)
    f2 <- f2[is.finite(f2) & f2 > 0]
    if (length(f2) < 2L) next
    for (qi in seq_along(q_values)) {
      q <- q_values[qi]
      Fqs[si, qi] <- if (abs(q) < 1e-8) exp(mean(log(f2)) / 2) else (mean(f2^(q / 2)))^(1 / q)
    }
  }

  logs <- log10(scales)
  h_q <- vapply(seq_along(q_values), function(qi) {
    fq <- Fqs[, qi]
    ok <- is.finite(fq) & fq > 0
    if (sum(ok) < 3L) return(NA_real_)
    unname(stats::lm.fit(cbind(1, logs[ok]), log10(fq[ok]))$coefficients[2])
  }, numeric(1))

  tau_q   <- q_values * h_q - 1
  alpha_h <- .central_diff(q_values, tau_q)        # Holder exponent = d tau / d q
  f_alpha <- q_values * alpha_h - tau_q

  structure(list(
    q_values = q_values, h_q = h_q, tau_q = tau_q,
    alpha = alpha_h, f_alpha = f_alpha,
    alpha_dfa = h_q[which.min(abs(q_values - 2))],
    width = if (any(is.finite(alpha_h))) max(alpha_h, na.rm = TRUE) - min(alpha_h, na.rm = TRUE) else NA_real_,
    scales = scales, both_ends = both_ends
  ), class = "actiRhythm_mfdfa")
}


# Detrended variance of each window at one scale (linear detrend, m = 1).
.mfdfa_variances <- function(Y, s, both_ends) {
  N <- length(Y)
  Ns <- floor(N / s)
  if (Ns < 1L) return(numeric(0))
  starts <- (0:(Ns - 1L)) * s
  if (both_ends) starts <- c(starts, N - (1:Ns) * s)
  Xs <- cbind(1, seq_len(s))
  vapply(starts, function(st) mean(stats::lm.fit(Xs, Y[(st + 1L):(st + s)])$residuals^2),
         numeric(1))
}

# Central-difference derivative dy/dx (forward/backward at the ends).
.central_diff <- function(x, y) {
  n <- length(x)
  d <- rep(NA_real_, n)
  if (n < 2L) return(d)
  d[1] <- (y[2] - y[1]) / (x[2] - x[1])
  d[n] <- (y[n] - y[n - 1L]) / (x[n] - x[n - 1L])
  if (n > 2L) {
    i <- 2:(n - 1L)
    d[i] <- (y[i + 1L] - y[i - 1L]) / (x[i + 1L] - x[i - 1L])
  }
  d
}


#' @export
print.actiRhythm_mfdfa <- function(x, ...) {
  cat("Multifractal Detrended Fluctuation Analysis\n\n")
  cat(sprintf("  DFA scaling exponent h(2): %s\n",
              formatC(x$alpha_dfa, format = "f", digits = 3)))
  cat(sprintf("  h(q) range:                [%s, %s]\n",
              formatC(min(x$h_q, na.rm = TRUE), format = "f", digits = 3),
              formatC(max(x$h_q, na.rm = TRUE), format = "f", digits = 3)))
  cat(sprintf("  Multifractal spectrum width: %s\n",
              formatC(x$width, format = "f", digits = 3)))
  cat(sprintf("  Scales: %d (%d-%d samples)\n",
              length(x$scales), min(x$scales), max(x$scales)))
  invisible(x)
}
