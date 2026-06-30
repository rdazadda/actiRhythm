# Bridge from raw files to the existing circadian battery, plus two named
# estimators (ABI, transition probability) that the raw-data literature uses.

#' Circadian Analysis Directly from Raw Acceleration
#'
#' From a raw accelerometer file (or raw data frame), computes a per-epoch raw
#' metric (ENMO or MAD) with auto-calibration, then runs the count-style circadian
#' analysis (\code{\link{circadian.rhythm}}) on it. The metrics (IS, IV, RA,
#' L5/M10, cosinor, DFA, periodograms, SRI, ...) run on the raw metric the same way
#' they run on counts, so raw data needs no new per-method code.
#'
#' @param x A path to a raw file (\code{.gt3x}, \code{.cwa}, \code{.bin}) or a raw
#'   data frame (see \code{\link{raw.metrics}} / \code{\link{example_raw}}).
#' @param metric Raw metric to analyse, \code{"ENMO"} (default) or \code{"MAD"}.
#' @param device Device brand or \code{"auto"} (default; used for file input).
#' @param epoch Epoch length in seconds (default 60).
#' @param calibrate Apply van Hees auto-calibration first (default \code{TRUE}).
#' @param tz Time zone (default \code{"UTC"}).
#' @param ... Passed to \code{\link{circadian.rhythm}}.
#'
#' @return The \code{\link{circadian.rhythm}} result computed on the raw metric.
#'
#' @seealso \code{\link{raw.metrics}}, \code{\link{circadian.rhythm}},
#'   \code{\link{example_raw}}
#'
#' @examples
#' # The full count-style battery on synthetic raw ENMO
#' \donttest{
#' cr <- circadian.raw(example_raw(days = 2), metric = "ENMO")
#' c(IS = cr$IS, IV = cr$IV, RA = cr$RA)
#' }
#'
#' @export
circadian.raw <- function(x, metric = c("ENMO", "MAD"), device = "auto",
                          epoch = 60, calibrate = TRUE, tz = "UTC", ...) {
  metric <- match.arg(metric)
  rm <- raw.metrics(x, device = device, epoch = epoch, metrics = metric,
                    calibrate = calibrate, tz = tz)
  circadian.rhythm(rm[[metric]], rm$time, epoch_length = epoch, ...)
}

#' Activity Balance Index
#'
#' The Activity Balance Index (Danilevicz et al. 2024), a 0 to 1 transform of a
#' detrended fluctuation analysis scaling exponent that peaks at the healthy
#' \eqn{\alpha = 1} (1/f) balance: \eqn{ABI(\alpha) = \exp(-|\alpha - 1| / e^{-2}) = \exp(-e^{2}\,|\alpha - 1|)}.
#'
#' @param x Either a numeric scaling exponent, or a fractal object with
#'   \code{alpha} (and optionally \code{alpha1}, \code{alpha2}) such as the result
#'   of the package's detrended fluctuation analysis.
#'
#' @return If \code{x} is numeric, the scalar ABI. If \code{x} is a fractal
#'   object, a list with \code{ABI_overall}, \code{ABI_short}, \code{ABI_long}.
#'
#' @references
#' \insertRef{danilevicz2024}{actiRhythm}
#'
#' @examples
#' activity.balance.index(1.0)   # perfect 1/f balance -> 1
#' activity.balance.index(0.7)
#'
#' @export
activity.balance.index <- function(x) {
  abi <- function(a) if (is.null(a) || !is.finite(a)) NA_real_ else exp(-abs(a - 1) / exp(-2))
  if (is.numeric(x)) return(abi(x[1]))
  list(ABI_overall = abi(x$alpha), ABI_short = abi(x$alpha1), ABI_long = abi(x$alpha2))
}
