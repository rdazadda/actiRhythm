#' Activity Counts from a Raw .gt3x File
#'
#' Compute activity counts from a raw ActiGraph \code{.gt3x} accelerometer file
#' using the agcounts implementation of the ActiGraph count algorithm (Neishabouri
#' 2022). Requires the \pkg{agcounts} and \pkg{read.gt3x} packages.
#'
#' @param path Path to a \code{.gt3x} file.
#' @param epoch Epoch length in seconds (default 60).
#' @param lfe Use the low-frequency extension filter (default \code{FALSE}).
#' @param tz Time zone for the timestamps (default \code{"UTC"}).
#'
#' @return Data frame with \code{time}, \code{axis1}, \code{axis2}, \code{axis3}
#'   and \code{vm}, one row per epoch.
#'
#' @references
#' Neishabouri A, et al. (2022). Quantification of acceleration as activity counts
#' in ActiGraph wearable. \emph{Scientific Reports}, 12:11958.
#'
#' @seealso \code{\link{axivity.counts}}, \code{\link{geneactiv.counts}},
#'   \code{\link{read.raw}}
#'
#' @export
gt3x.counts <- function(path, epoch = 60, lfe = FALSE, tz = "UTC") {
  if (!requireNamespace("agcounts", quietly = TRUE)) {
    stop("gt3x.counts() requires the 'agcounts' package: install.packages('agcounts')")
  }
  if (!requireNamespace("read.gt3x", quietly = TRUE)) {
    stop("gt3x.counts() requires the 'read.gt3x' package")
  }
  raw <- agcounts::agread(path, parser = "read.gt3x", tz = tz, verbose = FALSE)
  .raw_to_counts(raw, epoch, lfe, tz)
}
