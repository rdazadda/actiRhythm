# Flag epochs in zero-runs of at least `frame` length, given a logical activity
# mask `active` (TRUE = nonzero/worn). Returns a non-wear logical vector.
.nonwear_from_zeroruns <- function(active, frame) {
  n <- length(active)
  nw <- logical(n)
  r <- rle(!active); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  for (i in which(r$values)) if (r$lengths[i] >= frame) nw[starts[i]:ends[i]] <- TRUE
  nw
}

#' Choi (2011) Non-Wear Detection
#'
#' Classifies each epoch as wear or non-wear with the Choi et al. (2011)
#' algorithm: a run of consecutive zero-count epochs of at least \code{frame}
#' minutes is non-wear, tolerating a short nonzero spike only when it is flanked
#' by a fully zero window of \code{stream} minutes both before and after. The
#' returned mask can be passed as the \code{wear_time} argument to the
#' rest-activity and sleep functions.
#'
#' @param counts Numeric activity vector (vertical axis), minute epochs assumed.
#' @param epoch_length Epoch length in seconds (default 60). Window lengths are
#'   given in minutes and scaled to epochs by this value.
#' @param frame Minimum non-wear window in minutes (default 90).
#' @param spike_tolerance Maximum tolerated nonzero spike in minutes (default 2).
#' @param stream Flanking all-zero window required around a tolerated spike, in
#'   minutes (default 30).
#'
#' @return A logical vector, one per epoch: \code{TRUE} = wear, \code{FALSE} =
#'   non-wear. Never errors.
#'
#' @references
#' \insertRef{choi2011}{actiRhythm}
#'
#' @seealso \code{\link{detect.nonwear.troiano}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 600)
#' counts <- c(rep(200, 200), rep(0, 200), rep(200, 200))   # 200-min non-wear gap
#' table(detect.nonwear.choi(counts))
#'
#' @export
detect.nonwear.choi <- function(counts, epoch_length = 60, frame = 90,
                                spike_tolerance = 2, stream = 30) {
  x <- suppressWarnings(as.numeric(counts)); n <- length(x)
  if (n == 0L) return(logical(0))
  x[is.na(x)] <- 0
  per_min <- 60 / epoch_length
  fr <- round(frame * per_min); sp <- max(1L, round(spike_tolerance * per_min))
  st <- round(stream * per_min)
  active <- x > 0
  r <- rle(active); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  for (i in which(r$values)) {                  # each activity run
    if (r$lengths[i] <= sp) {                    # short enough to be a spike
      s <- starts[i]; e <- ends[i]
      up   <- if (s > 1L) sum(x[max(s - st, 1L):(s - 1L)]) else 0
      down <- if (e < n)  sum(x[(e + 1L):min(e + st, n)])   else 0
      if (up == 0 && down == 0) active[s:e] <- FALSE   # absorb into the zero run
    }
  }
  !.nonwear_from_zeroruns(active, fr)
}

#' Troiano (2008) Non-Wear Detection
#'
#' Classifies each epoch as wear or non-wear with the Troiano et al. (2008,
#' NHANES) algorithm: a run of at least \code{frame} minutes of zero counts is
#' non-wear, tolerating up to \code{spike_tolerance} consecutive nonzero minutes
#' only if every one of them is at or below \code{stoplevel} counts. Unlike Choi
#' it has no flanking-window requirement but applies a count ceiling on spikes.
#'
#' @param counts Numeric activity vector (vertical axis), minute epochs assumed.
#' @param epoch_length Epoch length in seconds (default 60).
#' @param frame Minimum non-wear window in minutes (default 60).
#' @param spike_tolerance Maximum tolerated nonzero spike in minutes (default 2).
#' @param stoplevel Count above which a spike ends the non-wear bout (default 100).
#'
#' @return A logical vector, one per epoch: \code{TRUE} = wear, \code{FALSE} =
#'   non-wear. Never errors.
#'
#' @references
#' \insertRef{troiano2008}{actiRhythm}
#'
#' @seealso \code{\link{detect.nonwear.choi}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 400)
#' counts <- c(rep(150, 100), rep(0, 200), rep(150, 100))
#' table(detect.nonwear.troiano(counts))
#'
#' @export
detect.nonwear.troiano <- function(counts, epoch_length = 60, frame = 60,
                                   spike_tolerance = 2, stoplevel = 100) {
  x <- suppressWarnings(as.numeric(counts)); n <- length(x)
  if (n == 0L) return(logical(0))
  x[is.na(x)] <- 0
  per_min <- 60 / epoch_length
  fr <- round(frame * per_min); sp <- max(1L, round(spike_tolerance * per_min))
  active <- x > 0
  r <- rle(active); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  for (i in which(r$values)) {
    idx <- starts[i]:ends[i]
    if (r$lengths[i] <= sp && max(x[idx]) <= stoplevel) active[idx] <- FALSE
  }
  !.nonwear_from_zeroruns(active, fr)
}
