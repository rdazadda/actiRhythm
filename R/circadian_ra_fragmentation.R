#' Rest-Activity Bout Fragmentation
#'
#' Summarizes how broken up the rest-activity rhythm is, from a per-epoch
#' rest/active state: the mean and median rest and active bout durations, the
#' number of state transitions, and transitions per day. These add a
#' bout-length view of fragmentation to the transition probabilities of
#' \code{\link{state.transitions}} (kRA/kAR) (Lim et al. 2011). This covers
#' rest-activity-rhythm fragmentation only; it omits the
#' sedentary-behaviour bout distribution (Gini, power law, hazard), which
#' is a physical-activity-epidemiology concern, not a circadian one.
#'
#' @param state Per-epoch state: a logical vector (TRUE = active) or a character
#'   vector where \code{"R"}/\code{"S"}/\code{"sleep"}/\code{"rest"} mark rest.
#' @param timestamps POSIXct timestamps, one per value.
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return An object of class \code{actiRhythm_rafrag}: mean/median rest and
#'   active bout durations (minutes), bout counts, transition count, and
#'   transitions per day. Never errors.
#'
#' @references
#' \insertRef{lim2011}{actiRhythm}
#'
#' @seealso \code{\link{state.transitions}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' rest.activity.fragmentation(h >= 7 & h < 23, ts)
#'
#' @export
rest.activity.fragmentation <- function(state, timestamps, epoch_length = 60) {
  if (length(state) != length(timestamps))
    stop("state and timestamps must have same length")
  active <- if (is.logical(state)) state
            else !(tolower(as.character(state)) %in% c("r", "s", "sleep", "rest"))
  ok <- !is.na(active) & !is.na(timestamps)
  active <- active[ok]; ts <- timestamps[ok]
  em <- epoch_length / 60

  na_out <- structure(list(
    mean_active_bout = NA_real_, median_active_bout = NA_real_,
    mean_rest_bout = NA_real_, median_rest_bout = NA_real_,
    n_active_bouts = 0L, n_rest_bouts = 0L, n_transitions = 0L,
    transitions_per_day = NA_real_), class = c("actiRhythm_rafrag", "list"))
  if (length(active) < 2L) return(na_out)

  r <- rle(active)
  active_bouts <- r$lengths[r$values] * em
  rest_bouts <- r$lengths[!r$values] * em
  span_days <- as.numeric(difftime(ts[length(ts)], ts[1], units = "days"))
  n_trans <- length(r$lengths) - 1L

  structure(list(
    mean_active_bout = if (length(active_bouts)) mean(active_bouts) else NA_real_,
    median_active_bout = if (length(active_bouts)) stats::median(active_bouts) else NA_real_,
    mean_rest_bout = if (length(rest_bouts)) mean(rest_bouts) else NA_real_,
    median_rest_bout = if (length(rest_bouts)) stats::median(rest_bouts) else NA_real_,
    n_active_bouts = length(active_bouts), n_rest_bouts = length(rest_bouts),
    n_transitions = n_trans,
    transitions_per_day = if (span_days > 0) n_trans / span_days else NA_real_),
    class = c("actiRhythm_rafrag", "list"))
}

#' @export
print.actiRhythm_rafrag <- function(x, ...) {
  cat("Rest-Activity Bout Fragmentation\n\n")
  cat(sprintf("  Active bouts: %d, mean %.0f min (median %.0f)\n",
              x$n_active_bouts, x$mean_active_bout, x$median_active_bout))
  cat(sprintf("  Rest bouts:   %d, mean %.0f min (median %.0f)\n",
              x$n_rest_bouts, x$mean_rest_bout, x$median_rest_bout))
  cat(sprintf("  Transitions:  %d (%.1f per day)\n\n",
              x$n_transitions, x$transitions_per_day))
  invisible(x)
}
