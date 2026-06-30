#' LIDS Ultradian Cycle Curve
#'
#' Plots the smoothed Locomotor Inactivity During Sleep (LIDS) signal across one
#' sleep period with the best-fit ultradian cosine overlaid, the companion figure
#' to \code{\link{lids}}. The cosine's period and Munich Rhythmicity Index
#' summarise the sleep-cycle oscillation. Returns a \code{ggplot} object and never
#' errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per count.
#' @param sleep_periods Data frame with \code{in_bed_time} and \code{out_bed_time}
#'   (as in \code{\link{lids}}).
#' @param period Which sleep period to plot (default 1).
#' @param ... Passed to \code{\link{lids}}.
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{winnebeck2018}{actiRhythm}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01 22:00", tz = "UTC"), by = 60, length.out = 480)
#' sp <- data.frame(in_bed_time = ts[1], out_bed_time = ts[480])
#' plot_lids(50 + 45 * cos(2 * pi * seq_along(ts) / 110), ts, sp)
#' }
#'
#' @export
plot_lids <- function(counts, timestamps, sleep_periods, period = 1L, ...) {
  empty <- function() .circ_empty_plot("Insufficient data for a LIDS fit",
                                        title = "LIDS ultradian cycle")
  li <- tryCatch(lids(counts, timestamps, sleep_periods, ...), error = function(e) NULL)
  if (is.null(li) || isTRUE(li$insufficient) || !length(li$fits)) return(empty())
  period <- max(1L, min(as.integer(period), length(li$fits)))
  f <- li$fits[[period]]
  t_h <- (seq_along(f$lids) - 1) * f$bin_min / 60

  ggplot2::ggplot(data.frame(t = t_h, lids = f$lids, fitted = f$fitted),
                  ggplot2::aes(.data$t)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$lids), colour = .circ_color("blue"),
                       linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted), colour = .circ_color("orange"),
                       linewidth = 0.9) +
    ggplot2::labs(title = "LIDS ultradian cycle",
                  subtitle = sprintf("period %.0f min, MRI %.2f (orange = fitted cosine)",
                                     f$period_min, f$MRI),
                  x = "Hours into sleep period", y = "LIDS (smoothed)") +
    .circ_theme()
}


#' Rest-Detector Comparison Strip
#'
#' Runs the four rest/sleep detectors on one recording and stacks their detected
#' rest bands over the activity series on a shared time axis, so the user can see
#' that the differing bout counts reflect different questions (main night vs every
#' bout vs latent state), not contradictory accuracy. Returns a \code{ggplot}
#' object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{sleep.changepoints}}, \code{\link{rest.periods}},
#'   \code{\link{rest.crespo}}, \code{\link{rest.hmm}}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' plot_rest_comparison(ifelse(h >= 23 | h < 7, 5, 300), ts)
#' }
#'
#' @export
plot_rest_comparison <- function(counts, timestamps) {
  empty <- function() .circ_empty_plot("Insufficient data for rest detection",
                                        title = "Rest detectors compared")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps)
  if (sum(ok) < 64L) return(empty())
  o <- order(timestamps[ok]); tt <- timestamps[ok][o]; xx <- x[ok][o]; n <- length(tt)

  inband <- function(df, on, off) {
    r <- logical(n)
    if (!is.null(df) && nrow(df))
      for (i in seq_len(nrow(df))) r[tt >= df[[on]][i] & tt <= df[[off]][i]] <- TRUE
    r
  }
  cp <- tryCatch(sleep.changepoints(counts, timestamps), error = function(e) NULL)
  rp <- tryCatch(rest.periods(counts, timestamps), error = function(e) NULL)
  rc <- tryCatch(rest.crespo(counts, timestamps), error = function(e) NULL)
  hm <- tryCatch(rest.hmm(counts, timestamps, seed = 1L), error = function(e) NULL)

  rest <- list(
    `sleep.changepoints` = if (!is.null(cp) && !isTRUE(cp$insufficient))
      inband(cp$sleep_episodes, "sleep_onset", "wake_onset") else logical(n),
    `rest.periods` = if (!is.null(rp) && !isTRUE(rp$insufficient))
      inband(rp$rest_periods, "onset", "offset") else logical(n),
    `rest.crespo` = if (!is.null(rc) && !isTRUE(rc$insufficient))
      inband(rc$rest_periods, "onset", "offset") else logical(n),
    `rest.hmm` = if (!is.null(hm) && !isTRUE(hm$insufficient) && length(hm$sleep_state) == n)
      hm$sleep_state == "S" else logical(n))

  runs <- function(v, label) {
    if (!any(v)) return(NULL)
    r <- rle(v); ends <- cumsum(r$lengths); starts <- c(1L, utils::head(ends, -1L) + 1L)
    k <- which(r$values)
    data.frame(track = label, xmin = tt[starts[k]], xmax = tt[ends[k]])
  }
  bands <- do.call(rbind, Map(runs, rest, names(rest)))
  levs <- c("activity", names(rest))
  act <- data.frame(track = factor("activity", levs), time = tt, value = xx)
  if (!is.null(bands)) bands$track <- factor(bands$track, levs)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = act, ggplot2::aes(.data$time, .data$value),
                       colour = "grey45", linewidth = 0.3)
  if (!is.null(bands))
    p <- p + ggplot2::geom_rect(data = bands,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = 1),
      fill = .circ_color("blue"))
  p +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$track), scales = "free_y",
                        switch = "y", drop = FALSE) +
    ggplot2::labs(title = "Rest detectors compared",
                  subtitle = "navy = detected rest; the count differs because the question does",
                  x = "Time", y = NULL) +
    .circ_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   strip.text.y.left = ggplot2::element_text(angle = 0))
}
