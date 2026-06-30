#' Circular Phase Rose
#'
#' Plots daily phase markers (acrophases, onsets, L5/M10 times) as a rose diagram
#' around the 24-hour clock, with the mean resultant vector drawn from the centre
#' and the Rayleigh and Hermans-Rasson results annotated. This is the companion
#' figure to \code{\link{phase.concentration}}. Sector area, not radius, encodes
#' the count (the radius uses a square-root scale) so a wide sector is not read as
#' a large one. Returns a \code{ggplot} object and never errors.
#'
#' @param times_h Numeric vector of clock times (hours), one per day.
#' @param period Period the times wrap on, in hours (default 24).
#' @param binwidth Sector width in hours (default 1).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{fisher1993}{actiRhythm}
#'
#' \insertRef{landler2019}{actiRhythm}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' plot_phase_rose(23 + stats::rnorm(30, 0, 1))
#' }
#'
#' @export
plot_phase_rose <- function(times_h, period = 24, binwidth = 1) {
  empty <- function() .circ_empty_plot("Need at least 3 phase markers",
                                        title = "Phase rose")
  h <- suppressWarnings(as.numeric(times_h)) %% period
  h <- h[is.finite(h)]
  if (length(h) < 3L) return(empty())

  pc <- phase.concentration(h, period)
  brk <- seq(0, period, by = binwidth)
  ymax <- max(table(cut(h, breaks = brk, include.lowest = TRUE)))
  sub <- sprintf("R = %.2f at %04.1f h.  Rayleigh p = %.3f, Hermans-Rasson p = %.3f",
                 pc$R, pc$mean_direction_h, pc$rayleigh_p, pc$hr_p)

  ggplot2::ggplot(data.frame(h = h), ggplot2::aes(.data$h)) +
    ggplot2::geom_histogram(breaks = brk, fill = .circ_color("blue"),
                            colour = "white", linewidth = 0.2) +
    ggplot2::geom_segment(x = pc$mean_direction_h, xend = pc$mean_direction_h,
                          y = 0, yend = pc$R^2 * ymax, colour = .circ_color("orange"),
                          linewidth = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "inches"))) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::scale_x_continuous(limits = c(0, period),
                                breaks = seq(0, period - binwidth, by = period / 6)) +
    ggplot2::scale_y_sqrt() +
    ggplot2::labs(title = "Phase markers around the clock",
                  subtitle = sub, x = "Clock hour", y = NULL) +
    .circ_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}


#' Sleep / Wake Change-Point Track
#'
#' Plots the activity series with the per-night sleep-onset and wake-onset change
#' points from \code{\link{sleep.changepoints}} marked and the detected sleep
#' episodes shaded, so each night's rest timing is visible against the raw counts.
#' Returns a \code{ggplot} object and never errors.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param ... Passed to \code{\link{sleep.changepoints}} (e.g. \code{thr}).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{chensun2024}{actiRhythm}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' plot_changepoints(ifelse(h >= 8 & h < 23, 300, 5), ts)
#' }
#'
#' @export
plot_changepoints <- function(counts, timestamps, ...) {
  empty <- function() .circ_empty_plot("Insufficient data for change-point detection",
                                        title = "Sleep / wake change points")
  cp <- tryCatch(sleep.changepoints(counts, timestamps, ...), error = function(e) NULL)
  if (is.null(cp) || isTRUE(cp$insufficient)) return(empty())

  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps)
  series <- data.frame(time = timestamps[ok], counts = x[ok])
  series <- series[order(series$time), ]

  p <- ggplot2::ggplot(series, ggplot2::aes(.data$time, .data$counts))
  if (nrow(cp$sleep_episodes))
    p <- p + ggplot2::geom_rect(data = cp$sleep_episodes, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$sleep_onset, xmax = .data$wake_onset, ymin = -Inf, ymax = Inf),
      fill = .circ_color("blue"), alpha = 0.15)
  p +
    ggplot2::geom_line(colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_vline(data = cp$changepoints, inherit.aes = FALSE,
      ggplot2::aes(xintercept = .data$time, colour = .data$type), linewidth = 0.5) +
    ggplot2::scale_colour_manual(name = NULL,
      values = c(`sleep onset` = .circ_color("blue"), `wake onset` = .circ_color("orange"))) +
    ggplot2::labs(title = "Sleep / wake change points",
                  subtitle = sprintf("%d episodes; mean rest %.1f h", cp$n_episodes,
                                     cp$mean_sleep_duration),
                  x = "Time", y = "Activity") +
    .circ_theme()
}
