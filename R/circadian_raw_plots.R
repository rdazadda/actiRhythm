#' z-Angle Sleep-Period-Time Detection
#'
#' Plots the z-angle of a raw recording with the HDCZA sleep-period-time window(s)
#' shaded and non-wear marked, the companion figure to \code{\link{rest.spt}}.
#' Computes the angle, non-wear, and SPT window from the raw input. Returns a
#' \code{ggplot} object and never errors.
#'
#' @param x A path to a raw file or a raw data frame (as in \code{\link{raw.metrics}}).
#' @param epoch_length Epoch length in seconds (default 5, matching \code{rest.spt}).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{vanhees2018}{actiRhythm}
#'
#' @examples
#' \donttest{
#' plot_spt(example_raw(days = 2))
#' }
#'
#' @export
plot_spt <- function(x, epoch_length = 5) {
  empty <- function() .circ_empty_plot("Insufficient data for SPT detection",
                                        title = "z-angle SPT window")
  m <- tryCatch(raw.metrics(x, epoch = epoch_length, metrics = "anglez"),
                error = function(e) NULL)
  if (is.null(m) || !nrow(m) || all(is.na(m$anglez))) return(empty())
  wear <- tryCatch(detect.nonwear.raw(x, epoch = epoch_length),
                   error = function(e) rep(TRUE, nrow(m)))
  spt  <- tryCatch(rest.spt(m$anglez, m$time, epoch_length = epoch_length, wear = wear),
                   error = function(e) NULL)
  series <- m[seq(1, nrow(m), by = max(1L, round(60 / epoch_length))), , drop = FALSE]

  runs <- function(v) {
    r <- rle(v); ends <- cumsum(r$lengths); starts <- c(1L, utils::head(ends, -1L) + 1L)
    k <- which(r$values)
    if (!length(k)) return(NULL)
    data.frame(xmin = m$time[starts[k]], xmax = m$time[ends[k]])
  }

  p <- ggplot2::ggplot()
  nw <- if (length(wear) == nrow(m)) runs(!wear) else NULL
  if (!is.null(nw))
    p <- p + ggplot2::geom_rect(data = nw, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -90, ymax = 90),
      fill = "grey80", alpha = 0.5)
  if (!is.null(spt) && nrow(spt))
    p <- p + ggplot2::geom_rect(data = spt, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$onset, xmax = .data$offset, ymin = -90, ymax = 90),
      fill = .circ_color("blue"), alpha = 0.22)
  p +
    ggplot2::geom_line(data = series, ggplot2::aes(.data$time, .data$anglez),
                       colour = "grey25", linewidth = 0.25) +
    ggplot2::scale_y_continuous(limits = c(-90, 90)) +
    ggplot2::labs(title = "z-angle sleep-period-time detection",
                  subtitle = sprintf("%d SPT window(s) shaded navy; grey = non-wear",
                                     if (is.null(spt)) 0L else nrow(spt)),
                  x = "Time", y = "z-angle (degrees)") +
    .circ_theme()
}


#' Raw-Acceleration Epoch Metrics Profile
#'
#' Plots the per-epoch ENMO, MAD, and z-angle from \code{\link{raw.metrics}} as a
#' faceted time series, a quality-control view of the gravity-preserving signals.
#' Returns a \code{ggplot} object and never errors.
#'
#' @param x A path to a raw file or a raw data frame (as in \code{\link{raw.metrics}}).
#' @param epoch_length Epoch length in seconds (default 60).
#' @param ... Passed to \code{\link{raw.metrics}} (e.g. \code{metrics}, \code{calibrate}).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{vanhees2013}{actiRhythm}
#'
#' \insertRef{vahaypya2015}{actiRhythm}
#'
#' \insertRef{vanhees2015}{actiRhythm}
#'
#' @examples
#' \donttest{
#' plot_raw_metrics(example_raw(days = 1))
#' }
#'
#' @export
plot_raw_metrics <- function(x, epoch_length = 60, ...) {
  empty <- function() .circ_empty_plot("Insufficient data for raw metrics",
                                        title = "Raw-acceleration metrics")
  m <- tryCatch(raw.metrics(x, epoch = epoch_length, ...), error = function(e) NULL)
  if (is.null(m) || !nrow(m)) return(empty())
  mets <- setdiff(names(m), "time")
  labs <- c(ENMO = "ENMO (mg)", MAD = "MAD (mg)", anglez = "z-angle (deg)")
  long <- do.call(rbind, lapply(mets, function(k) data.frame(
    time = m$time, metric = if (k %in% names(labs)) labs[[k]] else k, value = m[[k]])))
  long$metric <- factor(long$metric, levels = vapply(mets,
    function(k) if (k %in% names(labs)) labs[[k]] else k, character(1)))

  ggplot2::ggplot(long, ggplot2::aes(.data$time, .data$value)) +
    ggplot2::geom_line(colour = .circ_color("blue"), linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_y", ncol = 1) +
    ggplot2::labs(title = "Raw-acceleration epoch metrics", x = "Time", y = NULL) +
    .circ_theme()
}
