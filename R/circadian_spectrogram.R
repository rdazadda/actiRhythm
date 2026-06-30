#' Circadian Spectrogram (Period over Time)
#'
#' Slides a window across the recording and computes the chi-square
#' (Sokolove-Bushell) periodogram in each window, producing a period-by-time map
#' that shows how the dominant period and its strength drift across the recording
#' (non-stationarity, fragmentation, re-entrainment). A single global
#' periodogram or cosinor fit cannot show this.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param window_hours Sliding-window length in hours (default 72).
#' @param step_hours Step between successive windows in hours (default 6).
#' @param from,to Period search window in hours (default 18, 30).
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return An object of class \code{actiRhythm_spectrogram}: a long \code{data}
#'   frame (window centre time, period, power) and a \code{ggplot} heat map in
#'   \code{$plot}.
#'
#' @seealso \code{\link{chi.sq.periodogram}}, \code{\link{circadian.period}}
#'
#' @examples
#' \donttest{
#' t_hours <- seq(0, 8 * 24, by = 1 / 60)
#' ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
#' counts <- 100 + 80 * cos(2 * pi * t_hours / 24) + rnorm(length(t_hours), 0, 20)
#' circadian.spectrogram(counts, ts, step_hours = 24)$plot
#' }
#'
#' @export
circadian.spectrogram <- function(counts, timestamps, window_hours = 72,
                                  step_hours = 6, from = 18, to = 30,
                                  epoch_length = 60) {
  insuf <- function() structure(list(data = NULL,
    plot = .circ_empty_plot("Insufficient data", title = "Circadian Spectrogram"),
    window_hours = window_hours, step_hours = step_hours, n_windows = 0L,
    insufficient = TRUE), class = "actiRhythm_spectrogram")
  ts  <- as.POSIXct(timestamps)
  cnt <- as.numeric(counts)
  ok  <- is.finite(cnt) & !is.na(ts)
  ts  <- ts[ok]; cnt <- cnt[ok]
  if (length(cnt) < 10L) return(insuf())

  t0     <- min(ts)
  span_h <- as.numeric(difftime(max(ts), t0, units = "hours"))
  if (span_h < window_hours) return(insuf())

  starts <- seq(0, span_h - window_hours, by = step_hours)
  rows <- lapply(starts, function(s) {
    sel <- ts >= (t0 + s * 3600) & ts < (t0 + (s + window_hours) * 3600)
    pg <- tryCatch(chi.sq.periodogram(cnt[sel], ts[sel], from = from, to = to,
                                      epoch_length = epoch_length),
                   error = function(e) NULL)
    if (is.null(pg) || length(pg$scanned) == 0L) return(NULL)
    data.frame(center_time = t0 + (s + window_hours / 2) * 3600,
               window_start_h = s, period_h = pg$scanned, power = pg$Qp,
               rel = pg$Qp / pg$critical, stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, rows)
  if (is.null(long) || nrow(long) == 0L) return(insuf())

  peak <- do.call(rbind, by(long, long$center_time, function(d) d[which.max(d$rel), ]))
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$center_time, .data$period_h, fill = .data$rel)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Qp / threshold") +
    ggplot2::geom_line(data = peak, ggplot2::aes(.data$center_time, .data$period_h),
                       colour = "white", linewidth = 0.4, inherit.aes = FALSE) +
    ggplot2::geom_hline(yintercept = 24, linetype = "dashed", colour = "white", linewidth = 0.3) +
    ggplot2::labs(x = "Time", y = "Period (hours)",
                  title = "Circadian Spectrogram (period over time)") +
    .circ_theme()

  structure(list(data = long, plot = p, window_hours = window_hours,
                 step_hours = step_hours, n_windows = length(starts),
                 insufficient = FALSE),
            class = "actiRhythm_spectrogram")
}


#' @export
print.actiRhythm_spectrogram <- function(x, ...) {
  cat("Circadian Spectrogram\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Windows:      %d (%g h window, %g h step)\n",
              x$n_windows, x$window_hours, x$step_hours))
  cat(sprintf("  Period range: %g-%g h\n", min(x$data$period_h), max(x$data$period_h)))
  cat("  Use $plot for the heat map.\n")
  invisible(x)
}
