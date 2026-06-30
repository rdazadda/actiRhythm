#' Wavelet Power Scalogram
#'
#' Draws the Morlet wavelet power surface from \code{\link{circadian.wavelet}}:
#' time on the x-axis, period (hours, log scale) on the y, scale-rectified power as
#' the fill, with the cone of influence faded out, the per-time dominant-period
#' ridge traced, and the 24-hour reference marked. Returns a \code{ggplot} object
#' and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{torrence1998}{actiRhythm}
#'
#' \insertRef{leise2013}{actiRhythm}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 8 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' plot_wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#' }
#'
#' @export
plot_wavelet <- function(counts, timestamps, epoch_length = 60) {
  empty <- function() .circ_empty_plot("Insufficient data for the wavelet transform",
                                        title = "Wavelet scalogram")
  w <- tryCatch(circadian.wavelet(counts, timestamps, epoch_length = epoch_length),
                error = function(e) NULL)
  if (is.null(w) || isTRUE(w$insufficient) || !length(w$period_hours)) return(empty())

  np <- length(w$period_hours); nt <- length(w$times)
  surf <- data.frame(time = rep(w$times, each = np),
                     period = rep(w$period_hours, times = nt),
                     power = as.vector(w$power))
  pr   <- range(w$period_hours)
  coi  <- data.frame(time = w$times, period = pmin(pmax(w$coi_period_h, pr[1]), pr[2]))
  ridge <- data.frame(time = w$times, period = w$dominant_period)
  ridge <- ridge[is.finite(ridge$period), ]
  ymax <- pr[2]

  ggplot2::ggplot(surf, ggplot2::aes(.data$time, .data$period)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$power)) +
    ggplot2::scale_fill_viridis_c(name = "Power", option = "D") +
    ggplot2::geom_ribbon(data = coi,
      ggplot2::aes(x = .data$time, ymin = .data$period, ymax = ymax),
      inherit.aes = FALSE, fill = "white", alpha = 0.35) +
    ggplot2::geom_line(data = coi, ggplot2::aes(.data$time, .data$period),
                       colour = "white", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_line(data = ridge, ggplot2::aes(.data$time, .data$period),
                       colour = "white", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 24, linetype = "dotted", colour = "grey85", linewidth = 0.4) +
    ggplot2::scale_y_continuous(trans = "log2", breaks = c(3, 6, 12, 24, 48)) +
    ggplot2::labs(title = "Wavelet power scalogram",
                  subtitle = sprintf("dominant period %.1f h; faded region = cone of influence",
                                     w$peak_period_h),
                  x = "Time", y = "Period (hours, log)") +
    .circ_theme()
}


#' Empirical Mode Decomposition Stack
#'
#' Stacks the intrinsic mode functions and residual trend from
#' \code{\link{circadian.emd}}, finest at the top to the trend at the bottom, over
#' a shared time axis, with the circadian IMF highlighted. Returns a \code{ggplot}
#' object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param ensemble Ensemble size for EEMD (default 1 = plain EMD).
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{huang1998}{actiRhythm}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 8 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' plot_emd(100 + 60 * cos(2 * pi * th / 24) + 20 * cos(2 * pi * th / 8), ts, epoch_length = 600)
#' }
#'
#' @export
plot_emd <- function(counts, timestamps, ensemble = 1L, epoch_length = 60) {
  empty <- function() .circ_empty_plot("Insufficient data for EMD",
                                        title = "Empirical mode decomposition")
  e <- tryCatch(circadian.emd(counts, timestamps, ensemble = ensemble,
                              epoch_length = epoch_length),
                error = function(e) NULL)
  if (is.null(e) || isTRUE(e$insufficient) || !ncol(e$imfs)) return(empty())

  nimf <- ncol(e$imfs)
  mat  <- cbind(e$imfs, e$residual)
  labs <- c(paste0("IMF ", seq_len(nimf)), "trend")
  circ <- if (is.finite(e$circadian_imf)) paste0("IMF ", e$circadian_imf) else NA_character_
  long <- data.frame(time = rep(e$times, ncol(mat)),
                     value = as.vector(mat),
                     component = factor(rep(labs, each = length(e$times)), levels = labs))
  long$circadian <- !is.na(circ) & long$component == circ

  ggplot2::ggplot(long, ggplot2::aes(.data$time, .data$value, colour = .data$circadian)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = .circ_color("blue"),
                                            `TRUE` = .circ_color("orange")), guide = "none") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$component), scales = "free_y") +
    ggplot2::labs(title = "Empirical mode decomposition",
                  subtitle = if (!is.na(circ))
                    sprintf("%s (orange) is the circadian mode, period %.1f h", circ, e$circadian_period)
                  else "intrinsic mode functions, finest to trend",
                  x = "Time", y = NULL) +
    .circ_theme() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}


#' SSA w-Correlation Matrix
#'
#' Plots the weighted-correlation matrix of the leading SSA elementary components
#' from \code{\link{circadian.ssa}}. Bright off-diagonal 2x2 blocks are oscillatory
#' pairs (the grouping aid of Golyandina 2013); a single bright cell is the trend
#' and a diffuse block is noise. Returns a \code{ggplot} object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param window_hours SSA window length in hours (default 48).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{golyandina2013}{actiRhythm}
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 8 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' plot_ssa_wcor(100 + 60 * cos(2 * pi * th / 24), ts)
#' }
#'
#' @export
plot_ssa_wcor <- function(counts, timestamps, window_hours = 48) {
  empty <- function() .circ_empty_plot("Insufficient data for SSA",
                                        title = "SSA w-correlation")
  s <- tryCatch(circadian.ssa(counts, timestamps, window_hours = window_hours),
                error = function(e) NULL)
  if (is.null(s) || isTRUE(s$insufficient) || is.null(s$wcor) || !nrow(s$wcor)) return(empty())

  m <- abs(s$wcor); wc <- nrow(m)
  d <- data.frame(i = rep(seq_len(wc), times = wc),
                  j = rep(seq_len(wc), each = wc),
                  w = as.vector(m))

  ggplot2::ggplot(d, ggplot2::aes(.data$i, .data$j, fill = .data$w)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(name = expression("|" * rho^w * "|"),
                                 low = "white", high = "#1E3A5F", limits = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq_len(wc), expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_len(wc), expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "SSA w-correlation matrix",
                  subtitle = "bright 2x2 blocks = oscillatory pairs to group together",
                  x = "Elementary component", y = "Elementary component") +
    .circ_theme()
}
