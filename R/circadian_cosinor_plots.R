#' Cosinor Amplitude-Acrophase Confidence Ellipse
#'
#' Draws the joint confidence ellipse of the cosinor cosine and sine coefficients
#' on the amplitude-acrophase plane, over clock-hour spokes and amplitude rings.
#' The estimate is the vector from the pole, and the rhythm is detectable when the
#' ellipse excludes the pole (Bingham et al. 1982). Returns a \code{ggplot} object
#' and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Cosinor period in hours (default 24).
#' @param level Confidence level for the ellipse (default 0.95).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{bingham1982}{actiRhythm}
#'
#' @examples
#' set.seed(1)
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24) + rnorm(length(h), 0, 25))
#' plot_cosinor_ellipse(counts, ts)
#'
#' @export
plot_cosinor_ellipse <- function(counts, timestamps, period = 24, level = 0.95) {
  empty <- function() .circ_empty_plot("Confidence ellipse unavailable",
                                        title = "Cosinor confidence ellipse")
  cos <- tryCatch(cosinor.analysis(counts, timestamps, period = period),
                  error = function(e) NULL)
  if (is.null(cos)) return(empty())
  ell <- tryCatch(cosinor.confidence.ellipse(cos, level = level),
                  error = function(e) NULL)
  if (is.null(ell) || anyNA(ell$center) || anyNA(ell$ellipse$x) || nrow(ell$ellipse) < 3L)
    return(empty())

  cx <- unname(ell$center[1]); cy <- unname(ell$center[2])
  amp  <- sqrt(cx^2 + cy^2)
  acro <- (atan2(cy, cx) * period / (2 * pi)) %% period
  hit  <- isTRUE(ell$rhythm_detected)
  col  <- .circ_color(if (hit) "blue" else "orange")

  rmax  <- max(sqrt(ell$ellipse$x^2 + ell$ellipse$y^2), amp) * 1.15
  th    <- seq(0, 2 * pi, length.out = 100)
  rings <- pretty(c(0, rmax), n = 4); rings <- rings[rings > 0 & rings < rmax]
  ringd <- if (length(rings))
    do.call(rbind, lapply(rings, function(r) data.frame(x = r * cos(th), y = r * sin(th), r = r)))
  else data.frame(x = numeric(0), y = numeric(0), r = numeric(0))
  hrs    <- seq(0, period - 1e-9, by = period / 8)
  ang    <- hrs * 2 * pi / period
  spokes <- data.frame(x = rmax * cos(ang), y = rmax * sin(ang),
                       lab = sprintf("%02d:00", as.integer(round(hrs))))
  vec  <- data.frame(x = 0, y = 0, xend = cx, yend = cy)
  pole <- data.frame(x = 0, y = 0)

  ggplot2::ggplot() +
    ggplot2::geom_path(data = ringd,
      ggplot2::aes(.data$x, .data$y, group = .data$r),
      colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_segment(data = spokes,
      ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
      colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_text(data = spokes,
      ggplot2::aes(1.07 * .data$x, 1.07 * .data$y, label = .data$lab),
      size = 3, colour = "grey45") +
    ggplot2::geom_polygon(data = ell$ellipse,
      ggplot2::aes(.data$x, .data$y),
      fill = col, alpha = 0.18, colour = col, linewidth = 0.7) +
    ggplot2::geom_segment(data = vec,
      ggplot2::aes(.data$x, .data$y, xend = .data$xend, yend = .data$yend),
      colour = col, linewidth = 0.9,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "inches"))) +
    ggplot2::geom_point(data = pole, ggplot2::aes(.data$x, .data$y),
      shape = 3, size = 3.5, stroke = 0.9, colour = "grey25") +
    ggplot2::coord_fixed(xlim = c(-rmax, rmax) * 1.18, ylim = c(-rmax, rmax) * 1.18) +
    ggplot2::labs(
      title = "Cosinor confidence ellipse",
      subtitle = sprintf("amplitude %.1f, acrophase %.1f h, %s",
        amp, acro, if (hit) "rhythm detected" else "not detected"),
      x = expression(beta[1] == A * cos(phi)),
      y = expression(beta[2] == A * sin(phi))) +
    .circ_theme()
}


#' Single Cosine vs Multicomponent Cosinor Fit
#'
#' Overlays the averaged daily profile with the single 24-hour cosine and the
#' selected multi-harmonic cosinor (Cornelissen 2014), showing the structure a
#' single symmetric cosine cannot follow. The number of harmonics and the fit are
#' taken from \code{\link{cosinor.multicomponent}}; the curves are refit on the
#' hour-of-day profile so they align with the plotted points. Returns a
#' \code{ggplot} object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Fundamental period in hours (default 24).
#' @param max_harmonics Largest number of harmonics to consider (default 3).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{cornelissen2014}{actiRhythm}
#'
#' @examples
#' set.seed(2)
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' siesta <- pmax(0, 120 + 70 * cos(2 * pi * (h - 14) / 24) +
#'   55 * cos(2 * pi * 2 * (h - 14) / 24) + rnorm(length(ts), 0, 12))
#' plot_multicomponent(siesta, ts)
#'
#' @export
plot_multicomponent <- function(counts, timestamps, period = 24, max_harmonics = 3) {
  empty <- function() .circ_empty_plot("Insufficient data for the multicomponent fit",
                                        title = "Multicomponent cosinor")
  if (length(counts) != length(timestamps)) return(empty())
  prof <- .circ_hourly_profile(counts, timestamps)
  if (nrow(prof) < 4L) return(empty())
  mc <- tryCatch(cosinor.multicomponent(counts, timestamps, period = period,
                                        max_harmonics = max_harmonics),
                 error = function(e) NULL)
  if (is.null(mc) || isTRUE(mc$insufficient)) return(empty())
  K <- mc$n_harmonics

  ph <- prof$hour; py <- prof$activity
  grid <- seq(0, period, length.out = 481)
  fit_curve <- function(k) {
    cols_at <- function(x) cbind(1, do.call(cbind, lapply(seq_len(k), function(j)
      cbind(cos(2 * pi * j * x / period), sin(2 * pi * j * x / period)))))
    b <- stats::lm.fit(cols_at(ph), py)$coefficients
    data.frame(hour = grid, fit = as.numeric(cols_at(grid) %*% b))
  }
  lab_multi <- sprintf("%d harmonics", K)
  single <- fit_curve(1L); single$model <- "Single cosine"
  curves <- if (K > 1L) {
    multi <- fit_curve(K); multi$model <- lab_multi
    rbind(single, multi)
  } else single
  curves$model <- factor(curves$model, levels = c("Single cosine", lab_multi))
  cols <- c("Single cosine" = .circ_color("orange"), .circ_color("blue"))
  names(cols)[2] <- lab_multi

  sub <- if (K > 1L)
    sprintf("%d harmonics capture what one cosine misses (R-squared %.2f)", K, mc$r_squared)
  else sprintf("a single cosine is adequate here (R-squared %.2f)", mc$r_squared)

  ggplot2::ggplot() +
    ggplot2::geom_point(data = prof, ggplot2::aes(.data$hour, .data$activity),
                        colour = "grey55", size = 1.4) +
    ggplot2::geom_line(data = curves,
                       ggplot2::aes(.data$hour, .data$fit, colour = .data$model),
                       linewidth = 1) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, period, by = 6), limits = c(0, period)) +
    ggplot2::labs(title = "Single cosine vs multicomponent fit", subtitle = sub,
                  x = "Hour of day", y = "Activity") +
    .circ_theme() +
    ggplot2::theme(legend.position = "top")
}


#' Annotated Cosinor Parameter Schematic
#'
#' A teaching figure: one smooth cosine annotating the MESOR (midline), the
#' amplitude A and double amplitude 2A, the acrophase (clock time of the peak),
#' and the period. Useful as a legend for the cosinor parameters. Returns a
#' \code{ggplot} object and never errors.
#'
#' @param mesor,amplitude Midline and amplitude of the illustrated cosine.
#' @param acrophase Clock time of the peak, in hours.
#' @param period Period in hours (default 24).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot_cosinor_schematic()
#'
#' @export
plot_cosinor_schematic <- function(mesor = 100, amplitude = 50, acrophase = 16,
                                   period = 24) {
  t <- seq(0, period, length.out = 400)
  curve <- data.frame(t = t, y = mesor + amplitude * cos(2 * pi * (t - acrophase) / period))
  col <- .circ_color("blue")
  peak <- mesor + amplitude; trough <- mesor - amplitude
  arr <- ggplot2::arrow(length = ggplot2::unit(0.06, "inches"), ends = "both")

  ggplot2::ggplot(curve, ggplot2::aes(.data$t, .data$y)) +
    ggplot2::geom_hline(yintercept = mesor, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_line(colour = col, linewidth = 1) +
    ggplot2::annotate("point", x = acrophase, y = peak, colour = col, size = 2.6) +
    ggplot2::annotate("segment", x = acrophase, xend = acrophase, y = mesor, yend = peak,
                      colour = "grey25", arrow = arr) +
    ggplot2::annotate("text", x = acrophase + 0.4, y = (mesor + peak) / 2,
                      label = "A", hjust = 0, size = 4) +
    ggplot2::annotate("segment", x = acrophase, xend = period, y = peak, yend = peak,
                      linetype = "dotted", colour = "grey75") +
    ggplot2::annotate("segment", x = (acrophase + period / 2) %% period, xend = period,
                      y = trough, yend = trough, linetype = "dotted", colour = "grey75") +
    ggplot2::annotate("segment", x = period * 0.95, xend = period * 0.95, y = trough, yend = peak,
                      colour = "grey25", arrow = arr) +
    ggplot2::annotate("text", x = period * 0.95 - 0.4, y = mesor, label = "2A",
                      hjust = 1, size = 4) +
    ggplot2::annotate("segment", x = acrophase, xend = acrophase, y = trough - 0.18 * amplitude,
                      yend = mesor, linetype = "dotted", colour = "grey45") +
    ggplot2::annotate("text", x = acrophase, y = trough - 0.26 * amplitude, label = "acrophase",
                      hjust = 0.5, vjust = 1, size = 3.5) +
    ggplot2::annotate("text", x = 0.3, y = mesor + 0.06 * amplitude, label = "MESOR",
                      hjust = 0, vjust = 0, size = 3.5, colour = "grey35") +
    ggplot2::scale_x_continuous(breaks = seq(0, period, by = 6), limits = c(0, period)) +
    ggplot2::labs(title = "Cosinor parameters",
                  subtitle = sprintf("MESOR, amplitude A, double amplitude 2A, acrophase (period %gh)", period),
                  x = "Hour of day", y = "Activity") +
    .circ_theme()
}
