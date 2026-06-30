#' Average-Day Profile with the L5 and M10 Windows
#'
#' The averaged daily activity profile (mean across days with a one-standard-
#' deviation band), with the least-active 5-hour window L5 and the most-active
#' 10-hour window M10 marked on a window track beneath the curve. Returns a
#' \code{ggplot} object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param bin_min Profile resolution in minutes (default 30).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{vansomeren1999}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 7 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' plot_profile(pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24)), ts)
#'
#' @export
plot_profile <- function(counts, timestamps, bin_min = 30) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_profile().")
  empty <- function() .circ_empty_plot("Insufficient data for the average-day profile",
                                        title = "Average day")
  if (length(counts) != length(timestamps)) return(empty())
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps)
  x <- x[ok]; timestamps <- timestamps[ok]
  if (length(x) < 48L) return(empty())

  nslot <- 1440L %/% bin_min
  lt  <- as.POSIXlt(timestamps)
  sec <- lt$hour * 3600 + lt$min * 60 + lt$sec
  slot <- factor(pmin(sec %/% (bin_min * 60), nslot - 1L), levels = 0:(nslot - 1L))
  prof <- data.frame(
    hour = (0:(nslot - 1L)) * (bin_min / 60) + bin_min / 120,
    mean = as.numeric(tapply(x, slot, mean, na.rm = TRUE)),
    sd   = as.numeric(tapply(x, slot, stats::sd, na.rm = TRUE)))
  prof$sd[is.na(prof$sd)] <- 0
  if (all(is.na(prof$mean))) return(empty())
  prof$lo <- pmax(prof$mean - prof$sd, 0)
  prof$hi <- prof$mean + prof$sd

  blue <- .circ_color("blue"); orange <- .circ_color("orange")
  ytop <- max(prof$hi, na.rm = TRUE)
  track <- -0.06 * ytop
  ext <- tryCatch(activity.extrema(x, timestamps, windows = c(5, 10))$table,
                  error = function(e) NULL)
  segs <- NULL; labs <- NULL
  if (!is.null(ext)) {
    seg_of <- function(onset, width, label, col) {
      o <- onset %% 24; e <- o + width
      pieces <- if (e <= 24) data.frame(x = o, xend = e)
                else data.frame(x = c(o, 0), xend = c(24, e - 24))
      pieces$col <- col
      list(seg = pieces,
           lab = data.frame(x = (pieces$x[1] + pieces$xend[1]) / 2, label = label, col = col))
    }
    l5  <- ext[ext$window_h == 5, ]; m10 <- ext[ext$window_h == 10, ]
    if (nrow(l5))  { s <- seg_of(l5$L_onset_h, 5,  "L5",  blue);   segs <- rbind(segs, s$seg); labs <- rbind(labs, s$lab) }
    if (nrow(m10)) { s <- seg_of(m10$M_onset_h, 10, "M10", orange); segs <- rbind(segs, s$seg); labs <- rbind(labs, s$lab) }
  }

  p <- ggplot2::ggplot(prof, ggplot2::aes(.data$hour, .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         fill = blue, alpha = 0.13) +
    ggplot2::geom_area(fill = blue, alpha = 0.18) +
    ggplot2::geom_line(color = blue, linewidth = 1)
  if (!is.null(segs)) {
    p <- p +
      ggplot2::geom_segment(data = segs, inherit.aes = FALSE,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = track, yend = track, color = .data$col),
        linewidth = 3, lineend = "round") +
      ggplot2::geom_text(data = labs, inherit.aes = FALSE,
        ggplot2::aes(x = .data$x, y = track, label = .data$label, color = .data$col),
        vjust = -0.8, fontface = "bold", size = 3.6, show.legend = FALSE) +
      ggplot2::scale_color_identity()
  }
  p +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0.01, 0)) +
    ggplot2::labs(title = "Average day", subtitle = "Mean activity by time of day, with the L5 and M10 windows",
                  x = "Hour of day", y = "Activity") +
    .circ_theme()
}


#' Dichotomy Index Plot
#'
#' Cumulative distributions of the rest-span and active-span activity counts on a
#' log scale, marking the active-span median and the dichotomy index I<O (the
#' share of rest-span counts below that median). Returns a \code{ggplot} object
#' and never errors.
#'
#' @param counts Numeric activity vector.
#' @param rest Logical (TRUE = rest), or a character state vector as accepted by
#'   \code{\link{dichotomy.index}}. Same length as \code{counts}.
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{mormont2000}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' plot_dichotomy(ifelse(h >= 23 | h < 7, 5, 300), rest = h >= 23 | h < 7)
#'
#' @export
plot_dichotomy <- function(counts, rest) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_dichotomy().")
  empty <- function() .circ_empty_plot("Insufficient data for the dichotomy index",
                                        title = "Dichotomy index")
  di <- tryCatch(dichotomy.index(counts, rest), error = function(e) NULL)
  if (is.null(di) || !is.finite(di$IO)) return(empty())

  x <- suppressWarnings(as.numeric(counts))
  r <- if (is.logical(rest)) rest
       else tolower(as.character(rest)) %in% c("r", "s", "sleep", "rest", "true", "1")
  ok <- is.finite(x) & !is.na(r)
  d <- data.frame(count = x[ok], span = ifelse(r[ok], "Rest", "Active"))
  blue <- .circ_color("blue"); orange <- .circ_color("orange")
  cols <- c(Rest = blue, Active = orange)
  med <- di$active_median

  ggplot2::ggplot(d, ggplot2::aes(.data$count, color = .data$span)) +
    ggplot2::stat_ecdf(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = med, linetype = "dashed", color = "grey55") +
    ggplot2::annotate("point", x = max(med, 0), y = di$IO / 100, color = blue, size = 3) +
    ggplot2::annotate("text", x = max(med, 0), y = di$IO / 100,
                      label = sprintf("  I<O = %.0f%%", di$IO),
                      hjust = 0, vjust = -0.6, fontface = "bold", color = blue) +
    ggplot2::scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(v * 100, "%"), limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::labs(title = "Dichotomy index (I<O)",
                  subtitle = "Cumulative share of counts; the dashed line is the active-span median",
                  x = "Activity count (log)", y = "Cumulative share") +
    .circ_theme()
}


#' Rest-Activity Transition Curves
#'
#' The rest-to-active and active-to-rest transition hazards against bout length,
#' with a LOWESS fit and the sustained rates kRA and kAR marked. Returns a
#' \code{ggplot} object and never errors.
#'
#' @param counts Numeric activity vector.
#' @param threshold Counts at or above which an epoch is active (default 1).
#' @param frac LOWESS span (default 0.3).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{lim2011}{actiRhythm}
#'
#' @examples
#' set.seed(1)
#' plot_transitions(as.integer(stats::runif(8000) < 0.1) * 100)
#'
#' @export
plot_transitions <- function(counts, threshold = 1, frac = 0.3) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_transitions().")
  empty <- function() .circ_empty_plot("Insufficient data for transition curves",
                                        title = "Transition curves")
  st <- tryCatch(state.transitions(counts, threshold = threshold, frac = frac),
                 error = function(e) NULL)
  if (is.null(st) || isTRUE(st$insufficient) ||
      (is.null(st$rest_curve) && is.null(st$act_curve))) return(empty())

  lev <- c("Rest to active", "Active to rest")
  parts <- list()
  if (!is.null(st$rest_curve)) { rc <- st$rest_curve; rc$type <- lev[1]; parts[[length(parts) + 1L]] <- rc }
  if (!is.null(st$act_curve))  { ac <- st$act_curve;  ac$type <- lev[2]; parts[[length(parts) + 1L]] <- ac }
  d <- do.call(rbind, parts); d$type <- factor(d$type, levels = lev)
  smooth <- do.call(rbind, lapply(split(d, d$type), function(g) {
    if (nrow(g) < 3L) return(NULL)
    lo <- stats::lowess(g$lag, g$prob, f = frac)
    data.frame(lag = lo$x, prob = lo$y, type = g$type[1])
  }))
  klab <- data.frame(type = factor(lev, levels = lev), k = c(st$kRA, st$kAR),
                     label = sprintf("k = %.3f", c(st$kRA, st$kAR)))
  klab <- klab[klab$type %in% unique(d$type), ]
  blue <- .circ_color("blue"); orange <- .circ_color("orange")

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$lag, .data$prob, color = .data$type)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$weight), alpha = 0.45)
  if (!is.null(smooth))
    p <- p + ggplot2::geom_line(data = smooth, linewidth = 1)
  p +
    ggplot2::geom_hline(data = klab, ggplot2::aes(yintercept = .data$k),
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_text(data = klab, ggplot2::aes(x = Inf, y = .data$k, label = .data$label),
                       hjust = 1.1, vjust = -0.5, color = "grey35", size = 3.4,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ .data$type, scales = "free") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = stats::setNames(c(blue, orange), lev)) +
    ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
    ggplot2::labs(title = "Rest-activity transition hazards",
                  subtitle = "Per-epoch probability of ending a bout, by bout length",
                  x = "Bout length (epochs, log)", y = "Transition probability") +
    .circ_theme() +
    ggplot2::theme(legend.position = "none")
}


#' Multiscale IS and IV Profiles
#'
#' Interdaily stability and intradaily variability recomputed across epoch
#' lengths, with the averaged ISm and IVm marked. Returns a \code{ggplot} object
#' and never errors.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' \insertRef{goncalves2014}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' plot_multiscale(pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24)), ts)
#'
#' @export
plot_multiscale <- function(counts, timestamps) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_multiscale().")
  empty <- function() .circ_empty_plot("Insufficient data for the multiscale profiles",
                                        title = "Multiscale IS and IV")
  iv <- tryCatch(intradaily.variability.multiscale(counts, timestamps), error = function(e) NULL)
  is <- tryCatch(circadian.is.multiscale(counts, timestamps), error = function(e) NULL)
  if (is.null(iv) || is.null(is)) return(empty())

  lev <- c("IS (synchronisation)", "IV (fragmentation)")
  d <- rbind(
    data.frame(bin = is$table$bin_minutes, value = is$table$IS, metric = lev[1]),
    data.frame(bin = iv$table$bin_minutes, value = iv$table$IV, metric = lev[2]))
  d <- d[is.finite(d$value), ]; d$metric <- factor(d$metric, levels = lev)
  if (!nrow(d)) return(empty())
  means <- data.frame(metric = factor(lev, levels = lev), m = c(is$ISm, iv$IVm),
                      label = sprintf(c("ISm = %.2f", "IVm = %.2f"), c(is$ISm, iv$IVm)))
  blue <- .circ_color("blue"); orange <- .circ_color("orange")

  ggplot2::ggplot(d, ggplot2::aes(.data$bin, .data$value, color = .data$metric)) +
    ggplot2::geom_hline(data = means, ggplot2::aes(yintercept = .data$m),
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_text(data = means, ggplot2::aes(x = Inf, y = .data$m, label = .data$label),
                       hjust = 1.1, vjust = -0.5, color = "grey35", size = 3.4,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y") +
    ggplot2::scale_color_manual(values = stats::setNames(c(blue, orange), lev)) +
    ggplot2::labs(title = "Multiscale IS and IV",
                  subtitle = "Recomputed across epoch lengths and averaged",
                  x = "Epoch length (minutes)", y = NULL) +
    .circ_theme() +
    ggplot2::theme(legend.position = "none")
}
