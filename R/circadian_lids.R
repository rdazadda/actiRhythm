#' Locomotor Inactivity During Sleep (LIDS)
#'
#' Quantifies the ultradian (sleep-cycle) oscillation of inactivity during sleep,
#' following Winnebeck et al. (2018). Within each sleep period the activity is
#' transformed to LIDS (\code{100 / (activity + 1)}), smoothed with a centered
#' 30-minute moving average, and fit with an ultradian cosine over a scan of
#' candidate periods; the best period is the one maximising the Munich
#' Rhythmicity Index (MRI).
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct (or parseable) timestamps, one per count.
#' @param sleep_periods Data frame with \code{in_bed_time} and \code{out_bed_time}
#'   (same schema as \code{\link{social.jet.lag}}); one LIDS fit per period.
#' @param epoch_length Epoch length in seconds (default 60).
#' @param smooth_minutes Centered moving-average window in minutes (default 30).
#' @param period_min,period_max,period_step Period scan grid in minutes
#'   (defaults 30, 180, 5).
#'
#' @return An object of class \code{actiRhythm_lids}: a per-period data frame
#'   (period in minutes, MRI, amplitude, Pearson r, ...) and the mean period and
#'   MRI across periods.
#'
#' @references
#' Winnebeck EC, Fischer D, Leise T, Roenneberg T (2018). Dynamics and ultradian
#' structure of human sleep in real life. Current Biology, 28(1):49-59.
#'
#' @export
lids <- function(counts, timestamps, sleep_periods, epoch_length = 60,
                 smooth_minutes = 30, period_min = 30, period_max = 180,
                 period_step = 5) {
  if (!all(c("in_bed_time", "out_bed_time") %in% names(sleep_periods))) {
    stop("'sleep_periods' needs 'in_bed_time' and 'out_bed_time' columns")
  }
  to_ct <- function(z) if (inherits(z, "POSIXct")) z else as.POSIXct(as.character(z), tz = "UTC")
  ts <- to_ct(timestamps)
  ib <- to_ct(sleep_periods$in_bed_time)
  ob <- to_ct(sleep_periods$out_bed_time)
  grid_epochs <- seq(period_min, period_max, by = period_step) * 60 / epoch_length

  rows <- lapply(seq_len(nrow(sleep_periods)), function(i) {
    sel <- ts >= ib[i] & ts < ob[i]
    if (sum(sel) < 10L) return(NULL)
    fit <- .lids_fit(counts[sel], epoch_length, smooth_minutes, grid_epochs)
    if (is.null(fit)) return(NULL)
    data.frame(period = i,
               period_min = fit$period_epochs * epoch_length / 60,
               MRI = fit$MRI, amplitude = fit$amplitude,
               pearson_r = fit$pearson_r, offset = fit$offset,
               n_epochs = fit$n_epochs, stringsAsFactors = FALSE)
  })
  per <- do.call(rbind, rows)
  if (is.null(per)) {
    return(structure(list(periods = NULL, mean_period_min = NA_real_,
      mean_MRI = NA_real_, n_periods = 0L, epoch_length = epoch_length,
      insufficient = TRUE), class = "actiRhythm_lids"))
  }

  structure(list(
    periods = per,
    mean_period_min = mean(per$period_min),
    mean_MRI = mean(per$MRI),
    n_periods = nrow(per),
    epoch_length = epoch_length,
    insufficient = FALSE
  ), class = "actiRhythm_lids")
}


# Fit the LIDS ultradian cosine for one sleep period: transform, smooth, then
# scan candidate periods and keep the one with the highest MRI. At each fixed
# period the cosine is linear in (offset, cos, sin), so a plain WLS fit suffices.
.lids_fit <- function(activity, epoch_length, smooth_minutes, grid_epochs) {
  lids_raw <- 100 / (as.numeric(activity) + 1)
  win  <- max(1L, round(smooth_minutes * 60 / epoch_length))
  y    <- .centered_ma(lids_raw, win)
  y    <- y[is.finite(y)]
  n    <- length(y)
  if (n < 10L) return(NULL)
  x <- seq_len(n) - 1L

  best <- NULL
  for (Te in grid_epochs) {
    if (Te < 2 || Te > n) next
    w <- 2 * pi / Te
    X <- cbind(1, cos(w * x), sin(w * x))
    fit <- stats::lm.fit(X, y)
    if (any(!is.finite(fit$coefficients))) next
    a <- fit$coefficients[2]; b <- fit$coefficients[3]
    A <- sqrt(a^2 + b^2)
    yhat <- as.numeric(X %*% fit$coefficients)
    r <- suppressWarnings(stats::cor(y, yhat))
    if (!is.finite(r)) next
    MRI <- 2 * A * r
    if (is.null(best) || MRI > best$MRI) {
      best <- list(period_epochs = Te, amplitude = A, offset = fit$coefficients[1],
                   pearson_r = r, MRI = MRI, n_epochs = n)
    }
  }
  best
}

# Centered moving average with min_periods = 1 (edges average available points).
.centered_ma <- function(x, win) {
  if (win <= 1L) return(x)
  n <- length(x)
  half <- floor(win / 2)
  vapply(seq_len(n), function(i) {
    mean(x[max(1L, i - half):min(n, i + half)], na.rm = TRUE)
  }, numeric(1))
}


#' @export
print.actiRhythm_lids <- function(x, ...) {
  cat("Locomotor Inactivity During Sleep (LIDS)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Sleep periods:    %d\n", x$n_periods))
  cat(sprintf("  Mean LIDS period: %.1f min\n", x$mean_period_min))
  cat(sprintf("  Mean MRI:         %.3f\n", x$mean_MRI))
  invisible(x)
}
