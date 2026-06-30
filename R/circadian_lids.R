#' Locomotor Inactivity During Sleep (LIDS)
#'
#' Quantifies the ultradian (sleep-cycle) oscillation of inactivity during sleep,
#' following Winnebeck et al. (2018). Within each sleep period the activity is
#' summed into 10-minute bins, transformed to LIDS (\code{100 / (activity + 1)}),
#' smoothed with a centered 30-minute moving average, and fit with an ultradian
#' cosine over a scan of candidate periods; the best period is the interior peak
#' that maximises the Munich Rhythmicity Index (MRI).
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
#' \insertRef{winnebeck2018}{actiRhythm}
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
  grid_min <- seq(period_min, period_max, by = period_step)

  fits <- Filter(Negate(is.null), lapply(seq_len(nrow(sleep_periods)), function(i) {
    sel <- ts >= ib[i] & ts < ob[i]
    if (sum(sel) < 10L) return(NULL)
    f <- .lids_fit(counts[sel], epoch_length, smooth_minutes, grid_min)
    if (!is.null(f)) f$period <- i
    f
  }))
  if (!length(fits)) {
    return(structure(list(periods = NULL, fits = list(), mean_period_min = NA_real_,
      mean_MRI = NA_real_, n_periods = 0L, epoch_length = epoch_length,
      insufficient = TRUE), class = "actiRhythm_lids"))
  }

  per <- do.call(rbind, lapply(fits, function(f) data.frame(
    period = f$period, period_min = f$period_min, MRI = f$MRI, amplitude = f$amplitude,
    pearson_r = f$pearson_r, offset = f$offset, n_epochs = f$n_epochs,
    stringsAsFactors = FALSE)))

  structure(list(
    periods = per, fits = fits,
    mean_period_min = mean(per$period_min),
    mean_MRI = mean(per$MRI),
    n_periods = nrow(per),
    epoch_length = epoch_length,
    insufficient = FALSE
  ), class = "actiRhythm_lids")
}


# Fit the LIDS ultradian cosine for one sleep period. Winnebeck (2018) sums
# activity into 10-min bins, transforms to LIDS = 100/(activity+1) per bin, then
# smooths; at each fixed period the cosine is linear in (offset, cos, sin), so a
# plain least-squares fit suffices. The period scan keeps the highest-MRI
# interior peak (not a grid-boundary maximum).
.lids_fit <- function(activity, epoch_length, smooth_minutes, grid_min) {
  bin_w   <- max(1L, round(600 / epoch_length))       # epochs per 10-min bin
  bin_min <- bin_w * epoch_length / 60                # minutes per bin
  a  <- as.numeric(activity)
  nb <- length(a) %/% bin_w
  if (nb < 10L) return(NULL)
  binned <- vapply(seq_len(nb), function(j)
    sum(a[((j - 1L) * bin_w + 1L):(j * bin_w)]), numeric(1))
  lids_raw <- 100 / (binned + 1)
  win <- max(1L, round(smooth_minutes / bin_min))      # smoothing window in bins
  y   <- .centered_ma(lids_raw, win)
  y   <- y[is.finite(y)]
  n   <- length(y)
  if (n < 10L) return(NULL)
  x <- seq_len(n) - 1L

  grid <- grid_min / bin_min                           # candidate periods in bins
  mri <- amp <- off <- rr <- rep(NA_real_, length(grid))
  for (k in seq_along(grid)) {
    Te <- grid[k]
    if (Te < 2 || Te > n) next
    w <- 2 * pi / Te
    X <- cbind(1, cos(w * x), sin(w * x))
    fit <- stats::lm.fit(X, y)
    if (any(!is.finite(fit$coefficients))) next
    A <- sqrt(fit$coefficients[2]^2 + fit$coefficients[3]^2)
    r <- suppressWarnings(stats::cor(y, as.numeric(X %*% fit$coefficients)))
    if (!is.finite(r)) next
    mri[k] <- 2 * A * r; amp[k] <- A; off[k] <- fit$coefficients[1]; rr[k] <- r
  }
  ok <- which(is.finite(mri))
  if (!length(ok)) return(NULL)
  ng <- length(grid)
  interior <- ok[which(ok > 1L & ok < ng &
                 mri[ok] > mri[pmax(ok - 1L, 1L)] & mri[ok] > mri[pmin(ok + 1L, ng)])]
  cand <- if (length(interior)) interior else ok
  k <- cand[which.max(mri[cand])]
  w <- 2 * pi / grid[k]
  X <- cbind(1, cos(w * x), sin(w * x))
  yhat <- as.numeric(X %*% stats::lm.fit(X, y)$coefficients)
  list(period_min = grid[k] * bin_min, amplitude = amp[k], offset = off[k],
       pearson_r = rr[k], MRI = mri[k], n_epochs = n,
       lids = y, fitted = yhat, bin_min = bin_min)
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
