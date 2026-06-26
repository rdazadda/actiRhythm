# Gap-aware intradaily variability of a series binned to bin_min resolution,
# matching the consecutive-pair convention of .calculate.IS.IV.
.iv_at_bin <- function(counts, timestamps, bin_min) {
  t0 <- min(timestamps, na.rm = TRUE)
  bin <- floor(as.numeric(difftime(timestamps, t0, units = "mins")) / bin_min)
  agg <- tapply(counts, bin, mean, na.rm = TRUE)
  b <- as.integer(names(agg)); x <- as.numeric(agg)
  o <- order(b); b <- b[o]; x <- x[o]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  xbar <- mean(x, na.rm = TRUE)
  consec <- which(diff(b) == 1L)
  den <- sum((x - xbar)^2, na.rm = TRUE) / n
  if (!length(consec) || den == 0) return(NA_real_)
  num <- sum((x[consec + 1L] - x[consec])^2, na.rm = TRUE) / length(consec)
  num / den
}

#' Multiscale Intradaily Variability (IVm)
#'
#' Intradaily variability computed across a set of bin sizes and averaged, the
#' counterpart to multiscale interdaily stability (\code{\link{circadian.is.multiscale}}).
#' Fragmentation that is invisible at the hourly scale often shows at finer bins,
#' so the averaged IVm varies less across recordings than the single
#' hourly IV (Goncalves et al. 2014).
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param bin_minutes Bin sizes in minutes (default 5 to 60).
#'
#' @return An object of class \code{actiRhythm_ivm}: a per-bin \code{IV} table and
#'   the averaged \code{IVm}. Never errors; returns \code{NA} on insufficient data.
#'
#' @references
#' \insertRef{goncalves2014}{actiRhythm}
#'
#' @seealso \code{\link{circadian.is.multiscale}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' intradaily.variability.multiscale(ifelse(h >= 23 | h < 7, 5, 300), ts)
#'
#' @export
intradaily.variability.multiscale <- function(counts, timestamps,
                                              bin_minutes = c(5, 10, 15, 30, 60)) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  keep <- is.finite(suppressWarnings(as.numeric(counts))) & !is.na(timestamps)
  counts <- as.numeric(counts)[keep]; timestamps <- timestamps[keep]
  iv <- vapply(bin_minutes, function(b) .iv_at_bin(counts, timestamps, b), numeric(1))
  structure(list(
    table = data.frame(bin_minutes = bin_minutes, IV = iv),
    IVm = if (all(is.na(iv))) NA_real_ else mean(iv, na.rm = TRUE)),
    class = c("actiRhythm_ivm", "list"))
}

#' @export
print.actiRhythm_ivm <- function(x, ...) {
  cat("Multiscale Intradaily Variability\n\n")
  cat(sprintf("  IVm (averaged): %.3f\n\n", x$IVm))
  print(x$table, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' Generalized Least- and Most-Active Periods (MX / LX)
#'
#' The least-active \code{LX} and most-active \code{MX} periods for arbitrary
#' window lengths, generalizing L5 and M10 (Van Someren et al. 1999). For each
#' window the function returns the mean activity over the window and its onset
#' and midpoint clock times. A fixed L5/M10 pair does not give these phase markers.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param windows Window lengths in hours (default \code{c(5, 10)} = L5, M10).
#'
#' @return An object of class \code{actiRhythm_mxlx}: a data frame with, for each
#'   window, the least- and most-active mean level, onset hour, and midpoint hour.
#'   Never errors.
#'
#' @references
#' \insertRef{vansomeren1999}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' activity.extrema(ifelse(h >= 23 | h < 7, 5, 300), ts, windows = c(5, 8, 10))
#'
#' @export
activity.extrema <- function(counts, timestamps, windows = c(5, 10)) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  keep <- is.finite(x) & !is.na(timestamps)
  x <- x[keep]; ts <- timestamps[keep]

  # Average 24h profile on a 1-minute grid (minute-of-day mean across days).
  mod <- as.POSIXlt(ts)$hour * 60 + as.POSIXlt(ts)$min
  prof <- tapply(x, factor(mod, levels = 0:1439), mean, na.rm = TRUE)
  prof <- as.numeric(prof)
  na_out <- function() structure(list(table = data.frame(), profile = prof,
    insufficient = TRUE), class = c("actiRhythm_mxlx", "list"))
  if (all(is.na(prof))) return(na_out())
  prof[is.na(prof)] <- mean(prof, na.rm = TRUE)
  cprof <- c(prof, prof)                      # wrap midnight

  rows <- lapply(windows, function(w) {
    k <- as.integer(round(w * 60))
    if (k < 1L || k > 1440L) return(NULL)
    csum <- cumsum(c(0, cprof))
    wsum <- csum[(k + 1L):(k + 1440L)] - csum[1:1440]   # window sums starting each minute
    lo <- which.min(wsum); hi <- which.max(wsum)
    data.frame(window_h = w,
      L_mean = wsum[lo] / k, L_onset_h = (lo - 1) / 60, L_mid_h = ((lo - 1 + k / 2) %% 1440) / 60,
      M_mean = wsum[hi] / k, M_onset_h = (hi - 1) / 60, M_mid_h = ((hi - 1 + k / 2) %% 1440) / 60)
  })
  tab <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(tab)) return(na_out())
  structure(list(table = tab, profile = prof, insufficient = FALSE),
            class = c("actiRhythm_mxlx", "list"))
}

#' @export
print.actiRhythm_mxlx <- function(x, ...) {
  cat("Least / Most Active Periods (MX / LX)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  t <- x$table
  for (i in seq_len(nrow(t)))
    cat(sprintf("  L%g: %.0f at %04.1fh   M%g: %.0f at %04.1fh\n",
                t$window_h[i], t$L_mean[i], t$L_onset_h[i],
                t$window_h[i], t$M_mean[i], t$M_onset_h[i]))
  cat("\n")
  invisible(x)
}


#' Dichotomy Index (I < O)
#'
#' The fraction of in-rest activity counts that fall below the median of the
#' out-of-rest (active) counts, a rest/active separation index used in cancer
#' chronobiology (Mormont et al. 2000). A high I<O means rest is quiet relative
#' to the active day, marking a well-separated rhythm.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param rest A logical vector (TRUE = rest / in-bed), or a character vector of
#'   states where \code{"R"}/\code{"S"}/\code{"sleep"} mark rest. Same length as
#'   \code{counts}.
#'
#' @return An object of class \code{actiRhythm_dichotomy}: the index \code{IO}
#'   (percent), the active-period median, and epoch counts. Never errors.
#'
#' @references
#' \insertRef{mormont2000}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' counts <- ifelse(h >= 23 | h < 7, 5, 300)
#' dichotomy.index(counts, ts, rest = h >= 23 | h < 7)
#'
#' @export
dichotomy.index <- function(counts, timestamps, rest) {
  if (length(counts) != length(rest))
    stop("counts and rest must have same length")
  x <- suppressWarnings(as.numeric(counts))
  r <- if (is.logical(rest)) rest
       else tolower(as.character(rest)) %in% c("r", "s", "sleep", "rest", "true", "1")
  ok <- is.finite(x) & !is.na(r)
  x <- x[ok]; r <- r[ok]
  na_out <- structure(list(IO = NA_real_, active_median = NA_real_,
    n_rest = sum(r), n_active = sum(!r)), class = c("actiRhythm_dichotomy", "list"))
  if (!any(r) || !any(!r)) return(na_out)
  o_med <- stats::median(x[!r], na.rm = TRUE)
  structure(list(
    IO = 100 * mean(x[r] < o_med, na.rm = TRUE),
    active_median = o_med, n_rest = sum(r), n_active = sum(!r)),
    class = c("actiRhythm_dichotomy", "list"))
}

#' @export
print.actiRhythm_dichotomy <- function(x, ...) {
  cat("Dichotomy Index (I<O)\n\n")
  cat(sprintf("  I<O:             %.1f%%\n", x$IO))
  cat(sprintf("  Active median:   %.1f counts\n", x$active_median))
  cat(sprintf("  Rest / active epochs: %d / %d\n\n", x$n_rest, x$n_active))
  invisible(x)
}


#' Per-Day Nonparametric Metrics
#'
#' The intraday nonparametric metrics computed for each day separately, showing
#' within-recording drift that the pooled values hide. Interdaily stability is
#' a between-day measure and so is not a per-day quantity; the per-day
#' table reports L5, M10, their onset times, relative amplitude, and intradaily
#' variability (Goncalves et al. 2014).
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param min_hours Minimum hours of data a day needs to be reported (default 12).
#'
#' @return An object of class \code{actiRhythm_daily}: a \code{daily} data frame,
#'   one row per day. Never errors.
#'
#' @references
#' \insertRef{goncalves2014}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' circadian.daily(ifelse(h >= 23 | h < 7, 5, 300), ts)
#'
#' @export
circadian.daily <- function(counts, timestamps, min_hours = 12) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  keep <- !is.na(timestamps)
  x <- x[keep]; ts <- timestamps[keep]
  day <- as.Date(ts)
  rows <- lapply(sort(unique(day)), function(d) {
    d <- as.Date(d, origin = "1970-01-01")
    idx <- day == d; xd <- x[idx]; td <- ts[idx]
    if (sum(is.finite(xd)) < min_hours * 60) return(NULL)
    mod <- as.POSIXlt(td)$hour * 60 + as.POSIXlt(td)$min
    prof <- as.numeric(tapply(xd, factor(mod, levels = 0:1439), mean, na.rm = TRUE))
    if (all(is.na(prof))) return(NULL)
    prof[is.na(prof)] <- mean(prof, na.rm = TRUE)
    csum <- cumsum(c(0, c(prof, prof)))
    l5 <- csum[301:1740] - csum[1:1440]; m10 <- csum[601:2040] - csum[1:1440]
    L5 <- min(l5) / 300; M10 <- max(m10) / 600
    data.frame(date = d, L5 = L5, L5_onset_h = (which.min(l5) - 1) / 60,
      M10 = M10, M10_onset_h = (which.max(m10) - 1) / 60,
      RA = if (M10 + L5 > 0) (M10 - L5) / (M10 + L5) else NA_real_,
      IV = .iv_at_bin(xd, td, 60), total = sum(xd, na.rm = TRUE))
  })
  tab <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  structure(list(daily = if (is.null(tab)) data.frame() else tab,
    n_days = if (is.null(tab)) 0L else nrow(tab)),
    class = c("actiRhythm_daily", "list"))
}

#' @export
print.actiRhythm_daily <- function(x, ...) {
  cat(sprintf("Per-Day Nonparametric Metrics (%d days)\n\n", x$n_days))
  if (x$n_days) {
    d <- x$daily
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], round, 3)
    print(d, row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}
