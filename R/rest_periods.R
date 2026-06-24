# Centered rolling mean with a min-periods gate (O(n) via cumsums), reproducing
# pandas rolling(win, center = TRUE, min_periods = min_win).mean().
.rp_trend <- function(x, win, min_win) {
  n <- length(x)
  if (n == 0L) return(numeric(0))
  xx <- ifelse(is.na(x), 0, x)
  ok <- as.integer(!is.na(x))
  cs <- c(0, cumsum(xx)); cn <- c(0, cumsum(ok))
  h_lo <- floor(win / 2); h_hi <- win - 1L - h_lo
  idx <- seq_len(n)
  lo <- pmax(idx - h_lo, 1L); hi <- pmin(idx + h_hi, n)
  cnt <- cn[hi + 1L] - cn[lo]
  tr <- (cs[hi + 1L] - cs[lo]) / cnt
  tr[cnt < min_win] <- NA_real_
  tr
}

# Onset indices of below-threshold runs of at least min_seed epochs.
.rp_seeds <- function(sw, min_seed) {
  v <- ifelse(is.na(sw), -1L, sw)
  r <- rle(v)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  starts[r$values == 1L & r$lengths >= min_seed]
}

# Pearson correlation of a binary slice against each "k+1 ones then zeros"
# template, in closed form (no per-template loop).
.rp_corr <- function(s) {
  L <- length(s)
  if (L < 2L) return(numeric(0))
  mbar <- mean(s)
  ss <- sqrt(sum((s - mbar)^2))
  if (ss == 0) return(rep(NA_real_, L - 1L))
  cs <- cumsum(s)
  m <- seq_len(L - 1L)
  r <- (cs[m] - m * mbar) / (sqrt(m * (L - m)) * ss / sqrt(L))
  r[!is.finite(r)] <- NA_real_
  r
}

# Highest correlation peak that exceeds the next (n_succ + 1) values, the short
# above-threshold tolerance that keeps a bout from splitting on brief blips.
.rp_peak <- function(corr, n_succ) {
  L <- length(corr)
  if (L == 0L) return(NA_integer_)
  look <- n_succ + 1L
  cand <- integer(0)
  for (k in seq_len(L)) {
    hi <- min(k + look, L)
    if (k == L || isTRUE(all(corr[k] > corr[(k + 1L):hi], na.rm = TRUE))) cand <- c(cand, k)
  }
  cand <- cand[is.finite(corr[cand])]
  if (!length(cand)) return(NA_integer_)
  cand[which.max(corr[cand])]
}

#' Consolidated Rest-Period Detection (Roenneberg / MASDA)
#'
#' Detects all consolidated rest bouts across a recording with the Roenneberg
#' consolidation algorithm (the Munich Actimetry Sleep Detection Algorithm).
#' Each epoch is compared to a fraction of its own 24-hour activity trend, runs
#' below that threshold seed candidate rest bouts, and a correlation procedure
#' grows each seed into a consolidated bout. Unlike \code{\link{sleep.changepoints}}
#' (one nightly bout per cycle), this returns any number of bouts,
#' including daytime naps and fragmented or polyphasic rest.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param epoch_length Epoch length in seconds (default 60).
#' @param trend_period Moving-average window for the activity trend, in seconds
#'   (default 86400, 24 hours).
#' @param min_trend_period Minimum non-missing seconds required in the trend
#'   window (default 43200, 12 hours).
#' @param threshold Fraction of the local trend below which an epoch is candidate
#'   rest (default 0.15).
#' @param min_seed_period Minimum below-threshold run to seed a bout, in seconds
#'   (default 1800, 30 minutes).
#' @param max_test_period Maximum bout length the consolidation tests, in seconds
#'   (default 43200, 12 hours).
#' @param r_consec_below Above-threshold tolerance during consolidation, in
#'   seconds (default 1800, 30 minutes).
#' @param nap_max_minutes A non-main bout shorter than this is labelled a nap
#'   (default 180).
#'
#' @return An object of class \code{actiRhythm_roenneberg}: a \code{rest_periods}
#'   data frame (one row per consolidated bout, with onset, offset, duration, and
#'   a main/nap label), the per-epoch \code{rest} vector, the activity
#'   \code{trend}, and per-bout counts. The function never errors; with no
#'   resolvable bout it returns an empty \code{rest_periods}.
#'
#' @references
#' Roenneberg T, Keller LK, Fischer D, Matera JL, Vetter C, Winnebeck EC (2015).
#' Human activity and rest in situ. \emph{Methods in Enzymology}, 552:257-283.
#' \doi{10.1016/bs.mie.2014.11.028}
#'
#' Loock A-S, Khan Sullivan A, Reis C, et al. (2021). Validation of the Munich
#' Actimetry Sleep Detection Algorithm for estimating sleep-wake patterns from
#' activity recordings. \emph{Journal of Sleep Research}, 30(6):e13371.
#' \doi{10.1111/jsr.13371}
#'
#' Hammad G, Reyt M, Beliy N, et al. (2021). pyActigraphy: open-source python
#' package for actigraphy data visualization and analysis. \emph{PLOS
#' Computational Biology}, 17(10):e1009514. \doi{10.1371/journal.pcbi.1009514}
#'
#' @seealso \code{\link{sleep.changepoints}}, \code{\link{sleep.cole.kripke}}
#'
#' @examples
#' # Three nights of rest plus one daytime nap
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
#' day <- as.integer(format(ts, "%j"))
#' counts <- ifelse(h >= 23 | h < 7, 5, 300)
#' counts[h >= 14 & h < 15.5 & day == min(day) + 1L] <- 5   # a nap on day 2
#' rest.periods(counts, ts)
#'
#' @export
rest.periods <- function(counts, timestamps, epoch_length = 60,
                         trend_period = 86400, min_trend_period = 43200,
                         threshold = 0.15, min_seed_period = 1800,
                         max_test_period = 43200, r_consec_below = 1800,
                         nap_max_minutes = 180) {
  el <- as.numeric(epoch_length)
  trend_win <- round(trend_period / el)
  min_win   <- round(min_trend_period / el)
  min_seed  <- max(1L, round(min_seed_period / el))
  max_test  <- max(1L, round(max_test_period / el))
  n_succ    <- round(r_consec_below / el)

  empty_periods <- function() data.frame(
    bout = integer(0), onset = as.POSIXct(character(0)),
    offset = as.POSIXct(character(0)), onset_index = integer(0),
    offset_index = integer(0), duration_min = numeric(0),
    date = as.Date(character(0)), is_main = logical(0),
    type = factor(character(0), levels = c("main", "nap")),
    stringsAsFactors = FALSE)

  params <- list(epoch_length = el, trend_win = trend_win, min_win = min_win,
                 min_seed = min_seed, max_test = max_test, n_succ = n_succ,
                 threshold = threshold, nap_max_minutes = nap_max_minutes)
  out <- function(rp, rest, trend, thr) {
    nd <- if (nrow(rp)) length(unique(rp$date)) else 0L
    structure(list(
      rest_periods = rp, rest = rest, trend = trend, threshold_series = thr,
      n_bouts = nrow(rp),
      total_rest_min = sum(rest == 1, na.rm = TRUE) * el / 60,
      n_days = nd, mean_bouts_per_day = if (nd) nrow(rp) / nd else NA_real_,
      params = params, epoch_length = el),
      class = c("actiRhythm_roenneberg", "list"))
  }

  x <- suppressWarnings(as.numeric(counts))
  ts <- timestamps
  n <- length(x)
  if (n != length(ts)) stop("counts and timestamps must have same length")
  if (n < 2L) return(out(empty_periods(), rep(NA_integer_, n), rep(NA_real_, n), rep(NA_real_, n)))

  trend <- .rp_trend(x, trend_win, min_win)
  thr <- threshold * trend
  sw <- ifelse(!is.na(trend) & x <= thr, 1L, 0L)
  sw[is.na(trend)] <- NA_integer_

  seeds <- .rp_seeds(sw, min_seed)
  if (!length(seeds)) return(out(empty_periods(), sw, trend, thr))

  first <- seeds[1L]
  if (first > 1L) { z <- seq_len(first - 1L); sw[z][which(sw[z] == 1L)] <- 0L }

  last_off <- 0L; sot <- list()
  for (sd0 in seeds) {
    if (sd0 <= last_off) next
    if (last_off >= 1L && sd0 > last_off + 1L) {
      gap <- (last_off + 1L):(sd0 - 1L); sw[gap][which(sw[gap] == 1L)] <- 0L
    }
    hi <- min(sd0 + max_test - 1L, n)
    s <- sw[sd0:hi]; s[is.na(s)] <- 0L
    pk <- .rp_peak(.rp_corr(s), n_succ)
    if (is.na(pk)) next
    offset <- sd0 + (pk - 1L)
    sw[sd0:offset] <- 1L
    sot[[length(sot) + 1L]] <- c(sd0, offset)
    last_off <- offset
  }
  if (last_off >= 1L && last_off < n) {
    tl <- (last_off + 1L):n; sw[tl][which(sw[tl] == 1L)] <- 0L
  }

  if (!length(sot)) return(out(empty_periods(), sw, trend, thr))
  onsets <- vapply(sot, `[`, integer(1), 1L)
  offsets <- vapply(sot, `[`, integer(1), 2L)
  rp <- data.frame(
    bout = seq_along(onsets), onset = ts[onsets], offset = ts[offsets],
    onset_index = onsets, offset_index = offsets,
    duration_min = (offsets - onsets + 1L) * el / 60,
    date = as.Date(ts[onsets]), is_main = FALSE, stringsAsFactors = FALSE)
  for (dd in unique(rp$date)) {
    idx <- which(rp$date == dd)
    rp$is_main[idx[which.max(rp$duration_min[idx])]] <- TRUE
  }
  rp$type <- factor(ifelse(rp$duration_min < nap_max_minutes & !rp$is_main, "nap", "main"),
                    levels = c("main", "nap"))
  out(rp, sw, trend, thr)
}


#' @export
print.actiRhythm_roenneberg <- function(x, ...) {
  cat("Consolidated Rest Periods (Roenneberg/MASDA)\n\n")
  cat(sprintf("  Bouts detected:   %d  (%.1f per day over %d days)\n",
              x$n_bouts, x$mean_bouts_per_day, x$n_days))
  cat(sprintf("  Total rest:       %.1f h\n", x$total_rest_min / 60))
  if (x$n_bouts) {
    naps <- sum(x$rest_periods$type == "nap")
    cat(sprintf("  Main bouts / naps: %d / %d\n", x$n_bouts - naps, naps))
    cat("\n  First bouts:\n")
    e <- utils::head(x$rest_periods, 3)
    for (i in seq_len(nrow(e)))
      cat(sprintf("    %-4s %s -> %s  (%.0f min)\n", as.character(e$type[i]),
                  format(e$onset[i], "%m-%d %H:%M"),
                  format(e$offset[i], "%m-%d %H:%M"), e$duration_min[i]))
  }
  cat("\n  Reference: Roenneberg et al. (2015); Loock et al. (2021)\n\n")
  invisible(x)
}
