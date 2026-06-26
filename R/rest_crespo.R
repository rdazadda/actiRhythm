# Flag runs of consecutive zeros longer than zeta as invalid (device-off /
# non-wear), returned as a logical mask the length of x.
.cr_mask <- function(x, zeta) {
  r <- rle(x == 0)
  bad <- r$values & r$lengths > zeta
  inverse.rle(list(lengths = r$lengths, values = bad))
}

# Centered running median of odd window k, edges padded with padval so the ends
# are biased toward the pad (active, when padded with the series max).
.cr_med <- function(x, k, padval) {
  if (k %% 2L == 0L) k <- k + 1L
  if (length(x) < k) k <- if (length(x) %% 2L == 0L) length(x) - 1L else length(x)
  if (k < 3L) return(x)
  h <- (k - 1L) %/% 2L
  xp <- c(rep(padval, h), x, rep(padval, h))
  stats::runmed(xp, k = k, endrule = "median")[(h + 1L):(h + length(x))]
}

# Binary closing (fill 0-runs shorter than L) then opening (drop 1-runs shorter
# than L): morphology with a flat structuring element via run-length flipping.
.cr_morph <- function(y, L) {
  flip <- function(y, target) {
    r <- rle(y)
    r$values[r$values == target & r$lengths < L] <- 1L - target
    inverse.rle(r)
  }
  flip(flip(y, 0L), 1L)
}

#' Crespo Rest/Activity Period Detection
#'
#' Detects all consolidated rest bouts across a recording with the Crespo
#' algorithm (Crespo et al. 2012): a two-pass rank-order (median) filter and
#' mathematical-morphology pipeline that turns activity counts into a binary
#' rest/activity series and reads every rest bout from its transitions. Like
#' \code{\link{rest.periods}} it returns any number of bouts, naps and
#' fragmented rest included; it reaches that result by morphological filtering
#' rather than by Roenneberg-style consolidation, so the two differ in method.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param epoch_length Epoch length in seconds (default 60).
#' @param zeta,zeta_r,zeta_a Maximum valid consecutive-zero run, in epochs, for
#'   the global pre-conditioning, within rest segments, and within active
#'   segments (defaults 15, 30, 2). Longer zero runs are treated as non-wear.
#' @param t Quantile of the counts used to replace non-wear epochs (default
#'   0.33).
#' @param alpha Expected daily rest length, in seconds (default 28800, 8 hours);
#'   sets the threshold percentile (alpha / 24 h) and the minimum data required.
#' @param beta Filter and morphology scale, in seconds (default 3600, 1 hour);
#'   sets the rank-order filter window and the structuring-element size.
#'
#' @return An object of class \code{actiRhythm_crespo}: a \code{rest_periods}
#'   data frame (one row per bout, with onset, offset, and duration), the
#'   per-epoch \code{rest_state} ("R"/"A"), and per-bout counts. The function
#'   never errors; with no resolvable bout it returns an empty
#'   \code{rest_periods}.
#'
#' @references
#' \insertRef{crespo2012}{actiRhythm}
#'
#' \insertRef{hammad2021}{actiRhythm}
#'
#' @seealso \code{\link{rest.periods}}, \code{\link{sleep.changepoints}}
#'
#' @examples
#' # Two nights of rest with a daytime nap between them
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
#' h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
#' day <- as.integer(format(ts, "%j"))
#' counts <- ifelse(h >= 23 | h < 7, 5, 300)
#' counts[h >= 14 & h < 15 & day == min(day) + 1L] <- 5
#' rest.crespo(counts, ts)
#'
#' @export
rest.crespo <- function(counts, timestamps, epoch_length = 60, zeta = 15,
                        zeta_r = 30, zeta_a = 2, t = 0.33,
                        alpha = 8 * 3600, beta = 1 * 3600) {
  el <- as.numeric(epoch_length)
  L_w <- round(alpha / el) + 1L
  L_p <- round(beta / el) + 1L
  pct <- alpha / 86400

  empty_periods <- function() data.frame(
    bout = integer(0), onset = as.POSIXct(character(0)),
    offset = as.POSIXct(character(0)), onset_index = integer(0),
    offset_index = integer(0), duration_min = numeric(0),
    date = as.Date(character(0)), stringsAsFactors = FALSE)

  out <- function(rp, rest_state) {
    nd <- if (nrow(rp)) length(unique(rp$date)) else 0L
    structure(list(
      rest_periods = rp, rest_state = rest_state, n_rest_periods = nrow(rp),
      total_rest_min = sum(rest_state == "R", na.rm = TRUE) * el / 60,
      n_days = nd, mean_bouts_per_day = if (nd) nrow(rp) / nd else NA_real_,
      params = list(epoch_length = el, L_w = L_w, L_p = L_p, threshold_pct = pct,
                    zeta = zeta, zeta_r = zeta_r, zeta_a = zeta_a, t = t)),
      class = c("actiRhythm_crespo", "list"))
  }

  x <- suppressWarnings(as.numeric(counts))
  ts <- timestamps
  n <- length(x)
  if (n != length(ts)) stop("counts and timestamps must have same length")
  if (n < L_w || all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0)
    return(out(empty_periods(), rep(NA_character_, n)))
  x[is.na(x)] <- 0

  s_t <- as.numeric(stats::quantile(x, t, na.rm = TRUE))
  s_max <- max(x)

  # Pass 1: condition, median-filter, threshold, consolidate. The rank-order
  # filter runs at the beta (morphology) scale; a window near the rest length
  # would erase the rest itself. alpha sets the threshold percentile below.
  x1 <- x; x1[.cr_mask(x, zeta)] <- s_t
  xf1 <- .cr_med(x1, L_p, s_max)
  thr <- as.numeric(stats::quantile(xf1, pct, na.rm = TRUE))
  y1 <- as.integer(xf1 > thr)               # 1 = active, 0 = rest
  ye <- .cr_morph(y1, L_p)

  # Pass 2: re-mask rest and active segments separately, refilter, consolidate.
  mask2 <- logical(n)
  rest_seg <- ye == 0L
  if (any(rest_seg))  mask2[rest_seg]  <- .cr_mask(ifelse(rest_seg, x, 1), zeta_r)[rest_seg]
  if (any(!rest_seg)) mask2[!rest_seg] <- .cr_mask(ifelse(!rest_seg, x, 1), zeta_a)[!rest_seg]
  x2 <- x; x2[mask2] <- s_t                  # hold invalid epochs at the resting level
  xf2 <- .cr_med(x2, L_p, s_max)
  y2 <- as.integer(xf2 > thr)
  ba <- .cr_morph(y2, 2L * (L_p - 1L) + 1L)
  ba[1] <- 1L; ba[n] <- 1L                   # force active endpoints

  rest_state <- ifelse(ba == 0L, "R", "A")
  d <- diff(ba)
  starts <- which(d == -1L) + 1L             # rest onset epochs
  ends <- which(d == 1L)                     # rest offset epochs (last rest epoch)
  if (!length(starts) || !length(ends)) return(out(empty_periods(), rest_state))
  m <- min(length(starts), length(ends))
  starts <- starts[seq_len(m)]; ends <- ends[seq_len(m)]

  rp <- data.frame(
    bout = seq_len(m), onset = ts[starts], offset = ts[ends],
    onset_index = starts, offset_index = ends,
    duration_min = (ends - starts + 1L) * el / 60,
    date = as.Date(ts[starts]), stringsAsFactors = FALSE)
  out(rp, rest_state)
}


#' @export
print.actiRhythm_crespo <- function(x, ...) {
  cat("Crespo Rest/Activity Periods\n\n")
  cat(sprintf("  Rest bouts:   %d  (%.1f per day over %d days)\n",
              x$n_rest_periods, x$mean_bouts_per_day, x$n_days))
  cat(sprintf("  Total rest:   %.1f h\n", x$total_rest_min / 60))
  if (x$n_rest_periods) {
    cat("\n  First bouts:\n")
    e <- utils::head(x$rest_periods, 3)
    for (i in seq_len(nrow(e)))
      cat(sprintf("    %s -> %s  (%.0f min)\n",
                  format(e$onset[i], "%m-%d %H:%M"),
                  format(e$offset[i], "%m-%d %H:%M"), e$duration_min[i]))
  }
  cat("\n  Reference: Crespo et al. (2012)\n\n")
  invisible(x)
}
