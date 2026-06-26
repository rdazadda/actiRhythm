# Raw-native, diary-free sleep detection from the z-angle (van Hees 2015/2018).
# Posture-based: counts cannot reproduce this because the band-pass filter removes
# the gravity component the arm angle is built from.

# Consolidate a logical no-movement mask into the single longest sustained block:
# drop TRUE blocks <= min_ep, fill FALSE gaps < gap_ep, keep the longest TRUE run.
.spt_consolidate <- function(nomov, min_ep, gap_ep) {
  flip <- function(v, keep_val, lim, op) {
    r <- rle(v); e <- cumsum(r$lengths); s <- e - r$lengths + 1L
    for (i in which(r$values == keep_val & op(r$lengths, lim))) v[s[i]:e[i]] <- !keep_val
    v
  }
  nomov <- flip(nomov, TRUE,  min_ep, `<=`)   # remove short no-move blocks
  nomov <- flip(nomov, FALSE, gap_ep, `<`)    # fill short movement gaps
  r <- rle(nomov); e <- cumsum(r$lengths); s <- e - r$lengths + 1L
  tr <- which(r$values)
  if (!length(tr)) return(c(NA_integer_, NA_integer_))
  best <- tr[which.max(r$lengths[tr])]
  c(s[best], e[best])
}

#' Non-Wear Detection from Raw Acceleration
#'
#' Flags non-wear time from raw acceleration by the van Hees et al. (2011)
#' standard-deviation-and-range rule: a block is non-wear when, over a window
#' centred on it (default 60 minutes), at least two of the three axes have both a
#' standard deviation below \code{sd_crit} and a value range below
#' \code{range_crit}. A stationary, taken-off device reads as non-wear. This lets
#' the z-angle sleep detector tell device-off periods from real sleep (a still arm
#' that keeps micro-movement). Pass the result as the
#' \code{wear} argument of \code{\link{rest.spt}}.
#'
#' @param x A path to a raw file or a raw data frame (see \code{\link{raw.metrics}}).
#' @param device Device brand or \code{"auto"} (file input only).
#' @param epoch Output epoch length in seconds for the returned mask (default 5,
#'   matching \code{\link{rest.spt}}).
#' @param block Internal classification block length in seconds (default 300).
#' @param window Window in seconds over which the SD and range are taken (default 3600).
#' @param sd_crit Per-axis SD threshold in g (default 0.013).
#' @param range_crit Per-axis range threshold in g (default 0.050).
#' @param tz Time zone (default \code{"UTC"}).
#'
#' @return A logical vector, one per epoch: \code{TRUE} = wear, \code{FALSE} =
#'   non-wear. Never errors.
#'
#' @references
#' \insertRef{vanhees2011}{actiRhythm}
#'
#' @seealso \code{\link{rest.spt}}, \code{\link{raw.metrics}}
#'
#' @examples
#' \donttest{
#' raw  <- example_raw(days = 2, device_off = 1)   # two worn days + one device-off day
#' mean(detect.nonwear.raw(raw, epoch = 60))       # fraction of epochs worn
#' }
#'
#' @export
detect.nonwear.raw <- function(x, device = "auto", epoch = 5, block = 300,
                               window = 3600, sd_crit = 0.013, range_crit = 0.050,
                               tz = "UTC") {
  r <- if (is.character(x)) .raw_xyz(x, device, tz) else .as_raw_xyz(x)
  X <- cbind(r$x, r$y, r$z); fs <- r$fs; n <- nrow(X)
  n_out <- n %/% round(fs * epoch)
  if (n_out < 1L) return(logical(0))
  spb <- round(fs * block); nb <- n %/% spb
  if (nb < 2L) return(rep(TRUE, n_out))                  # too short -> assume worn
  g  <- rep(seq_len(nb), each = spb)[seq_len(nb * spb)]
  Xt <- X[seq_len(nb * spb), , drop = FALSE]
  S  <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, sum)), numeric(nb))
  SS <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, function(v) sum(v^2))), numeric(nb))
  MN <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, min)), numeric(nb))
  MX <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, max)), numeric(nb))
  cS <- apply(rbind(0, S), 2, cumsum); cSS <- apply(rbind(0, SS), 2, cumsum)
  wb <- max(1L, round(window / block)); half <- wb %/% 2L
  score <- integer(nb)
  for (i in seq_len(nb)) {
    # exact `window`-long span centred on block i, with GGIR's fixed first/last
    # windows at the recording edges (g.getmeta detect_nonwear_clipping, 2013).
    if (i <= half)            { lo <- 1L;                    hi <- min(nb, wb) }
    else if (i > nb - half)   { lo <- max(1L, nb - wb + 1L); hi <- nb }
    else                      { lo <- i - half;              hi <- min(nb, i - half + wb - 1L) }
    N <- (hi - lo + 1L) * spb
    cnt <- 0L
    for (a in 1:3) {
      s <- cS[hi + 1L, a] - cS[lo, a]; ss <- cSS[hi + 1L, a] - cSS[lo, a]
      sdv <- sqrt(max(ss / N - (s / N)^2, 0))
      rng <- max(MX[lo:hi, a]) - min(MN[lo:hi, a])
      if (sdv < sd_crit && rng < range_crit) cnt <- cnt + 1L
    }
    score[i] <- cnt
  }
  # GGIR bridging: a single-axis block between two multi-axis blocks is non-wear.
  for (i in which(score == 1L))
    if (i > 1L && i < nb && score[i - 1L] > 1L && score[i + 1L] > 1L) score[i] <- 2L
  nonwear <- score >= 2L
  blk <- pmin(floor((seq_len(n_out) - 1L) * epoch / block) + 1L, nb)
  !nonwear[blk]
}

#' Sleep-Period-Time Window from the z-Angle (HDCZA)
#'
#' Detects the main sleep-period-time (SPT) window per day directly from the
#' z-angle, with no sleep diary, using the van Hees et al. (2018) Heuristic
#' algorithm based on the Distribution of Change in Z-Angle (HDCZA): the absolute
#' 5-second change in arm angle is smoothed over 5 minutes, thresholded at a
#' fraction of its own distribution, and the longest sustained low-change block
#' per noon-to-noon day becomes the SPT. This is raw-native: there is no
#' count-based equivalent.
#'
#' @param anglez Numeric z-angle (degrees) per epoch, e.g. from
#'   \code{raw.metrics(..., metrics = "anglez")} at a short epoch.
#' @param timestamps POSIXct timestamps, one per epoch.
#' @param epoch_length Epoch length in seconds (default 5, the validated value).
#' @param pct Percentile of the change distribution for the threshold (default 10).
#' @param mult Multiplier on that percentile (default 15).
#' @param clamp Lower/upper clamp (degrees) on the threshold (default
#'   \code{c(0.13, 0.50)}).
#' @param min_block Minimum SPT block length in minutes (default 30).
#' @param max_gap Maximum movement gap to bridge within the SPT, minutes
#'   (default 60).
#' @param algo \code{"HDCZA"} (wrist, change-in-angle) or \code{"HorAngle"} (a
#'   hip variant thresholding the absolute angle at 60 degrees).
#' @param wear Optional logical wear mask, one per epoch (e.g. from
#'   \code{\link{detect.nonwear.raw}}); non-wear epochs are excluded from the SPT
#'   so a stationary, taken-off device is not scored as sleep.
#'
#' @return An object of class \code{actiRhythm_spt}: a data frame with one row per
#'   day (\code{date}, \code{onset}, \code{offset}, \code{duration} in hours, and
#'   the \code{threshold} used). Never errors; returns no rows if no day has a
#'   detectable window.
#'
#' @references
#' \insertRef{vanhees2018}{actiRhythm}
#'
#' @seealso \code{\link{sib.vanhees}}, \code{\link{sleep.from.spt}},
#'   \code{\link{raw.metrics}}
#'
#' @examples
#' # One still night inside an active day yields one SPT window
#' ts <- seq(as.POSIXct("2024-01-01 12:00", tz = "UTC"), by = 5, length.out = 17280)
#' h <- as.numeric(format(ts, "%H")); night <- h >= 23 | h < 7
#' set.seed(1)
#' anglez <- ifelse(night, -60, -30) + rnorm(17280, 0, ifelse(night, 0.02, 20))
#' rest.spt(anglez, ts, epoch_length = 5)
#'
#' @export
rest.spt <- function(anglez, timestamps, epoch_length = 5, pct = 10, mult = 15,
                     clamp = c(0.13, 0.50), min_block = 30, max_gap = 60,
                     algo = c("HDCZA", "HorAngle"), wear = NULL) {
  algo <- match.arg(algo)
  n <- length(anglez)
  empty <- structure(data.frame(date = as.Date(character(0)),
    onset = as.POSIXct(character(0)), offset = as.POSIXct(character(0)),
    duration = numeric(0), threshold = numeric(0)), class = c("actiRhythm_spt", "data.frame"))
  if (n < 2L || length(timestamps) != n) return(empty)
  if (!is.null(wear) && length(wear) != n) wear <- NULL   # ignore a mismatched mask
  per_min <- 60 / epoch_length
  k1 <- max(1L, round(5 * per_min)); if (k1 %% 2L == 0L) k1 <- k1 + 1L
  if (algo == "HDCZA") {
    dz <- abs(c(0, diff(anglez)))
    x  <- if (n >= k1) stats::runmed(dz, k1, endrule = "median") else dz
  } else {
    x <- abs(anglez)
  }
  day <- as.Date(timestamps - 12 * 3600)        # noon-to-noon day
  min_ep <- round(min_block * per_min); gap_ep <- round(max_gap * per_min)
  rows <- list()
  for (d in unique(day)) {
    idx <- which(day == d); if (length(idx) < min_ep) next
    xd <- x[idx]
    if (algo == "HDCZA") {
      thr <- min(max(stats::quantile(xd, pct / 100, names = FALSE) * mult,
                     clamp[1]), clamp[2])
    } else thr <- 60
    nomov <- xd < thr
    if (!is.null(wear)) nomov <- nomov & wear[idx]    # non-wear cannot be sleep
    blk <- .spt_consolidate(nomov, min_ep, gap_ep)
    if (anyNA(blk)) next
    on <- timestamps[idx[blk[1]]]; off <- timestamps[idx[blk[2]]]
    rows[[length(rows) + 1L]] <- data.frame(date = as.Date(d, origin = "1970-01-01"),
      onset = on, offset = off,
      duration = as.numeric(difftime(off, on, units = "hours")), threshold = thr)
  }
  if (!length(rows)) return(empty)
  structure(do.call(rbind, rows), class = c("actiRhythm_spt", "data.frame"))
}

#' Sustained-Inactivity-Bout Sleep Scoring from the z-Angle
#'
#' Scores each epoch sleep (\code{"S"}) or wake (\code{"W"}) from the z-angle by
#' the van Hees et al. (2015) sustained-inactivity-bout rule: an interval with no
#' arm-posture change exceeding \code{angle_thresh} degrees for at least
#' \code{time_thresh} minutes is sustained inactivity (sleep). The output has the
#' same shape as \code{\link{sleep.cole.kripke}}, so it feeds
#' \code{\link{sleep.regularity.index}}, \code{\link{lids}} and the rest consumers
#' directly. Intersect it with \code{\link{rest.spt}} via \code{\link{sleep.from.spt}}.
#'
#' @param anglez Numeric z-angle (degrees) per epoch.
#' @param angle_thresh Posture-change threshold in degrees (default 5).
#' @param time_thresh Minimum sustained-inactivity duration in minutes (default 5).
#' @param epoch_length Epoch length in seconds (default 5).
#'
#' @return Character vector of \code{"S"}/\code{"W"}, one per epoch. A day with
#'   fewer than 10 posture changes is scored all \code{"S"} (sustained
#'   inactivity, to be gated by the SPT window and wear time).
#'
#' @references
#' \insertRef{vanhees2015}{actiRhythm}
#'
#' @seealso \code{\link{rest.spt}}, \code{\link{sleep.from.spt}}
#'
#' @examples
#' # A long still stretch scores as sustained inactivity (sleep)
#' set.seed(1)
#' anglez <- c(rnorm(3000, -60, 0.02), rnorm(3000, 0, 20))   # still, then active
#' table(sib.vanhees(anglez, epoch_length = 5))
#'
#' @export
sib.vanhees <- function(anglez, angle_thresh = 5, time_thresh = 5, epoch_length = 5) {
  n <- length(anglez)
  if (n < 2L) return(rep("W", n))
  postch <- which(abs(c(0, diff(anglez))) > angle_thresh)
  if (length(postch) < 10L) return(rep("S", n))
  gap_ep <- time_thresh * (60 / epoch_length)
  sib <- rep("W", n)
  q <- which(diff(postch) > gap_ep)
  for (i in q) sib[postch[i]:postch[i + 1L]] <- "S"
  sib
}

#' Sleep Parameters from SPT and Sustained-Inactivity Bouts
#'
#' Intersects the per-epoch sustained-inactivity sleep score
#' (\code{\link{sib.vanhees}}) with the sleep-period-time window
#' (\code{\link{rest.spt}}) to derive per-night sleep parameters: total sleep
#' time, onset and wake, wake-after-sleep-onset, efficiency, and awakenings.
#'
#' @param spt An \code{actiRhythm_spt} object from \code{\link{rest.spt}}.
#' @param sib A character \code{"S"}/\code{"W"} vector from \code{\link{sib.vanhees}}.
#' @param timestamps POSIXct timestamps, one per epoch (matching \code{sib}).
#' @param epoch_length Epoch length in seconds (default 5).
#'
#' @return An object of class \code{actiRhythm_sleep}: a data frame with one row
#'   per night (\code{date}, sleep \code{onset}, \code{offset}, \code{tst} hours,
#'   \code{waso} minutes, \code{efficiency}, \code{n_awakenings}, \code{mid_sleep}).
#'
#' @references
#' \insertRef{vanhees2018}{actiRhythm}
#'
#' @seealso \code{\link{rest.spt}}, \code{\link{sib.vanhees}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01 12:00", tz = "UTC"), by = 5, length.out = 17280)
#' h <- as.numeric(format(ts, "%H")); night <- h >= 23 | h < 7
#' set.seed(1)
#' anglez <- ifelse(night, -60, -30) + rnorm(17280, 0, ifelse(night, 0.02, 20))
#' spt <- rest.spt(anglez, ts, epoch_length = 5)
#' sib <- sib.vanhees(anglez, epoch_length = 5)
#' sleep.from.spt(spt, sib, ts, epoch_length = 5)
#'
#' @export
sleep.from.spt <- function(spt, sib, timestamps, epoch_length = 5) {
  empty <- structure(data.frame(date = as.Date(character(0)),
    onset = as.POSIXct(character(0)), offset = as.POSIXct(character(0)),
    tst = numeric(0), waso = numeric(0), efficiency = numeric(0),
    n_awakenings = integer(0), mid_sleep = as.POSIXct(character(0))),
    class = c("actiRhythm_sleep", "data.frame"))
  if (!nrow(spt) || !length(sib)) return(empty)
  rows <- list()
  for (i in seq_len(nrow(spt))) {
    inw <- which(timestamps >= spt$onset[i] & timestamps <= spt$offset[i])
    if (!length(inw)) next
    s <- sib[inw]; sleep_ep <- which(s == "S")
    if (!length(sleep_ep)) next
    span    <- sleep_ep[1]:sleep_ep[length(sleep_ep)]
    on_s    <- timestamps[inw[sleep_ep[1]]]
    wake_s  <- timestamps[inw[sleep_ep[length(sleep_ep)]]]
    waso_ep <- sum(s[span] == "W")
    rows[[length(rows) + 1L]] <- data.frame(date = spt$date[i],
      onset = on_s, offset = wake_s,
      tst = sum(s == "S") * epoch_length / 3600,
      waso = waso_ep * epoch_length / 60,
      efficiency = sum(s == "S") / length(inw),
      n_awakenings = sum(rle(s[span] == "W")$values),
      mid_sleep = on_s + as.numeric(difftime(wake_s, on_s, units = "secs")) / 2)
  }
  if (!length(rows)) return(empty)
  structure(do.call(rbind, rows), class = c("actiRhythm_sleep", "data.frame"))
}
