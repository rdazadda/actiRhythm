# Single change point in a segment by the least-squares (mean-shift) cost:
# the split k minimising the within-segment residual sum of squares. Vectorised.
.circa_cp_detect <- function(seg, min_seg = 5L) {
  n <- length(seg)
  if (n < 2L * min_seg + 1L) return(NA_integer_)
  cs <- cumsum(seg); cs2 <- cumsum(seg^2)
  tot <- cs[n]; tot2 <- cs2[n]
  k  <- seq.int(min_seg, n - min_seg)
  sl <- cs[k];        nl <- k
  sr <- tot - sl;     nr <- n - k
  cost <- (cs2[k] - sl^2 / nl) + ((tot2 - cs2[k]) - sr^2 / nr)
  k[which.min(cost)]
}

#' Change-Point Detection of Sleep and Wake Onsets
#'
#' Locates the sleep-onset and wake-onset time of each circadian cycle with a
#' cosinor-anchored mean-shift change point. A fixed 24-hour cosinor bounds each
#' rest and active span roughly (the cosinor anchoring follows CircaCP, Chen and
#' Sun 2024), and the precise transition inside each bound is then placed with a
#' single least-squares mean-shift change point on the raw counts. The result is a
#' per-night sleep-onset / wake-onset table.
#' A single rest-activity transition rate (such as
#' \code{\link{state.transitions}}) cannot localise this timing.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Cosinor period in minutes (default 1440, one day).
#' @param thr Dichotomisation threshold on the range-scaled cosine, in
#'   \eqn{[0, 1]} (default 0.2, approximating CircaCP's lower-20\% sleep cut): the
#'   fitted curve above \code{thr} is the rough active span, below it the rough
#'   rest span.
#' @param window_minutes Half-width of the search window, in minutes, in which
#'   each rough boundary is refined to a change point (default 240).
#'
#' @return An object of class \code{actiRhythm_changepoints}: the cosinor
#'   summary, a \code{changepoints} data frame (time and type, "sleep onset" or
#'   "wake onset"), a \code{sleep_episodes} data frame (sleep onset, wake onset,
#'   and duration in hours), and the mean sleep duration. The function never
#'   errors; on insufficient data it returns the structure with
#'   \code{insufficient = TRUE}.
#'
#' @references
#' \insertRef{chensun2024}{actiRhythm}
#'
#' @seealso \code{\link{state.transitions}}, \code{\link{sleep.regularity.index}}
#'
#' @examples
#' # Five days with a clear active day / restful night
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' counts <- ifelse(h >= 8 & h < 23, 300, 5) + pmax(0, stats::rnorm(length(ts), 0, 5))
#' sleep.changepoints(counts, ts)
#'
#' @export
sleep.changepoints <- function(counts, timestamps, period = 1440, thr = 0.2,
                               window_minutes = 240) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  if (!inherits(timestamps, "POSIXct"))
    stop("timestamps must be POSIXct")

  na_out <- function() structure(list(
    n = 0L, span_days = NA_real_,
    cosinor = c(mesor = NA_real_, amplitude = NA_real_, acrophase_hours = NA_real_),
    threshold = thr,
    changepoints = data.frame(time = as.POSIXct(character(0)), type = character(0)),
    sleep_episodes = data.frame(sleep_onset = as.POSIXct(character(0)),
                                wake_onset = as.POSIXct(character(0)),
                                duration_hours = numeric(0)),
    n_episodes = 0L, mean_sleep_duration = NA_real_, insufficient = TRUE),
    class = c("actiRhythm_changepoints", "list"))

  ts <- timestamps
  y  <- suppressWarnings(as.numeric(counts))
  keep <- is.finite(y) & !is.na(ts)
  ts <- ts[keep]; y <- y[keep]
  o <- order(ts); ts <- ts[o]; y <- y[o]
  n <- length(y)
  if (n < 2L * 1440L) return(na_out())
  span_days <- as.numeric(difftime(ts[n], ts[1], units = "days"))
  if (span_days < 2 || stats::sd(y) == 0) return(na_out())

  t_min <- as.numeric(difftime(ts, ts[1], units = "mins"))
  omega <- 2 * pi / period
  X <- cbind(1, cos(omega * t_min), sin(omega * t_min))
  fit <- tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(fit$coefficients))) return(na_out())
  b <- fit$coefficients
  amp <- sqrt(b[2]^2 + b[3]^2)
  acro_hours <- ((atan2(b[3], b[2]) / omega) %% period) / 60
  fitted <- as.numeric(X %*% b)

  rng <- range(fitted)
  if (diff(rng) == 0) return(na_out())
  Fn <- (fitted - rng[1]) / diff(rng)
  D <- as.integer(Fn > thr)          # 1 = active, 0 = rest
  E <- diff(D)
  B <- which(E != 0)                 # rough cosinor boundaries
  if (length(B) < 2L) return(na_out())

  refined <- integer(length(B)); types <- character(length(B))
  for (i in seq_along(B)) {
    lo <- max(1L, B[i] - window_minutes); hi <- min(n, B[i] + window_minutes)
    k <- .circa_cp_detect(y[lo:hi])
    refined[i] <- if (is.na(k)) B[i] else lo + k - 1L
    types[i] <- if (E[B[i]] == 1L) "wake onset" else "sleep onset"
  }
  ord <- order(refined)
  refined <- refined[ord]; types <- types[ord]
  dedup <- c(TRUE, diff(refined) > 0)
  refined <- refined[dedup]; types <- types[dedup]

  cps <- data.frame(time = ts[refined], type = types, stringsAsFactors = FALSE)

  # Pair each sleep onset with the next wake onset into a sleep episode.
  episodes <- data.frame(sleep_onset = as.POSIXct(character(0)),
                         wake_onset = as.POSIXct(character(0)),
                         duration_hours = numeric(0))
  sot_idx <- which(types == "sleep onset")
  for (si in sot_idx) {
    wi <- which(types == "wake onset" & refined > refined[si])
    if (!length(wi)) next
    wi <- wi[1]
    dur <- as.numeric(difftime(ts[refined[wi]], ts[refined[si]], units = "hours"))
    episodes <- rbind(episodes, data.frame(
      sleep_onset = ts[refined[si]], wake_onset = ts[refined[wi]],
      duration_hours = dur))
  }

  structure(list(
    n = n, span_days = span_days,
    cosinor = c(mesor = unname(b[1]), amplitude = unname(amp),
                acrophase_hours = unname(acro_hours)),
    threshold = thr, changepoints = cps, sleep_episodes = episodes,
    n_episodes = nrow(episodes),
    mean_sleep_duration = if (nrow(episodes)) mean(episodes$duration_hours) else NA_real_,
    insufficient = FALSE), class = c("actiRhythm_changepoints", "list"))
}


#' @export
print.actiRhythm_changepoints <- function(x, ...) {
  cat("Change-Point Sleep/Wake Detection\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data for detection\n\n"); return(invisible(x)) }
  cat(sprintf("  Span:           %.1f days (%d epochs)\n", x$span_days, x$n))
  cat(sprintf("  Cosinor acrophase: %.1f h\n", x$cosinor["acrophase_hours"]))
  cat(sprintf("  Change points:  %d (%d sleep episodes)\n",
              nrow(x$changepoints), x$n_episodes))
  if (x$n_episodes) {
    cat(sprintf("  Mean sleep duration: %.1f h\n", x$mean_sleep_duration))
    cat("\n  First sleep episodes:\n")
    e <- utils::head(x$sleep_episodes, 3)
    for (i in seq_len(nrow(e)))
      cat(sprintf("    sleep %s  ->  wake %s  (%.1f h)\n",
                  format(e$sleep_onset[i], "%m-%d %H:%M"),
                  format(e$wake_onset[i], "%m-%d %H:%M"), e$duration_hours[i]))
  }
  cat("\n  Reference: Chen and Sun (2024)\n\n")
  invisible(x)
}
