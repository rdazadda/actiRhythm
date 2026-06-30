.circ_mean_h <- function(h, period = 24) {
  a <- 2 * pi * (h %% period) / period
  (atan2(mean(sin(a), na.rm = TRUE), mean(cos(a), na.rm = TRUE)) %% (2 * pi)) / (2 * pi) * period
}
.circ_sd_h <- function(h, period = 24) {
  a <- 2 * pi * (h %% period) / period
  R <- sqrt(mean(cos(a), na.rm = TRUE)^2 + mean(sin(a), na.rm = TRUE)^2)
  sqrt(-2 * log(R)) / (2 * pi) * period
}

#' Curve Registration of Daily Activity Profiles
#'
#' Aligns each day's 24-hour activity profile on its active-phase landmark (the
#' M10 centre), separating the horizontal phase variation (how the timing shifts
#' day to day) from the vertical amplitude variation (the registered mean profile,
#' sharper than the plain average because phase jitter no longer blurs it). The
#' per-day landmark times are a scale-invariant phase marker (the M10-window
#' centre), unchanged by any rescaling of the counts (Ramsay and Silverman 2005).
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param n_grid Bins per day for the profile (default 144 = 10-minute bins).
#' @param period Period in hours (default 24).
#'
#' @return An object of class \code{actiRhythm_registration}: a per-day landmark
#'   table (L5 and M10 centre hours), the circular-mean landmarks, the phase
#'   variability (circular SD of the M10 landmark), and the registered mean
#'   profile. Never errors.
#'
#' @references
#' \insertRef{ramsay2005}{actiRhythm}
#'
#' \insertRef{vansomeren1999}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' curve.registration(ifelse(h >= 23 | h < 7, 5, 300), ts)
#'
#' @export
curve.registration <- function(counts, timestamps, n_grid = 144L, period = 24) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  day <- as.Date(ts)
  binw <- 1440 / n_grid
  k5 <- round(5 * 60 / binw); k10 <- round(10 * 60 / binw)

  profs <- list(); L5c <- numeric(0); M10c <- numeric(0); dts <- character(0)
  udays <- sort(unique(day))
  for (di in seq_along(udays)) {
    d <- udays[di]
    idx <- day == d; xd <- x[idx]; td <- ts[idx]
    if (sum(is.finite(xd)) < n_grid * binw * 0.3) next
    mod <- as.POSIXlt(td)$hour * 60 + as.POSIXlt(td)$min
    g <- pmin(floor(mod / binw) + 1L, n_grid)
    p <- as.numeric(tapply(xd, factor(g, levels = seq_len(n_grid)), mean, na.rm = TRUE))
    if (mean(is.na(p)) > 0.5) next
    p[is.na(p)] <- mean(p, na.rm = TRUE)
    cs <- cumsum(c(0, c(p, p)))
    l5 <- cs[(k5 + 1):(k5 + n_grid)] - cs[1:n_grid]
    m10 <- cs[(k10 + 1):(k10 + n_grid)] - cs[1:n_grid]
    L5c <- c(L5c, ((which.min(l5) - 1) * binw / 60 + 2.5) %% 24)
    M10c <- c(M10c, ((which.max(m10) - 1) * binw / 60 + 5) %% 24)
    profs[[length(profs) + 1L]] <- p; dts <- c(dts, as.character(d))
  }
  na_out <- structure(list(landmarks = data.frame(), mean_L5 = NA_real_,
    mean_M10 = NA_real_, phase_sd = NA_real_, registered_mean = numeric(0),
    n_days = 0L, insufficient = TRUE), class = c("actiRhythm_registration", "list"))
  if (length(profs) < 2L) return(na_out)

  target <- .circ_mean_h(M10c, period)
  # circular-shift each day so its M10 lands on the target, then average
  reg <- vapply(seq_along(profs), function(i) {
    shift <- round((M10c[i] - target) * 60 / binw)
    p <- profs[[i]]
    p[((seq_len(n_grid) - 1L + shift) %% n_grid) + 1L]
  }, numeric(n_grid))
  registered_mean <- rowMeans(reg)

  structure(list(
    landmarks = data.frame(date = as.Date(dts), L5_center_h = L5c, M10_center_h = M10c),
    mean_L5 = .circ_mean_h(L5c, period), mean_M10 = target,
    phase_sd = .circ_sd_h(M10c, period), registered_mean = registered_mean,
    n_days = length(profs), insufficient = FALSE),
    class = c("actiRhythm_registration", "list"))
}

#' @export
print.actiRhythm_registration <- function(x, ...) {
  cat("Curve Registration (landmark)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data (need >= 2 days)\n\n"); return(invisible(x)) }
  cat(sprintf("  Days registered:   %d\n", x$n_days))
  cat(sprintf("  Mean L5 centre:    %05.2f h\n", x$mean_L5))
  cat(sprintf("  Mean M10 centre:   %05.2f h\n", x$mean_M10))
  cat(sprintf("  Phase variability: %.2f h (circular SD of M10)\n\n", x$phase_sd))
  invisible(x)
}


#' Residual Circadian Spectrum
#'
#' Removes the fitted cosinor mean and estimates the spectrum of what is left, the
#' residual circadian spectrum, integrating it into frequency bands. It measures
#' the ultradian and noise structure that the cosinor does not fit (Krafty et
#' al. 2019), so two recordings with the same 24-hour rhythm but different
#' residual fragmentation are told apart.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Cosinor period in hours (default 24).
#' @param bands Named list of \code{c(low_hours, high_hours)} period bands for the
#'   residual power; defaults to ultradian (2-8 h) and high-frequency (0.5-2 h).
#'
#' @return An object of class \code{actiRhythm_rcs}: the residual variance, a
#'   per-band power table, and the spectrum. Never errors.
#'
#' @references
#' \insertRef{krafty2019}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' residual.spectrum(100 + 80 * cos(2 * pi * th / 24) + 20 * sin(2 * pi * th / 4), ts,
#'                   period = 24)
#'
#' @export
residual.spectrum <- function(counts, timestamps, period = 24,
                              bands = list(ultradian = c(2, 8), high_freq = c(0.5, 2))) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  o <- order(ts); x <- x[o]; ts <- ts[o]
  n <- length(x)
  dt_h <- as.numeric(stats::median(diff(as.numeric(ts)))) / 3600
  na_out <- structure(list(residual_var = NA_real_, bands = data.frame(),
    spectrum = data.frame(), insufficient = TRUE), class = c("actiRhythm_rcs", "list"))
  if (n < 32L || stats::sd(x) == 0 || !is.finite(dt_h) || dt_h <= 0) return(na_out)

  th <- as.numeric(difftime(ts, ts[1], units = "hours"))
  X <- cbind(1, cos(2 * pi * th / period), sin(2 * pi * th / period))
  b <- stats::lm.fit(X, x)$coefficients
  resid <- x - as.numeric(X %*% b)

  pg <- stats::spec.pgram(resid, taper = 0, detrend = FALSE, plot = FALSE)
  freq_per_h <- pg$freq / dt_h                 # cycles per hour
  period_h <- 1 / freq_per_h
  spec <- pg$spec
  total <- sum(spec)

  brow <- lapply(names(bands), function(nm) {
    bd <- bands[[nm]]; sel <- period_h >= bd[1] & period_h <= bd[2]
    p <- if (any(sel)) sum(spec[sel]) else 0
    data.frame(band = nm, low_h = bd[1], high_h = bd[2], power = p,
               fraction = if (total > 0) p / total else NA_real_)
  })
  structure(list(residual_var = stats::var(resid), bands = do.call(rbind, brow),
    spectrum = data.frame(period_h = period_h, power = spec), insufficient = FALSE),
    class = c("actiRhythm_rcs", "list"))
}

#' @export
print.actiRhythm_rcs <- function(x, ...) {
  cat("Residual Circadian Spectrum\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Residual variance: %.1f\n\n", x$residual_var))
  b <- x$bands; b$fraction <- round(b$fraction, 3); b$power <- round(b$power, 1)
  print(b, row.names = FALSE)
  cat("\n")
  invisible(x)
}
