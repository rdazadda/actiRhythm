# FFT-based Morlet continuous wavelet transform (Torrence & Compo 1998).
# Returns the complex transform, scales, and Fourier periods (in epochs).
.cwt_morlet <- function(x, dj, omega0) {
  N <- length(x)
  ff <- 4 * pi / (omega0 + sqrt(2 + omega0^2))   # Fourier factor (Torrence-Compo)
  s0 <- 2 / ff                                    # smallest scale about 2 epochs
  smax <- (N / 3) / ff
  J <- max(1L, floor(log2(smax / s0) / dj))
  scales <- s0 * 2^((0:J) * dj)
  xhat <- stats::fft(x)
  k <- c(0:floor(N / 2), -(rev(seq_len(ceiling(N / 2) - 1))))
  w <- 2 * pi * k / N
  W <- matrix(0 + 0i, length(scales), N)
  for (i in seq_along(scales)) {
    s <- scales[i]
    norm <- sqrt(2 * pi * s) * pi^(-0.25)
    daughter <- norm * exp(-0.5 * (s * w - omega0)^2) * (w > 0)
    W[i, ] <- stats::fft(xhat * daughter, inverse = TRUE) / N
  }
  list(W = W, scales = scales, period = ff * scales, ff = ff)
}

# Torrence-Compo (1998) cone of influence: the largest reliable Fourier period
# (hours) at each time column, rising from a sliver at the edges to N/2 in the
# middle. Power at periods above this curve is edge-contaminated.
.torrence_coi <- function(N, ff, ep_h) {
  d <- pmin(seq_len(N) - 1L, N - seq_len(N))
  d[d == 0] <- 1e-5
  (ff / sqrt(2)) * d * ep_h
}

#' Continuous Wavelet Transform of the Activity Rhythm
#'
#' Runs a Morlet continuous wavelet transform on the activity series and returns
#' the time-frequency power surface across circadian and ultradian periods and a
#' dominant-period-over-time track (Torrence and Compo 1998; Leise 2013). Unlike
#' the sliding chi-square spectrogram, it localizes period drift at every time
#' point, so a lengthening or fragmenting rhythm shows up directly.
#'
#' @param counts Numeric activity vector (a coarse epoch, e.g. 10-minute bins, is
#'   recommended for speed; see the example).
#' @param timestamps POSIXct timestamps, one per value.
#' @param dj Scale resolution in voices per octave step (default 1/12).
#' @param omega0 Morlet central frequency (default 6).
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return An object of class \code{actiRhythm_wavelet}: the period grid (hours),
#'   the time-averaged global power spectrum, the per-time dominant period, the
#'   overall peak period, the power matrix, and the cone of influence
#'   (\code{coi_period_h}, the largest reliable period at each time, with the
#'   \code{in_coi} mask of edge-affected cells). The global spectrum and dominant
#'   period are computed outside the cone so edge effects do not bias them. A
#'   \code{significant} logical matrix flags cells whose power exceeds the 95\%
#'   confidence level against an AR(1) red-noise background (with the per-scale
#'   threshold \code{sig_power} and the lag-1 autocorrelation \code{phi}), which
#'   separates a real rhythm from background. Never errors.
#'
#' @references
#' \insertRef{torrence1998}{actiRhythm}
#'
#' \insertRef{leise2013}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' circadian.wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#'
#' @export
circadian.wavelet <- function(counts, timestamps, dj = 1 / 12, omega0 = 6,
                              epoch_length = 60) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  o <- order(ts); x <- x[o]; ts <- ts[o]
  N <- length(x)
  ep_h <- epoch_length / 3600
  na_out <- structure(list(period_hours = numeric(0), global_power = numeric(0),
    dominant_period = numeric(0), peak_period_h = NA_real_, power = matrix(nrow = 0, ncol = 0),
    coi_period_h = numeric(0), coi_multiplier = NA_real_, in_coi = matrix(nrow = 0, ncol = 0),
    significant = matrix(nrow = 0, ncol = 0), sig_power = numeric(0), phi = NA_real_,
    times = ts, insufficient = TRUE), class = c("actiRhythm_wavelet", "list"))
  if (N < 64L || stats::sd(x) == 0) return(na_out)

  xc <- x - mean(x)
  cw <- .cwt_morlet(xc, dj, omega0)
  period_h <- cw$period * ep_h
  power <- (Mod(cw$W)^2) / cw$scales            # divide-by-scale bias rectification
  coi_h <- .torrence_coi(N, cw$ff, ep_h)        # cone of influence (hours, per time)
  in_coi <- outer(period_h, coi_h, ">")         # TRUE where the cell is edge-affected
  masked <- power; masked[in_coi] <- NA
  global <- rowMeans(masked, na.rm = TRUE)
  dom <- period_h[apply(masked, 2, function(col)
    if (all(is.na(col))) NA_integer_ else which.max(col))]
  peak <- period_h[which.max(global)]
  # Torrence-Compo (1998) red-noise significance: power above the background
  # spectrum of an AR(1) process at 95%, counted only outside the cone of influence.
  phi <- max(0, min(stats::acf(xc, lag.max = 1, plot = FALSE)$acf[2], 0.99))
  Pk  <- (1 - phi^2) / (1 + phi^2 - 2 * phi * cos(2 * pi / cw$period))
  sig_power   <- stats::var(xc) * Pk * stats::qchisq(0.95, 2) / 2 / cw$scales
  significant <- sweep(power, 1, sig_power, ">") & !in_coi

  structure(list(period_hours = period_h, global_power = global,
    dominant_period = dom, peak_period_h = peak, power = power,
    coi_period_h = coi_h, coi_multiplier = cw$ff / sqrt(2), in_coi = in_coi,
    significant = significant, sig_power = sig_power, phi = phi,
    times = ts, insufficient = FALSE), class = c("actiRhythm_wavelet", "list"))
}

#' @export
print.actiRhythm_wavelet <- function(x, ...) {
  cat("Continuous Wavelet Transform (Morlet)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Periods scanned: %.1f to %.1f h (%d scales)\n",
              min(x$period_hours), max(x$period_hours), length(x$period_hours)))
  cat(sprintf("  Peak global period: %.2f h\n", x$peak_period_h))
  cat(sprintf("  Dominant period: median %.2f h (IQR %.2f-%.2f)\n",
              stats::median(x$dominant_period, na.rm = TRUE),
              stats::quantile(x$dominant_period, 0.25, na.rm = TRUE),
              stats::quantile(x$dominant_period, 0.75, na.rm = TRUE)))
  cat(sprintf("  Cone of influence: %.0f%% of cells edge-reliable\n",
              100 * mean(!x$in_coi)))
  cat(sprintf("  Significant power: %.1f%% of reliable cells (95%% vs AR(1) red noise, phi=%.2f)\n\n",
              100 * sum(x$significant) / max(1, sum(!x$in_coi)), x$phi))
  invisible(x)
}

#' Wavelet Cone of Influence
#'
#' The Morlet cone of influence (Torrence and Compo 1998): the largest period, in
#' hours, at each time point below which the wavelet power is free of edge
#' effects. Power at periods above this curve, near the start and end of the
#' recording, is unreliable. \code{\link{circadian.wavelet}} returns this curve in
#' its result; this helper computes it directly for plotting or masking.
#'
#' @param timestamps POSIXct timestamps, one per epoch.
#' @param omega0 Morlet central frequency (default 6).
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return Numeric vector, one per epoch, of the maximum reliable period (hours).
#'
#' @references
#' \insertRef{torrence1998}{actiRhythm}
#'
#' @seealso \code{\link{circadian.wavelet}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' range(wavelet.coi(ts, epoch_length = 600))
#'
#' @export
wavelet.coi <- function(timestamps, omega0 = 6, epoch_length = 60) {
  ff <- 4 * pi / (omega0 + sqrt(2 + omega0^2))
  .torrence_coi(length(timestamps), ff, epoch_length / 3600)
}


# One level of a circular a-trous (undecimated, MODWT) Haar filter.
.modwt_haar_level <- function(v, j) {
  N <- length(v); shift <- 2^(j - 1L)
  vp <- v[((seq_len(N) - 1L - shift) %% N) + 1L]   # circularly shift by 2^(j-1)
  list(W = (v - vp) / 2, V = (v + vp) / 2)
}

#' Ultradian Wavelet Band Power
#'
#' Partitions the activity variance into dyadic period bands with an undecimated
#' (MODWT) Haar wavelet transform, isolating ultradian bands such as the
#' about-90-minute, about-4-hour, and about-8-hour rhythms (Percival and Walden
#' 2000; Leise 2013). Reports the energy, power, and fraction of total variance
#' in each band.
#'
#' @param counts Numeric activity vector (minute epochs recommended).
#' @param timestamps POSIXct timestamps, one per value.
#' @param bands Named list of \code{c(low_hours, high_hours)} period bands
#'   (default about 90 min, about 4 h, about 8 h).
#' @param epoch_length Epoch length in seconds (default 60).
#'
#' @return An object of class \code{actiRhythm_bandpower}: a per-band table
#'   (energy, power, fraction) and the per-level wavelet variance. Never errors.
#'
#' @references
#' \insertRef{percival2000}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' ultradian.bandpower(100 + 50 * sin(2 * pi * th / 4), ts)   # a 4-hour rhythm
#'
#' @export
ultradian.bandpower <- function(counts, timestamps,
                                bands = list(`90min` = c(1, 2), `4h` = c(2, 6),
                                             `8h` = c(6, 12)),
                                epoch_length = 60) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  x[!is.finite(x)] <- 0
  ep_h <- epoch_length / 3600
  N <- length(x)
  na_out <- structure(list(bands = data.frame(), levels = data.frame(),
    insufficient = TRUE), class = c("actiRhythm_bandpower", "list"))
  if (N < 16L || stats::sd(x) == 0) return(na_out)

  xc <- x - mean(x)
  J <- min(14L, floor(log2(N)) - 1L)
  v <- xc; energy <- numeric(J)
  for (j in seq_len(J)) {
    d <- .modwt_haar_level(v, j)
    energy[j] <- sum(d$W^2); v <- d$V
  }
  smooth_energy <- sum(v^2)
  total <- sum(xc^2)
  # level j is the Fourier period band [2^j, 2^{j+1}] epochs
  lo_h <- 2^(seq_len(J)) * ep_h; hi_h <- 2^(seq_len(J) + 1L) * ep_h
  levels <- data.frame(level = seq_len(J), period_lo_h = lo_h, period_hi_h = hi_h,
                       energy = energy, variance = energy / N)

  brow <- lapply(names(bands), function(nm) {
    b <- bands[[nm]]
    sel <- hi_h > b[1] & lo_h < b[2]            # levels overlapping the band
    e <- sum(energy[sel])
    data.frame(band = nm, low_h = b[1], high_h = b[2], energy = e,
               power = e / N, fraction = e / total)
  })
  structure(list(bands = do.call(rbind, brow), levels = levels,
    smooth_fraction = smooth_energy / total, insufficient = FALSE),
    class = c("actiRhythm_bandpower", "list"))
}

#' @export
print.actiRhythm_bandpower <- function(x, ...) {
  cat("Ultradian Wavelet Band Power\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  b <- x$bands; b$fraction <- round(b$fraction, 3); b$power <- round(b$power, 1)
  print(b[c("band", "low_h", "high_h", "power", "fraction")], row.names = FALSE)
  cat("\n")
  invisible(x)
}
