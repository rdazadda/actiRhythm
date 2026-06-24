# Raw-acceleration epoch metrics (GGIR-style), computed on the native-rate
# gravity-preserving signal from .raw_xyz(): ENMO, MAD and the z-angle. These are
# the raw-NATIVE inputs the agcounts band-pass filter cannot produce.

#' Auto-Calibrate Raw Acceleration to the Unit Gravity Sphere
#'
#' Estimates per-axis gain and offset corrections by the van Hees et al. (2014)
#' method: during non-movement windows the acceleration vector should lie on the
#' 1 g sphere, so gain/offset are fit by iteratively projecting non-movement
#' window means onto the closest point of the unit sphere. Apply as
#' \code{(raw - offset) * scale}.
#'
#' @param xyz A data frame or matrix of raw acceleration in g (columns x, y, z).
#' @param fs Sample rate in Hz.
#' @param sphere_crit Minimum coverage (g) each axis must span on both sides of
#'   zero for a stable fit (default 0.3).
#' @param sd_crit Per-axis rolling SD (g) below which a 10-second window counts as
#'   non-movement (default 0.013).
#' @param max_iter Maximum refinement iterations (default 1000).
#' @param tol Convergence tolerance on the calibration error (default 1e-9).
#'
#' @return A list with \code{scale} and \code{offset} (length-3), the calibration
#'   error before and after (\code{cal_error_start}, \code{cal_error_end}, mean
#'   absolute deviation from 1 g), \code{npoints}, and a \code{calibrated} flag.
#'   When there is too little non-movement data the identity correction is
#'   returned with \code{calibrated = FALSE}.
#'
#' @references
#' van Hees VT, et al. (2014). Autocalibration of accelerometer data for
#' free-living physical activity assessment using local gravity and temperature.
#' \emph{Journal of Applied Physiology}, 117(7):738-744.
#' \doi{10.1152/japplphysiol.00421.2014}
#'
#' @examples
#' # Recover a known per-axis gain and offset from non-movement windows
#' set.seed(1)
#' u <- matrix(rnorm(40 * 3), 40, 3); u <- u / sqrt(rowSums(u^2))   # sphere directions
#' raw <- do.call(rbind, lapply(seq_len(40), function(i)
#'   matrix(rep(u[i, ] / c(1.03, 0.97, 1.01) + c(0.04, -0.03, 0.02), each = 300),
#'          300, 3) + rnorm(900, 0, 0.004)))
#' auto.calibrate(data.frame(x = raw[, 1], y = raw[, 2], z = raw[, 3]), fs = 30)$scale
#'
#' @export
auto.calibrate <- function(xyz, fs, sphere_crit = 0.3, sd_crit = 0.013,
                           max_iter = 1000, tol = 1e-9) {
  X <- as.matrix(xyz[, 1:3]); n <- nrow(X)
  ident <- function(reason) list(scale = c(1, 1, 1), offset = c(0, 0, 0),
    cal_error_start = NA_real_, cal_error_end = NA_real_, npoints = 0L,
    calibrated = FALSE, reason = reason)
  win <- max(1L, round(fs * 10)); nw <- n %/% win
  if (nw < 10L) return(ident("insufficient data for calibration"))
  g  <- rep(seq_len(nw), each = win)[seq_len(nw * win)]
  Xt <- X[seq_len(nw * win), , drop = FALSE]
  mean_w <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, mean)), numeric(nw))
  sd_w   <- vapply(1:3, function(a) as.numeric(tapply(Xt[, a], g, stats::sd)), numeric(nw))
  nomov  <- rowSums(sd_w < sd_crit) == 3L & rowSums(abs(mean_w) < 2) == 3L
  D <- mean_w[nomov, , drop = FALSE]
  rad0 <- sqrt(rowSums(D^2))                            # keep only near-gravity points;
  D <- D[is.finite(rad0) & rad0 > 0.3 & rad0 < 1.7, , drop = FALSE]  # drop idle-sleep zeros
  if (nrow(D) < 10L) return(ident("too few non-movement windows"))
  if (any(apply(D, 2, min) > -sphere_crit) || any(apply(D, 2, max) < sphere_crit))
    return(ident("insufficient sphere coverage on all axes"))

  scale <- c(1, 1, 1); offset <- c(0, 0, 0)
  cal_start <- mean(abs(sqrt(rowSums(D^2)) - 1))
  cur <- D; prev <- Inf
  for (it in seq_len(max_iter)) {
    rad     <- pmax(sqrt(rowSums(cur^2)), 1e-8)
    closest <- cur / rad
    w       <- pmin(1 / pmax(abs(rad - 1), 0.01), 100)
    err     <- mean(abs(rad - 1))
    if (abs(prev - err) < tol) break
    prev <- err
    for (a in 1:3) {
      fit <- stats::lm.wfit(cbind(1, cur[, a]), closest[, a], w)
      ik <- unname(fit$coefficients[1]); sk <- unname(fit$coefficients[2])
      scale[a]  <- scale[a] * sk
      offset[a] <- offset[a] - ik / scale[a]
    }
    cur <- sweep(sweep(D, 2, offset, "-"), 2, scale, "*")
  }
  cal_end <- mean(abs(sqrt(rowSums(cur^2)) - 1))
  list(scale = scale, offset = offset, cal_error_start = cal_start,
       cal_error_end = cal_end, npoints = nrow(D),
       calibrated = is.finite(cal_end) && cal_end < 0.02, reason = "ok")
}

# Aggregate native-rate calibrated g (matrix x/y/z) to per-epoch raw metrics.
.raw_metrics_compute <- function(X, fs, epoch, metrics, start, tz) {
  n <- nrow(X); spe <- round(fs * epoch); ne <- n %/% spe
  if (ne < 1L) stop("recording shorter than one ", epoch, "s epoch", call. = FALSE)
  g  <- rep(seq_len(ne), each = spe)[seq_len(ne * spe)]
  Xt <- X[seq_len(ne * spe), , drop = FALSE]
  en <- sqrt(rowSums(Xt^2))                       # per-sample Euclidean norm (g)
  out <- data.frame(time = seq(start, by = epoch, length.out = ne))
  if ("ENMO" %in% metrics)
    out$ENMO <- as.numeric(tapply(pmax(en - 1, 0), g, mean)) * 1000   # mg
  if ("MAD" %in% metrics) {
    em <- as.numeric(tapply(en, g, mean))[g]
    out$MAD <- as.numeric(tapply(abs(en - em), g, mean)) * 1000       # mg
  }
  if ("anglez" %in% metrics) {
    k <- round(fs * 5); if (k %% 2 == 0) k <- k + 1L                  # odd window
    x5 <- stats::runmed(Xt[, 1], k, endrule = "median")
    y5 <- stats::runmed(Xt[, 2], k, endrule = "median")
    z5 <- stats::runmed(Xt[, 3], k, endrule = "median")
    out$anglez <- as.numeric(tapply(atan2(z5, sqrt(x5^2 + y5^2)) * 180 / pi, g, mean))
  }
  out
}

#' Raw-Acceleration Epoch Metrics (ENMO, MAD, z-Angle)
#'
#' Reads a raw accelerometer file and returns per-epoch raw activity and posture
#' metrics, the gravity-preserving signals that counts cannot represent: ENMO
#' (Euclidean Norm Minus One, a raw activity metric), MAD (Mean Amplitude
#' Deviation), and the z-angle (arm/posture angle). Auto-calibration
#' (van Hees 2014) is applied first by default. Requires the relevant raw reader
#' (\pkg{read.gt3x} for \code{.gt3x}, \pkg{GGIRread} for \code{.cwa}/\code{.bin}).
#'
#' @param x A path to a raw file (\code{.gt3x}, \code{.cwa}, \code{.bin}) or a raw
#'   data frame with \code{x}/\code{y}/\code{z} columns in g and an \code{fs}
#'   attribute (e.g. from \code{\link{example_raw}} or your own device).
#' @param device One of \code{"auto"}, \code{"gt3x"}, \code{"axivity"},
#'   \code{"geneactiv"} (default \code{"auto"}, inferred from the extension; used
#'   only when \code{x} is a file path).
#' @param epoch Epoch length in seconds (default 60).
#' @param metrics Which metrics to return; any of \code{"ENMO"}, \code{"MAD"},
#'   \code{"anglez"} (default all three).
#' @param calibrate Apply van Hees auto-calibration first (default \code{TRUE}).
#' @param tz Time zone for the timestamps (default \code{"UTC"}).
#'
#' @return A data frame with \code{time} and the requested metrics (ENMO and MAD
#'   in mg, anglez in degrees), one row per epoch. The calibration result is
#'   attached as the \code{"calibration"} attribute.
#'
#' @references
#' van Hees VT, et al. (2013). Estimation of daily energy expenditure in pregnant
#' and non-pregnant women using a wrist-worn tri-axial accelerometer.
#' \emph{PLoS ONE}, 8(4):e61691. \doi{10.1371/journal.pone.0061691}
#'
#' Vaha-Ypya H, et al. (2015). A universal, accurate intensity-based
#' classification of different physical activities using raw data of
#' accelerometer. \emph{Clinical Physiology and Functional Imaging}, 35(1):64-70.
#' \doi{10.1111/cpf.12127}
#'
#' @seealso \code{\link{auto.calibrate}}, \code{\link{circadian.raw}},
#'   \code{\link{rest.spt}}, \code{\link{example_raw}}
#'
#' @examples
#' # On a synthetic raw recording (no file needed); pass a path for a real file
#' \donttest{
#' m <- raw.metrics(example_raw(days = 1), epoch = 60)
#' head(m)
#' }
#'
#' @export
raw.metrics <- function(x, device = "auto", epoch = 60,
                        metrics = c("ENMO", "MAD", "anglez"),
                        calibrate = TRUE, tz = "UTC") {
  metrics <- match.arg(metrics, c("ENMO", "MAD", "anglez"), several.ok = TRUE)
  r <- if (is.character(x)) .raw_xyz(x, device, tz) else .as_raw_xyz(x)
  X <- cbind(r$x, r$y, r$z)
  cal <- NULL
  if (isTRUE(calibrate)) {
    cal <- auto.calibrate(data.frame(x = r$x, y = r$y, z = r$z), r$fs)
    if (isTRUE(cal$calibrated))
      X <- sweep(sweep(X, 2, cal$offset, "-"), 2, cal$scale, "*")
  }
  out <- .raw_metrics_compute(X, r$fs, epoch, metrics, r$start, tz)
  attr(out, "calibration") <- cal
  out
}

# Extract x/y/z + sample rate + start from a raw acceleration data frame (the
# non-file input accepted by raw.metrics / circadian.raw).
.as_raw_xyz <- function(df) {
  if (!is.data.frame(df))
    stop("raw input must be a file path or a data frame with x/y/z columns", call. = FALSE)
  nm <- tolower(names(df))
  pick <- function(k) {
    i <- which(nm == k)
    if (!length(i)) stop("raw data frame needs a '", k, "' column", call. = FALSE)
    df[[i[1]]]
  }
  tcol <- if ("time" %in% nm) df[[which(nm == "time")[1]]] else NULL
  fs <- attr(df, "fs")
  if (is.null(fs)) {
    if (is.null(tcol)) stop("raw data frame needs an 'fs' attribute or a 'time' column", call. = FALSE)
    fs <- 1 / as.numeric(stats::median(diff(as.numeric(tcol))))
  }
  start <- if (!is.null(tcol)) as.POSIXct(tcol[1]) else as.POSIXct("2024-01-01 12:00", tz = "UTC")
  list(x = pick("x"), y = pick("y"), z = pick("z"), fs = fs, start = start)
}

# Deterministic pseudo-noise in [-0.5, 0.5] (fract-sin hash), so example_raw is
# reproducible without touching the global RNG state.
.det_noise <- function(n, k) {
  v <- sin((seq_len(n) + k * 7919) * 12.9898) * 43758.5453
  (v - floor(v)) - 0.5
}

#' Synthetic Raw Acceleration Recording
#'
#' Generates a deterministic synthetic triaxial raw acceleration recording (in g)
#' with a day/night posture cycle, daytime movement and posture changes, still
#' (sleeping) nights, and a slight built-in miscalibration. It is the file-free
#' stand-in for a real raw file in the examples: feed it to
#' \code{\link{raw.metrics}}, \code{\link{circadian.raw}}, or the z-angle sleep
#' pipeline. Raw acceleration is far too large to ship as data, so it is generated
#' on demand rather than bundled.
#'
#' @param days Recording length in days (default 2, giving two nights).
#' @param fs Sample rate in Hz (default 30).
#' @param device_off Days of a still, taken-off device to append at the end
#'   (default 0). Use it to exercise \code{\link{detect.nonwear.raw}} and the
#'   non-wear gate of \code{\link{rest.spt}}.
#'
#' @return A data frame with a POSIXct \code{time} column and \code{x}/\code{y}/
#'   \code{z} acceleration in g, with the sample rate in the \code{"fs"} attribute.
#'
#' @seealso \code{\link{raw.metrics}}, \code{\link{circadian.raw}},
#'   \code{\link{rest.spt}}
#'
#' @examples
#' raw <- example_raw(days = 1)
#' str(raw)
#'
#' @export
example_raw <- function(days = 2, fs = 30, device_off = 0) {
  n  <- as.integer(round(days * 86400 * fs))
  start <- as.POSIXct("2024-01-01 12:00", tz = "UTC")
  ti <- as.numeric(start) + (seq_len(n) - 1) / fs
  hr <- (ti %% 86400) / 3600                              # hour of day (UTC)
  night <- hr >= 23 | hr < 7
  wob <- ifelse(night, 0, 0.35 * sin(2 * pi * ti / 40))   # daytime posture changes
  gx <- ifelse(night, 0.20, 0.10)
  gy <- ifelse(night, 0.30, -0.95) + wob
  gz <- ifelse(night, -0.93, 0.25) - wob
  gm <- sqrt(gx^2 + gy^2 + gz^2); gx <- gx / gm; gy <- gy / gm; gz <- gz / gm
  mv <- ifelse(night, 0.05, 0.10)                         # sleep keeps micro-movement (worn)
  x <- gx + mv * .det_noise(n, 1)
  y <- gy + mv * .det_noise(n, 2)
  z <- gz + mv * .det_noise(n, 3)
  if (device_off > 0) {                                   # still, taken-off device
    no <- as.integer(round(device_off * 86400 * fs))
    ti <- c(ti, ti[n] + seq_len(no) / fs)
    x  <- c(x, 0.02  + 5e-4 * .det_noise(no, 4))
    y  <- c(y, 0.05  + 5e-4 * .det_noise(no, 5))
    z  <- c(z, 0.998 + 5e-4 * .det_noise(no, 6))
  }
  x <- x / 1.02 + 0.03; y <- y / 0.98 - 0.02; z <- z / 1.01 + 0.015  # slight miscalibration
  df <- data.frame(time = as.POSIXct(ti, origin = "1970-01-01", tz = "UTC"),
                   x = x, y = y, z = z)
  attr(df, "fs") <- fs
  df
}
