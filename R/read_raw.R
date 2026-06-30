# Multi-brand raw accelerometer ingest. Each brand reader normalizes its file to
# a triaxial acceleration frame (g), which is converted to ActiGraph-equivalent
# activity counts with the agcounts implementation of the count algorithm (Neishabouri 2022).
# Cross-brand counts are an APPROXIMATION, not native ActiGraph output. See the
# exported functions' documentation.

# Convert a raw frame (a data.frame with a POSIXct `time` column + X/Y/Z in g, at
# a sample rate agcounts accepts) to time/axis1/axis2/axis3/vm counts.
.raw_to_counts <- function(raw, epoch, lfe, tz) {
  dur <- as.numeric(difftime(raw$time[nrow(raw)], raw$time[1], units = "secs"))
  if (!is.finite(dur) || dur < epoch)
    stop("recording (", round(dur, 1), "s) is shorter than one ", epoch,
         "s epoch", call. = FALSE)
  cc <- agcounts::calculate_counts(raw, epoch = epoch, lfe_select = lfe,
                                   tz = tz, verbose = FALSE)
  data.frame(time = cc$time, axis1 = cc$Axis1, axis2 = cc$Axis2,
             axis3 = cc$Axis3, vm = cc$Vector.Magnitude)
}

# Sample rates agcounts accepts natively; any other rate is resampled to 30 Hz.
.agcounts_rates <- c(30, 40, 50, 60, 70, 80, 90, 100)

# Linearly resample triaxial g from fs to 30 Hz (for non-decade rates such as the
# 25 Hz of some Axivity configs or the 85.7 Hz of GENEActiv).
.resample_to_30hz <- function(xyz, fs) {
  n <- nrow(xyz); t_old <- (seq_len(n) - 1) / fs
  t_new <- seq(0, (n - 1) / fs, by = 1 / 30)
  data.frame(X = stats::approx(t_old, xyz$X, t_new, rule = 2)$y,
             Y = stats::approx(t_old, xyz$Y, t_new, rule = 2)$y,
             Z = stats::approx(t_old, xyz$Z, t_new, rule = 2)$y)
}

# Build the time/X/Y/Z frame agcounts wants from raw triaxial g + native sample
# rate + start time, resampling to 30 Hz when the rate is not one agcounts takes.
.build_raw <- function(xyz, fs, start, tz) {
  if (!any(abs(fs - .agcounts_rates) < 1e-6)) {
    xyz <- .resample_to_30hz(xyz, fs); fs <- 30
  }
  data.frame(time = seq(start, by = 1 / fs, length.out = nrow(xyz)),
             X = xyz$X, Y = xyz$Y, Z = xyz$Z)
}

# Infer the device brand from a file extension.
.detect_device <- function(path) {
  ext <- tolower(sub("^.*\\.", "", basename(path)))
  switch(ext, gt3x = "gt3x", cwa = "axivity", bin = "geneactiv",
    stop("cannot infer device from extension '", ext, "'; set 'device' explicitly"))
}

# Carry the last gravity vector forward through all-zero (idle-sleep) runs, the
# convention agcounts/GGIR use. read.gt3x imputes idle-sleep gaps as zeros; left
# as zeros they collapse the z-angle to the calibration-offset angle, so LOCF them.
# Vectorized: cummax of the last non-zero index, then index-assign. Leading zeros
# (no prior gravity) stay zero.
.locf_idle_sleep <- function(x, y, z) {
  zero <- x == 0 & y == 0 & z == 0
  if (!any(zero)) return(list(x = x, y = y, z = z))
  src <- seq_along(x); src[zero] <- 0L; src <- cummax(src)
  fill <- zero & src > 0L
  x[fill] <- x[src[fill]]; y[fill] <- y[src[fill]]; z[fill] <- z[src[fill]]
  list(x = x, y = y, z = z)
}

# Native-rate triaxial acceleration accessor (g), shared by the raw-metric and
# raw-sleep functions. Returns the unresampled x/y/z plus the native sample rate
# and the recording start time. It does NOT band-pass to counts, so the gravity
# / posture (DC) component that ENMO and the z-angle need is preserved.
.raw_xyz <- function(path, device = "auto", tz = "UTC") {
  if (device == "auto") device <- .detect_device(path)
  if (device == "gt3x") {
    if (!requireNamespace("read.gt3x", quietly = TRUE))
      stop("reading .gt3x requires the 'read.gt3x' package")
    mat <- read.gt3x::read.gt3x(path, asDataFrame = FALSE, imputeZeroes = TRUE)
    m  <- unclass(mat)
    t0 <- as.numeric(attr(mat, "start_time")) + attr(mat, "time_index")[1] / 100
    g  <- .locf_idle_sleep(m[, "X"], m[, "Y"], m[, "Z"])   # carry gravity through idle-sleep
    return(list(x = g$x, y = g$y, z = g$z, fs = attr(mat, "sample_rate"),
                start = as.POSIXct(t0, origin = "1970-01-01", tz = tz)))
  }
  if (!requireNamespace("GGIRread", quietly = TRUE))
    stop("reading raw ", device, " files requires the 'GGIRread' package")
  if (device == "axivity") {
    d <- GGIRread::readAxivity(path, start = 1, end = 1e6, desiredtz = tz)
    return(list(x = d$data$x, y = d$data$y, z = d$data$z, fs = d$header$frequency,
                start = as.POSIXct(d$data$time[1], origin = "1970-01-01", tz = tz)))
  }
  if (device == "geneactiv") {
    d  <- GGIRread::readGENEActiv(path, start = 1, end = 1e6, desiredtz = tz)
    dd <- if (!is.null(d$data.out)) d$data.out else d$data
    return(list(x = dd$x, y = dd$y, z = dd$z, fs = d$header$SampleRate,
                start = as.POSIXct(dd$time[1], origin = "1970-01-01", tz = tz)))
  }
  stop("unsupported device '", device, "'")
}

#' Activity Counts from a Raw Axivity .cwa File
#'
#' Reads a raw Axivity (\code{.cwa}) accelerometer file and converts it to
#' ActiGraph-equivalent activity counts via the agcounts implementation of the
#' ActiGraph count algorithm (Neishabouri 2022). Requires the \pkg{GGIRread} and
#' \pkg{agcounts} packages.
#'
#' @section Cross-brand counts: These are an \emph{approximation} of ActiGraph
#'   counts, not native ActiGraph output. Axivity-to-count conversion has been
#'   directly validated (Brond et al. 2017). The result is appropriate for the
#'   relative and normalized circadian metrics in this package (IS, IV, RA, L5,
#'   M10, SRI); it should \emph{not} be used to apply ActiGraph intensity
#'   cut-points or to compare absolute counts across device brands.
#'
#' @param path Path to a \code{.cwa} file.
#' @param epoch Epoch length in seconds (default 60).
#' @param lfe Use the low-frequency extension filter (default \code{FALSE}).
#' @param tz Time zone for the timestamps (default \code{"UTC"}).
#'
#' @return Data frame with \code{time}, \code{axis1}, \code{axis2}, \code{axis3}
#'   and \code{vm}, one row per epoch (the same shape as \code{\link{gt3x.counts}}).
#'
#' @references
#' \insertRef{neishabouri2022}{actiRhythm}
#'
#' \insertRef{brond2017}{actiRhythm}
#'
#' @seealso \code{\link{geneactiv.counts}}, \code{\link{gt3x.counts}},
#'   \code{\link{read.raw}}
#'
#' @export
axivity.counts <- function(path, epoch = 60, lfe = FALSE, tz = "UTC") {
  if (!requireNamespace("GGIRread", quietly = TRUE))
    stop("axivity.counts() requires the 'GGIRread' package: install.packages('GGIRread')")
  if (!requireNamespace("agcounts", quietly = TRUE))
    stop("axivity.counts() requires the 'agcounts' package: install.packages('agcounts')")
  d  <- GGIRread::readAxivity(path, start = 1, end = 1e6, desiredtz = tz)
  fs <- d$header$frequency
  start <- as.POSIXct(d$data$time[1], origin = "1970-01-01", tz = tz)
  raw <- .build_raw(data.frame(X = d$data$x, Y = d$data$y, Z = d$data$z), fs, start, tz)
  .raw_to_counts(raw, epoch, lfe, tz)
}

#' Activity Counts from a Raw GENEActiv .bin File
#'
#' Reads a raw GENEActiv (\code{.bin}) accelerometer file and converts it to
#' ActiGraph-equivalent activity counts via the agcounts implementation of the
#' ActiGraph count algorithm (Neishabouri 2022). Requires the \pkg{GGIRread} and
#' \pkg{agcounts} packages.
#'
#' @section Cross-brand counts: These are an \emph{approximation} of ActiGraph
#'   counts, not native ActiGraph output, and GENEActiv-to-count conversion is
#'   \emph{not} empirically validated (only theoretically motivated by the shared
#'   filter). GENEActiv records at about 85.7 Hz, which is resampled to 30 Hz before
#'   the count filter. Appropriate for the relative and normalized circadian
#'   metrics here (IS, IV, RA, L5, M10, SRI); do \emph{not} apply ActiGraph
#'   cut-points or compare absolute counts across brands.
#'
#' @param path Path to a \code{.bin} file.
#' @param epoch Epoch length in seconds (default 60).
#' @param lfe Use the low-frequency extension filter (default \code{FALSE}).
#' @param tz Time zone for the timestamps (default \code{"UTC"}).
#'
#' @return Data frame with \code{time}, \code{axis1}, \code{axis2}, \code{axis3}
#'   and \code{vm}, one row per epoch.
#'
#' @references
#' \insertRef{neishabouri2022}{actiRhythm}
#'
#' \insertRef{brond2017}{actiRhythm}
#'
#' @seealso \code{\link{axivity.counts}}, \code{\link{gt3x.counts}},
#'   \code{\link{read.raw}}
#'
#' @export
geneactiv.counts <- function(path, epoch = 60, lfe = FALSE, tz = "UTC") {
  if (!requireNamespace("GGIRread", quietly = TRUE))
    stop("geneactiv.counts() requires the 'GGIRread' package: install.packages('GGIRread')")
  if (!requireNamespace("agcounts", quietly = TRUE))
    stop("geneactiv.counts() requires the 'agcounts' package: install.packages('agcounts')")
  d  <- GGIRread::readGENEActiv(path, start = 1, end = 1e6, desiredtz = tz)
  dd <- if (!is.null(d$data.out)) d$data.out else d$data
  fs <- d$header$SampleRate
  start <- as.POSIXct(dd$time[1], origin = "1970-01-01", tz = tz)
  raw <- .build_raw(data.frame(X = dd$x, Y = dd$y, Z = dd$z), fs, start, tz)
  .raw_to_counts(raw, epoch, lfe, tz)
}

#' Activity Counts from a Raw Accelerometer File (Any Supported Brand)
#'
#' A single entry point that dispatches to the brand-specific reader by file
#' extension (or an explicit \code{device}) and returns ActiGraph-equivalent
#' counts: ActiGraph \code{.gt3x}, Axivity \code{.cwa}, or GENEActiv \code{.bin}.
#'
#' @param path Path to a raw file.
#' @param device One of \code{"auto"} (infer from the extension), \code{"gt3x"},
#'   \code{"axivity"}, \code{"geneactiv"}.
#' @param epoch Epoch length in seconds (default 60).
#' @param lfe Use the low-frequency extension filter (default \code{FALSE}).
#' @param tz Time zone for the timestamps (default \code{"UTC"}).
#'
#' @return Data frame with \code{time}, \code{axis1}, \code{axis2}, \code{axis3}
#'   and \code{vm}, one row per epoch. Counts from non-ActiGraph devices are an
#'   approximation; see \code{\link{axivity.counts}} / \code{\link{geneactiv.counts}}.
#'
#' @seealso \code{\link{gt3x.counts}}, \code{\link{axivity.counts}},
#'   \code{\link{geneactiv.counts}}
#'
#' @export
read.raw <- function(path, device = c("auto", "gt3x", "axivity", "geneactiv"),
                     epoch = 60, lfe = FALSE, tz = "UTC") {
  device <- match.arg(device)
  if (device == "auto") device <- .detect_device(path)
  switch(device,
    gt3x      = gt3x.counts(path, epoch, lfe, tz),
    axivity   = axivity.counts(path, epoch, lfe, tz),
    geneactiv = geneactiv.counts(path, epoch, lfe, tz))
}
