#' Export a Full Circadian Analysis to an Excel Workbook
#'
#' Runs the complete actiRhythm circadian analysis on one activity time series and
#' writes a multi-sheet \code{.xlsx} workbook: a one-row \strong{Summary} of every
#' metric (nonparametric IS/IV/RA/L5/M10, cosinor, the rhythmicity F-test, the
#' Lomb-Scargle period with its bootstrap confidence interval, the chi-square
#' periodogram, DFA, multifractal DFA, multiscale entropy, rest-activity state
#' transitions, and, when sleep is supplied, the Sleep Regularity Index,
#' social jet lag and LIDS), plus detail sheets for the hourly profile, both
#' periodograms, the fluctuation and multifractal spectra, the transition curves,
#' the LIDS fits, and a data dictionary.
#'
#' @param activity Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param file Output \code{.xlsx} path; if \code{NULL} the workbook object is
#'   returned without writing.
#' @param sleep_state Optional per-epoch sleep/wake vector (turns on the SRI).
#' @param sleep_periods Optional sleep-period data frame with
#'   \code{in_bed_time}/\code{out_bed_time} (turns on social jet lag and LIDS).
#' @param wear_time Optional logical wear-time mask.
#' @param epoch_length Epoch length in seconds (default 60).
#' @param include_period_ci Whether to bootstrap the period CI (default TRUE).
#' @param n_boot Bootstrap replicates for the period CI (default 200).
#'
#' @return The \code{openxlsx} workbook object, invisibly.
#'
#' @examples
#' \donttest{
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' act <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24) + rnorm(length(ts), 0, 20))
#' circadian.workbook(act, ts, file = tempfile(fileext = ".xlsx"), include_period_ci = FALSE)
#' }
#'
#' @export
circadian.workbook <- function(activity, timestamps, file = NULL,
                               sleep_state = NULL, sleep_periods = NULL,
                               wear_time = NULL, epoch_length = 60,
                               include_period_ci = TRUE, n_boot = 200) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("circadian.workbook() requires the 'openxlsx' package: install.packages('openxlsx')")
  }
  a <- .circadian.analyze(activity, timestamps, sleep_state, sleep_periods,
                          wear_time, epoch_length, include_period_ci, n_boot)
  summary <- a$summary
  rar <- a$results$rar; per <- a$results$per; chi <- a$results$chi
  dfa <- a$results$dfa; mf  <- a$results$mf;  mse <- a$results$mse
  st  <- a$results$st;  lid <- a$results$lid

  wb  <- openxlsx::createWorkbook()
  add <- function(name, df) {
    openxlsx::addWorksheet(wb, name)
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) df <- data.frame(Note = "No data available")
    openxlsx::writeData(wb, name, df)
  }

  add("Summary", summary)
  add("Hourly Profile", rar$hourly_profile)
  if (!is.null(per) && length(per$scanned)) {
    add("Lomb-Scargle Periodogram", data.frame(period_h = per$scanned, power = per$power))
  }
  if (!is.null(chi) && length(chi$scanned)) {
    add("Chi-square Periodogram", data.frame(period_h = chi$scanned, Qp = chi$Qp, critical = chi$critical))
  }
  if (!is.null(dfa) && length(dfa$scales)) {
    add("Fractal DFA", data.frame(scale = dfa$scales, fluctuation = dfa$fluctuations,
                                  log10_scale = log10(dfa$scales),
                                  log10_fluctuation = log10(dfa$fluctuations)))
  }
  if (!is.null(mf) && length(mf$q_values)) {
    add("MF-DFA", data.frame(q = mf$q_values, h_q = mf$h_q, tau_q = mf$tau_q,
                             alpha = mf$alpha, f_alpha = mf$f_alpha))
  }
  if (!is.null(mse) && length(mse$scales)) {
    add("Multiscale Entropy", data.frame(scale = mse$scales, sample_entropy = mse$mse))
  }
  if (!is.null(st)) {
    add("Rest-Active Transitions", st$rest_curve)
    add("Active-Rest Transitions", st$act_curve)
  }
  if (!is.null(lid)) add("LIDS", lid$periods)
  add("Data Dictionary", .wb_dictionary())

  if (!is.null(file)) openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(wb)
}


# Run the full circadian analysis on one series; returns the one-row Summary plus
# the individual result objects. Shared by circadian.workbook() and circadian.batch().
.circadian.analyze <- function(activity, timestamps, sleep_state = NULL, sleep_periods = NULL,
                               wear_time = NULL, epoch_length = 60,
                               include_period_ci = TRUE, n_boot = 200) {
  saf <- function(expr) tryCatch(expr, error = function(e) NULL)
  gv  <- function(x, d = NA_real_) if (is.null(x) || length(x) == 0) d else x[[1]]

  rar <- saf(circadian.rhythm(activity, timestamps, sleep_state = sleep_state,
                              wear_time = wear_time, epoch_length = epoch_length))
  cos <- saf(cosinor.analysis(activity, timestamps, wear_time = wear_time))
  rhy <- saf(rhythmicity.test(activity, timestamps, cosinor_result = cos))
  per <- saf(circadian.period(activity, timestamps))
  pci <- if (include_period_ci) saf(period.ci(activity, timestamps, n_boot = n_boot)) else NULL
  dfa <- saf(fractal.dfa(activity))
  mf  <- saf(mfdfa(activity))
  mse <- saf(multiscale.entropy(activity))
  chi <- saf(chi.sq.periodogram(activity, timestamps, epoch_length = epoch_length))
  st  <- saf(state.transitions(activity))
  sri <- if (!is.null(sleep_state))   saf(sleep.regularity.index(sleep_state, timestamps, epoch_length = epoch_length)) else NULL
  sjl <- if (!is.null(sleep_periods)) saf(social.jet.lag(sleep_periods)) else NULL
  lid <- if (!is.null(sleep_periods)) saf(lids(activity, timestamps, sleep_periods, epoch_length = epoch_length)) else NULL

  summary <- data.frame(
    n_epochs = length(activity), epoch_length_s = epoch_length,
    L5 = gv(rar$L5), M10 = gv(rar$M10), RA = gv(rar$RA),
    IS = gv(rar$IS), IV = gv(rar$IV), phi = gv(rar$phi),
    cosinor_mesor = gv(cos$mesor), cosinor_amplitude = gv(cos$amplitude),
    cosinor_acrophase = gv(cos$acrophase), cosinor_r_squared = gv(cos$r_squared),
    rhythm_F = gv(rhy$F), rhythm_df2 = gv(rhy$df2), rhythm_p_value = gv(rhy$p_value),
    percent_rhythm = gv(rhy$percent_rhythm), rhythmic = gv(rhy$rhythmic, NA),
    period_tau = gv(per$tau), period_peak_power = gv(per$peak_power), period_p_value = gv(per$p_value),
    period_ci_lower = gv(pci$ci_lower), period_ci_upper = gv(pci$ci_upper), period_ci_se = gv(pci$se),
    chisq_period = gv(chi$period), chisq_Qp_peak = gv(chi$Qp_peak), chisq_p_value = gv(chi$p_value),
    dfa_alpha = gv(dfa$alpha), dfa_alpha1 = gv(dfa$alpha1), dfa_alpha2 = gv(dfa$alpha2),
    mfdfa_h2 = gv(mf$alpha_dfa), mfdfa_width = gv(mf$width),
    mse_area = gv(mse$area), mse_slope = gv(mse$slope),
    kRA = gv(st$kRA), kAR = gv(st$kAR), pRA = gv(st$pRA), pAR = gv(st$pAR),
    SRI = gv(sri),
    social_jet_lag_hours = gv(sjl$social_jet_lag_hours), MSW = gv(sjl$MSW), MSF = gv(sjl$MSF),
    lids_period_min = gv(lid$mean_period_min), lids_MRI = gv(lid$mean_MRI),
    stringsAsFactors = FALSE
  )
  list(summary = summary,
       results = list(rar = rar, cos = cos, rhy = rhy, per = per, pci = pci,
                      dfa = dfa, mf = mf, mse = mse, chi = chi, st = st,
                      sri = sri, sjl = sjl, lid = lid))
}


#' Batch Circadian Analysis of Multiple AGD Files
#'
#' Reads a set of ActiGraph \code{.agd} files, runs the full circadian analysis on
#' each, and returns one combined summary row per file, the package equivalent of
#' uploading a batch of files. Files that fail to read or analyse are reported in an
#' \code{error} column rather than aborting the run. Optionally writes a single
#' Excel workbook whose Summary sheet has one row per file.
#'
#' @param files Character vector of \code{.agd} paths, or a single directory (all
#'   its \code{.agd} files are used).
#' @param file Optional output \code{.xlsx} path for a combined-summary workbook.
#' @param metric Activity metric: \code{"axis1"} (default) or \code{"vm"} (vector
#'   magnitude).
#' @param include_period_ci Bootstrap the period CI per file (default \code{FALSE}
#'   for batch speed).
#' @param n_boot Bootstrap replicates when \code{include_period_ci} is \code{TRUE}.
#' @param epoch_length Epoch length in seconds; \code{NULL} (default) infers it per
#'   file from the timestamps.
#' @param verbose Print per-file progress (default \code{TRUE}).
#'
#' @return A data frame with one row per file (a \code{file} column, an
#'   \code{error} column, and every summary metric). Returned invisibly when a
#'   workbook \code{file} is written.
#'
#' @seealso \code{\link{circadian.workbook}}
#'
#' @examples
#' \donttest{
#' # Pass a folder to analyse a whole batch; a single file is used here.
#' batch <- circadian.batch(example_agd(), verbose = FALSE)
#' batch[, c("file", "IS", "IV", "RA", "rhythm_p_value")]
#' }
#'
#' @export
circadian.batch <- function(files, file = NULL, metric = c("axis1", "vm"),
                            include_period_ci = FALSE, n_boot = 200,
                            epoch_length = NULL, verbose = TRUE) {
  metric <- match.arg(metric)
  if (length(files) == 1L && dir.exists(files)) {
    files <- list.files(files, pattern = "[.]agd$", full.names = TRUE)
  }
  if (!length(files)) stop("no .agd files supplied")

  rows <- lapply(seq_along(files), function(i) {
    f <- files[i]
    if (verbose) message(sprintf("[%d/%d] %s", i, length(files), basename(f)))
    tryCatch({
      agd <- agd.counts(read.agd(f, verbose = FALSE))
      act <- if (metric == "vm" && all(c("axis2", "axis3") %in% names(agd))) {
        sqrt(agd$axis1^2 + agd$axis2^2 + agd$axis3^2)
      } else agd$axis1
      el <- if (!is.null(epoch_length)) epoch_length else {
        d <- as.numeric(difftime(agd$timestamp[2], agd$timestamp[1], units = "secs"))
        if (is.finite(d) && d > 0) d else 60
      }
      a <- .circadian.analyze(act, agd$timestamp, epoch_length = el,
                              include_period_ci = include_period_ci, n_boot = n_boot)
      cbind(file = basename(f), error = NA_character_, a$summary, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(file = basename(f), error = conditionMessage(e), stringsAsFactors = FALSE)
    })
  })
  out <- .bind_rows_fill(rows)

  if (!is.null(file)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("writing a workbook requires the 'openxlsx' package")
    }
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary");         openxlsx::writeData(wb, "Summary", out)
    openxlsx::addWorksheet(wb, "Data Dictionary"); openxlsx::writeData(wb, "Data Dictionary", .wb_dictionary())
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    return(invisible(out))
  }
  out
}

# rbind a list of data frames with differing columns (column union, NA-filled).
.bind_rows_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(data.frame())
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] <- NA
    d[, cols, drop = FALSE]
  }))
}


# Short column -> definition table for the workbook's Data Dictionary sheet.
.wb_dictionary <- function() {
  rows <- list(
    c("IS / IV / RA", "Interdaily stability, intradaily variability, relative amplitude"),
    c("L5 / M10", "Least-active 5 h and most-active 10 h mean activity"),
    c("cosinor_*", "Single-component cosinor MESOR, amplitude, acrophase, R-squared"),
    c("rhythm_F / rhythm_p_value / percent_rhythm", "Halberg zero-amplitude F-test and percent rhythm"),
    c("period_tau / period_ci_*", "Lomb-Scargle endogenous period and its bootstrap confidence interval"),
    c("chisq_*", "Chi-square (Sokolove-Bushell) periodogram period, peak Q, and p-value"),
    c("dfa_alpha* / mfdfa_*", "Detrended-fluctuation scaling exponent(s) and multifractal spectrum width"),
    c("mse_*", "Multiscale entropy area and slope"),
    c("kRA / kAR / pRA / pAR", "Rest-activity state transition rates"),
    c("SRI", "Sleep Regularity Index (requires sleep_state)"),
    c("social_jet_lag_* / MSW / MSF", "Social jet lag and mid-sleep on work/free days (requires sleep_periods)"),
    c("lids_*", "Locomotor Inactivity During Sleep period and Munich Rhythmicity Index (requires sleep_periods)")
  )
  data.frame(Column = vapply(rows, `[`, "", 1), Definition = vapply(rows, `[`, "", 2),
             stringsAsFactors = FALSE)
}


#' Save a actiRhythm Plot to a File
#'
#' Convenience wrapper around \code{ggplot2::ggsave} for the \code{plot_*}
#' circadian figures; the output format is taken from the file extension
#' (\code{.png}, \code{.pdf}, \code{.svg}).
#'
#' @param plot A ggplot object (e.g. from \code{\link{plot_actogram}}).
#' @param file Output path; its extension sets the format.
#' @param width,height Size in inches (defaults 8 x 5).
#' @param dpi Resolution for raster formats (default 300).
#' @return The \code{file} path, invisibly.
#' @export
save.circadian.plot <- function(plot, file, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi)
  invisible(file)
}
