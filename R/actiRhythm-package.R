#' actiRhythm: Circadian Rest-Activity Rhythm Analysis
#'
#' Quantifies the circadian rest-activity rhythm of a single subject from activity
#' counts: nonparametric measures, cosinor models with a rhythmicity test, period
#' estimation, fractal and nonlinear measures, the Sleep Regularity Index, and
#' phase metrics. Analyses run on an activity-count vector and its timestamps, and
#' a built-in reader loads 'ActiGraph' '.agd' files directly.
#'
#' @importFrom stats acf aggregate approx lm.fit na.pass pf qt sd
#' @importFrom ggplot2 .data
#' @importFrom grDevices colorRampPalette
#' @importFrom scales squish
#' @importFrom utils globalVariables
#' @importFrom Rdpack reprompt
#' @keywords internal
"_PACKAGE"

# ggplot2 aes() variables referenced inside plot methods (non-standard evaluation).
globalVariables(c("hour", "mean_counts", "sd_counts", "RA"))
