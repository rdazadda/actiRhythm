#' Rest-Activity State Transition Rates (kRA, kAR)
#'
#' Computes the rest-to-activity and activity-to-rest transition rates from a
#' binarized activity series, following the survival-curve method used by
#' pyActigraphy (Lim et al. 2011). Thresholds the series into rest/active epochs,
#' builds the per-lag transition probability (hazard) of each bout type from the
#' bout-length survival curve, and takes a single rate over the LOWESS
#' "sustained" plateau of that curve.
#'
#' @param counts Numeric activity vector; \code{NA} are dropped.
#' @param threshold Activity level at or above which an epoch is "active"
#'   (default 1, i.e. any non-zero count).
#' @param frac LOWESS smoother span for the sustained-region search (default 0.3).
#' @param iter LOWESS robustifying iterations (default 0).
#'
#' @return An object of class \code{actiRhythm_transitions}: a list with
#'   \code{kRA}/\code{kAR} (sustained rest-to-active / active-to-rest rates over
#'   the LOWESS plateau), \code{pRA}/\code{pAR} (overall rest-to-active /
#'   active-to-rest rate, the reciprocal mean bout length), bout counts, and the
#'   two transition curves. The plateau search follows the pyActigraphy
#'   implementation of Lim et al. (2011).
#'
#' @references
#' \insertRef{lim2011}{actiRhythm}
#'
#' @examples
#' set.seed(1)
#' counts <- as.integer(stats::runif(5000) < 0.1) * 100
#' state.transitions(counts)
#'
#' @export
state.transitions <- function(counts, threshold = 1, frac = 0.3, iter = 0) {
  active <- as.integer(counts >= threshold)
  active <- active[!is.na(active)]
  if (length(active) < 4L) {
    return(structure(list(kRA = NA_real_, kAR = NA_real_, pRA = NA_real_, pAR = NA_real_,
      threshold = threshold, n_rest_bouts = 0L, n_act_bouts = 0L,
      rest_curve = NULL, act_curve = NULL, insufficient = TRUE),
      class = "actiRhythm_transitions"))
  }

  r <- rle(active)
  rest_bouts <- r$lengths[r$values == 0L]
  act_bouts  <- r$lengths[r$values == 1L]

  rest_curve <- .transition_curve(rest_bouts)
  act_curve  <- .transition_curve(act_bouts)

  structure(list(
    kRA = .sustain_rate(rest_curve, frac, iter),
    kAR = .sustain_rate(act_curve,  frac, iter),
    pRA = if (length(rest_bouts)) length(rest_bouts) / sum(rest_bouts) else NA_real_,
    pAR = if (length(act_bouts))  length(act_bouts)  / sum(act_bouts)  else NA_real_,
    threshold    = threshold,
    n_rest_bouts = length(rest_bouts),
    n_act_bouts  = length(act_bouts),
    rest_curve   = rest_curve,
    act_curve    = act_curve,
    insufficient = FALSE
  ), class = "actiRhythm_transitions")
}


# Bout-length survival curve -> per-lag transition probability (hazard), with the
# gap correction and sqrt-count weights used by pyActigraphy.
.transition_curve <- function(bouts) {
  bouts <- bouts[is.finite(bouts) & bouts > 0]
  if (length(bouts) < 3L) return(NULL)
  tab  <- table(bouts)
  t    <- as.integer(names(tab))
  n_at <- as.numeric(tab)
  if (length(t) < 2L) return(NULL)

  Nt   <- rev(cumsum(rev(n_at)))        # number of bouts of length >= t[i]
  dNt  <- Nt - c(Nt[-1], 0)             # number terminating exactly at t[i]
  prob <- (dNt / Nt)[-length(t)] / diff(t)   # gap-corrected per-epoch hazard
  w    <- sqrt(Nt[-length(Nt)] + Nt[-1])
  data.frame(lag = t[-length(t)], prob = prob, weight = w)
}

# Single rate over the LOWESS-defined sustained plateau of a transition curve.
.sustain_rate <- function(curve, frac = 0.3, iter = 0) {
  if (is.null(curve) || nrow(curve) == 0L) return(NA_real_)
  if (nrow(curve) < 3L) return(stats::weighted.mean(curve$prob, curve$weight))

  sm <- stats::lowess(curve$lag, curve$prob, f = frac, iter = iter)$y
  s  <- stats::sd(curve$prob)
  if (!is.finite(s) || s == 0) return(stats::weighted.mean(curve$prob, curve$weight))

  below <- abs(curve$prob - sm) < s
  rr <- rle(below)
  on <- which(rr$values)
  if (!length(on)) return(stats::weighted.mean(curve$prob, curve$weight))

  best   <- on[which.max(rr$lengths[on])]
  ends   <- cumsum(rr$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  span   <- starts[best]:ends[best]
  stats::weighted.mean(curve$prob[span], curve$weight[span])
}


#' @export
print.actiRhythm_transitions <- function(x, ...) {
  cat("Rest-Activity State Transitions\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Threshold:  >= %g counts = active\n", x$threshold))
  cat(sprintf("  kRA (rest->active): %s   (%d rest bouts)\n",
              formatC(x$kRA, format = "f", digits = 4), x$n_rest_bouts))
  cat(sprintf("  kAR (active->rest): %s   (%d active bouts)\n",
              formatC(x$kAR, format = "f", digits = 4), x$n_act_bouts))
  cat(sprintf("  pRA / pAR:          %s / %s\n",
              formatC(x$pRA, format = "f", digits = 4),
              formatC(x$pAR, format = "f", digits = 4)))
  invisible(x)
}


#' Rest-Active Transition Probabilities
#'
#' Maximum-likelihood and Bayesian estimates of the rest-to-active and
#' active-to-rest transition probabilities (Danilevicz et al. 2024): the
#' transitions out of a state divided by its epochs at risk.
#'
#' @param counts Numeric activity vector; \code{NA} are dropped.
#' @param threshold Counts at or above which an epoch is active (default 1).
#' @param eps Bayesian Beta pseudo-count added to the transition count and the
#'   epochs at risk (default 0.5).
#'
#' @return A list with the maximum-likelihood and Bayesian \code{tp_ra} (rest to
#'   active) and \code{tp_ar} (active to rest), the active bout count, and the
#'   mean active bout length.
#'
#' @references
#' \insertRef{danilevicz2024}{actiRhythm}
#'
#' @seealso \code{\link{state.transitions}}, \code{\link{activity.balance.index}}
#'
#' @examples
#' counts <- c(rep(0, 50), rep(100, 20), rep(0, 40), rep(80, 30), rep(0, 60))
#' transition.probability(counts)
#'
#' @export
transition.probability <- function(counts, threshold = 1, eps = 0.5) {
  active <- suppressWarnings(as.numeric(counts)) >= threshold
  active <- active[!is.na(active)]
  r <- rle(active); v <- r$values; nv <- length(v)
  if (!nv) return(list(tp_ar_mle = NA_real_, tp_ra_mle = NA_real_,
    tp_ar_bayes = NA_real_, tp_ra_bayes = NA_real_,
    n_active_bouts = 0L, mean_active_bout = NA_real_))
  T_active <- sum(active); T_rest <- sum(!active)
  n_ar <- if (nv > 1L) sum(v[-nv] & !v[-1]) else 0L
  n_ra <- if (nv > 1L) sum(!v[-nv] & v[-1]) else 0L
  rest_risk <- T_rest   - as.integer(!v[nv])
  act_risk  <- T_active - as.integer(v[nv])
  mle   <- function(n, risk) if (risk > 0) n / risk else NA_real_
  bayes <- function(n, risk) if (risk + eps > 0) (n + eps) / (risk + eps) else NA_real_
  list(tp_ar_mle = mle(n_ar, act_risk), tp_ra_mle = mle(n_ra, rest_risk),
       tp_ar_bayes = bayes(n_ar, act_risk), tp_ra_bayes = bayes(n_ra, rest_risk),
       n_active_bouts = sum(v), mean_active_bout = if (any(v)) mean(r$lengths[v]) else NA_real_)
}
