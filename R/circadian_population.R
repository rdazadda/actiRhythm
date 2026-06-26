#' Population-Mean Cosinor (Bingham)
#'
#' Pools single-subject cosinor fits into a group-mean rhythm with confidence
#' intervals, following Bingham et al. (1982). Fits each subject's averaged
#' 24-hour profile with the same weighted-least-squares engine as
#' \code{\link{cosinor.analysis}}, averages the linearized cos/sin coefficients
#' across subjects, and returns the group MESOR, amplitude, and acrophase with
#' Bingham confidence intervals.
#'
#' @param activity Numeric activity vector with all subjects stacked together.
#' @param timestamps POSIXct timestamps, one per value.
#' @param subject Subject identifier, one per value.
#' @param group Optional group identifier, one per value; when supplied a
#'   population cosinor is returned for each group.
#' @param period Rhythm period in hours (default 24).
#' @param level Confidence level (default 0.95).
#' @param min_valid_hours Minimum profile hours for a subject to be included
#'   (default 12).
#'
#' @return An object of class \code{actiRhythm_population_cosinor} (group MESOR,
#'   amplitude, acrophase with Bingham CIs and a \code{conf_interval_valid}
#'   flag), or a named list of them (class
#'   \code{actiRhythm_population_cosinor_list}) when \code{group} is supplied.
#'
#' @references
#' \insertRef{bingham1982}{actiRhythm}
#'
#' @examples
#' set.seed(1)
#' hrs <- 0:23
#' act <- ts <- subj <- NULL
#' for (i in 1:6) {
#'   y <- 100 + 40 * cos(2 * pi * (hrs - (8 + i / 3)) / 24) + rnorm(24, 0, 4)
#'   act <- c(act, y); subj <- c(subj, rep(paste0("S", i), 24))
#'   ts <- c(ts, as.POSIXct("2024-01-01", tz = "UTC") + hrs * 3600)
#' }
#' population.cosinor(act, as.POSIXct(ts, tz = "UTC", origin = "1970-01-01"), subj)
#'
#' @export
population.cosinor <- function(activity, timestamps, subject, group = NULL,
                               period = 24, level = 0.95, min_valid_hours = 12) {
  n <- length(activity)
  if (length(timestamps) != n || length(subject) != n) {
    stop("'activity', 'timestamps' and 'subject' must have the same length")
  }
  if (!is.null(group) && length(group) != n) {
    stop("'group' must have the same length as 'activity'")
  }

  subject <- as.character(subject)
  grp_of  <- if (is.null(group)) NULL else as.character(group)

  coefs <- .population.coefs(activity, timestamps, subject, grp_of, period, min_valid_hours)
  if (is.null(coefs) || nrow(coefs) < 2) {
    stop("need at least 2 subjects with >= ", min_valid_hours, " valid profile hours")
  }
  n_dropped <- length(unique(subject)) - nrow(coefs)

  run_group <- function(cf) {
    res <- .bingham.population(cf$mesor, cf$beta1, cf$beta2, period, level)
    res$n_subjects <- nrow(cf)
    res$subjects   <- cf
    res$period     <- period
    res$level      <- level
    class(res) <- "actiRhythm_population_cosinor"
    res
  }

  if (is.null(group)) {
    out <- run_group(coefs)
    out$n_dropped <- n_dropped
    return(out)
  }
  out <- lapply(split(coefs, coefs$group), run_group)
  class(out) <- "actiRhythm_population_cosinor_list"
  out
}


# Bingham population-mean cosinor on linearized per-subject coefficients.
# MESOR/amplitude/acrophase are in actiRhythm's native (beta1, beta2,
# positive-hours) convention; the acrophase CI is evaluated in cosinor2's
# (gamma = -beta2, negative-radian) convention and mapped back to hours.
.bingham.population <- function(M, b1, b2, period = 24, level = 0.95) {
  omega <- 2 * pi / period
  k <- length(M)
  Mbar <- mean(M); B <- mean(b1); G <- mean(b2)
  A <- sqrt(B^2 + G^2)
  acrophase <- (atan2(G, B) / omega) %% period

  tcrit <- stats::qt(1 - (1 - level) / 2, df = k - 1)
  sdM <- stats::sd(M)
  ci_mesor <- Mbar + c(-1, 1) * tcrit * sdM / sqrt(k)

  ci_amplitude <- c(NA_real_, NA_real_)
  ci_acrophase <- c(NA_real_, NA_real_)
  ci_valid <- FALSE

  if (k >= 3 && A > 0) {
    sdB <- stats::sd(b1); sdG <- stats::sd(b2); covBG <- stats::cov(b1, b2)
    denom <- A^2 * k
    c22 <- (sdB^2 * B^2 + 2 * covBG * B * G + sdG^2 * G^2) / denom
    c33 <- (sdB^2 * G^2 - 2 * covBG * B * G + sdG^2 * B^2) / denom
    c23 <- -((-(sdB^2 - sdG^2) * B * G + covBG * (B^2 - G^2)) / denom)  # gamma = -beta2

    ci_amplitude <- c(max(0, A - tcrit * sqrt(c22)), A + tcrit * sqrt(c22))

    den      <- A^2 - c22 * tcrit^2
    radicand <- A^2 - ((c22 * c33 - c23^2) * tcrit^2) / c33
    if (is.finite(c33) && c33 > 0 && is.finite(radicand) && radicand >= 0 && den > 0) {
      root   <- tcrit * sqrt(c33) * sqrt(radicand)
      bounds <- atan2(-G, B) + atan(c(c23 * tcrit^2 + root, c23 * tcrit^2 - root) / den)
      ci_acrophase <- sort((-bounds / omega) %% period)
      ci_valid <- TRUE
    }
  }

  list(
    mesor = Mbar, amplitude = A, acrophase = acrophase,
    ci_mesor = ci_mesor, ci_amplitude = ci_amplitude, ci_acrophase = ci_acrophase,
    conf_interval_valid = ci_valid,
    rhythm_detected = ci_valid && ci_amplitude[1] > 0,
    k = k
  )
}


#' @export
print.actiRhythm_population_cosinor <- function(x, ...) {
  cat("Population-Mean Cosinor (Bingham)\n\n")
  cat(sprintf("  Subjects:   %d\n", x$n_subjects))
  cat(sprintf("  Period:     %g h\n", x$period))
  cat(sprintf("  MESOR:      %.2f  [%.2f, %.2f]\n",
              x$mesor, x$ci_mesor[1], x$ci_mesor[2]))
  if (isTRUE(x$conf_interval_valid)) {
    cat(sprintf("  Amplitude:  %.2f  [%.2f, %.2f]\n",
                x$amplitude, x$ci_amplitude[1], x$ci_amplitude[2]))
    cat(sprintf("  Acrophase:  %.2f h  [%.2f, %.2f]\n",
                x$acrophase, x$ci_acrophase[1], x$ci_acrophase[2]))
  } else {
    cat(sprintf("  Amplitude:  %.2f  (CI undefined; need >= 3 subjects / clearer rhythm)\n",
                x$amplitude))
    cat(sprintf("  Acrophase:  %.2f h\n", x$acrophase))
  }
  invisible(x)
}


#' @export
print.actiRhythm_population_cosinor_list <- function(x, ...) {
  for (g in names(x)) {
    cat("== Group:", g, "==\n")
    print(x[[g]])
    cat("\n")
  }
  invisible(x)
}


# Per-subject linearized cosinor coefficients (same engine as cosinor.analysis).
.population.coefs <- function(activity, timestamps, subject, grp_of, period, min_valid_hours) {
  subject <- as.character(subject)
  rows <- lapply(unique(subject), function(s) {
    sel  <- subject == s
    prof <- .ext.hourly.profile(activity[sel], timestamps[sel])
    if (is.null(prof) || length(prof$t) < min_valid_hours) return(NULL)
    fit <- .ext.ordinary.cosinor(prof$t, prof$y, prof$n, period = period)
    if (is.null(fit)) return(NULL)
    data.frame(subject = s,
               group = if (is.null(grp_of)) "all" else grp_of[which(sel)[1]],
               mesor = fit$mesor, beta1 = fit$beta1, beta2 = fit$beta2,
               amplitude = fit$amplitude, acrophase = fit$acrophase,
               n_hours = length(prof$t), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


#' Compare Cosinor Rhythms Between Two Groups
#'
#' Tests whether the rest-activity rhythm differs between two groups. Fits each
#' subject's cosinor with the same engine as \code{\link{cosinor.analysis}}. An
#' overall multivariate test (Hotelling's \eqn{T^2} on the joint MESOR, cosine,
#' and sine parameters, the Bingham et al. (1982) population-cosinor comparison)
#' gives one omnibus p-value for whether the rhythms differ at all, accounting
#' for the correlation between the cosine and sine components. The per-parameter
#' two-sample (Welch) t-tests then break that result down for the MESOR,
#' amplitude, and acrophase. The acrophase test is circular-aware, unwrapping the
#' per-subject acrophases about their common circular mean.
#'
#' @param activity Numeric activity vector with all subjects stacked.
#' @param timestamps POSIXct timestamps, one per value.
#' @param subject Subject identifier, one per value.
#' @param group Group identifier (exactly two levels), one per value.
#' @param period Rhythm period in hours (default 24).
#' @param level Confidence level for the difference CIs (default 0.95).
#' @param min_valid_hours Minimum profile hours for a subject to be included.
#'
#' @return An object of class \code{actiRhythm_cosinor_compare}: a \code{joint}
#'   list (the omnibus Hotelling \eqn{T^2}, its F statistic, degrees of freedom
#'   and p-value), a \code{tests} data frame (one row per parameter, with each
#'   group's estimate, their difference, the t statistic, degrees of freedom,
#'   p-value and CI), and the per-subject coefficients.
#'
#' @references
#' \insertRef{bingham1982}{actiRhythm}
#'
#' @seealso \code{\link{population.cosinor}}
#'
#' @examples
#' set.seed(1); hrs <- 0:23
#' act <- ts <- subj <- grp <- NULL
#' for (g in c("A", "B")) for (i in 1:5) {
#'   acro <- if (g == "A") 8 else 11
#'   y <- 100 + 40 * cos(2 * pi * (hrs - acro) / 24) + rnorm(24, 0, 4)
#'   act <- c(act, y); subj <- c(subj, rep(paste0(g, i), 24)); grp <- c(grp, rep(g, 24))
#'   ts <- c(ts, as.POSIXct("2024-01-01", tz = "UTC") + hrs * 3600)
#' }
#' cosinor.compare(act, as.POSIXct(ts, tz = "UTC", origin = "1970-01-01"), subj, grp)
#'
#' @export
cosinor.compare <- function(activity, timestamps, subject, group, period = 24,
                            level = 0.95, min_valid_hours = 12) {
  n <- length(activity)
  if (length(timestamps) != n || length(subject) != n || length(group) != n) {
    stop("'activity', 'timestamps', 'subject' and 'group' must have the same length")
  }
  grp_of <- as.character(group)
  if (length(unique(grp_of)) != 2L) stop("'group' must have exactly two levels")

  coefs <- .population.coefs(activity, timestamps, as.character(subject), grp_of,
                             period, min_valid_hours)
  if (is.null(coefs)) stop("no subject had >= ", min_valid_hours, " valid profile hours")
  gl <- sort(unique(coefs$group))
  g1 <- coefs[coefs$group == gl[1], , drop = FALSE]
  g2 <- coefs[coefs$group == gl[2], , drop = FALSE]
  if (nrow(g1) < 2L || nrow(g2) < 2L) stop("each group needs >= 2 subjects with a valid fit")

  acro <- .center.hours(coefs$acrophase, period)   # circular-aware acrophase
  a1 <- acro[coefs$group == gl[1]]; a2 <- acro[coefs$group == gl[2]]

  tests <- rbind(
    .two.sample.test("mesor",     g1$mesor,     g2$mesor,     level),
    .two.sample.test("amplitude", g1$amplitude, g2$amplitude, level),
    .two.sample.test("acrophase", a1,           a2,           level)
  )
  structure(list(
    joint = .bingham.joint.test(g1, g2), tests = tests,
    groups = gl, n1 = nrow(g1), n2 = nrow(g2),
    period = period, level = level, subjects = coefs
  ), class = "actiRhythm_cosinor_compare")
}


# Bingham (1982) omnibus comparison: Hotelling's T^2 on the joint (MESOR, cosine,
# sine) parameter vectors of two groups with pooled within-group covariance.
.bingham.joint.test <- function(g1, g2) {
  X1 <- as.matrix(g1[, c("mesor", "beta1", "beta2")])
  X2 <- as.matrix(g2[, c("mesor", "beta1", "beta2")])
  k1 <- nrow(X1); k2 <- nrow(X2); K <- k1 + k2
  na <- list(valid = FALSE, T2 = NA_real_, statistic = NA_real_,
             df1 = 3L, df2 = K - 4L, p_value = NA_real_)
  if (K - 4L < 1L) return(na)
  Sp <- ((k1 - 1) * stats::cov(X1) + (k2 - 1) * stats::cov(X2)) / (K - 2)
  d <- colMeans(X1) - colMeans(X2)
  tryCatch({
    T2 <- (k1 * k2 / K) * as.numeric(t(d) %*% solve(Sp) %*% d)
    Fstat <- T2 * (K - 4) / (3 * (K - 2))
    list(valid = TRUE, T2 = T2, statistic = Fstat, df1 = 3L, df2 = K - 4L,
         p_value = stats::pf(Fstat, 3, K - 4L, lower.tail = FALSE))
  }, error = function(e) na)
}


.two.sample.test <- function(parameter, x1, x2, level) {
  tt <- stats::t.test(x1, x2, conf.level = level)
  data.frame(parameter = parameter,
             estimate1 = mean(x1), estimate2 = mean(x2),
             difference = mean(x1) - mean(x2),
             statistic = unname(tt$statistic), df = unname(tt$parameter),
             p_value = tt$p.value, ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
             stringsAsFactors = FALSE)
}

# Unwrap clock hours about their circular mean so a t-test is meaningful even
# when acrophases straddle the 0/period boundary.
.center.hours <- function(h, period = 24) {
  ang <- 2 * pi * h / period
  mu  <- atan2(mean(sin(ang)), mean(cos(ang))) * period / (2 * pi)
  mu + (((h - mu + period / 2) %% period) - period / 2)
}


#' @export
print.actiRhythm_cosinor_compare <- function(x, ...) {
  cat("Cosinor Comparison Between Two Groups\n\n")
  cat(sprintf("  Groups:  %s (n=%d)  vs  %s (n=%d)\n", x$groups[1], x$n1, x$groups[2], x$n2))
  cat(sprintf("  Period:  %g h\n\n", x$period))
  if (isTRUE(x$joint$valid)) {
    cat(sprintf("  Joint (Bingham/Hotelling T2):  F(%d,%d) = %.2f   p = %s\n\n",
                x$joint$df1, x$joint$df2, x$joint$statistic,
                format.pval(x$joint$p_value, digits = 3)))
  }
  t <- x$tests
  for (i in seq_len(nrow(t))) {
    cat(sprintf("  %-9s diff = %+.2f   t(%.1f) = %.2f   p = %s\n",
                t$parameter[i], t$difference[i], t$df[i], t$statistic[i],
                format.pval(t$p_value[i], digits = 3)))
  }
  invisible(x)
}
