#' Multicomponent Cosinor with Model Selection
#'
#' Fits the multi-component (harmonic) cosinor of Cornelissen (2014) with one to
#' several harmonics of the fundamental period and picks the number of harmonics by
#' an information criterion (AIC or BIC, a package choice), so it captures a bimodal
#' or asymmetric daily shape without your choosing the order by hand. The single
#' cosinor is the one-harmonic special case.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param period Fundamental period in hours (default 24).
#' @param max_harmonics Largest number of harmonics to consider (default 3).
#' @param criterion \code{"AIC"} (default) or \code{"BIC"} for model selection.
#'
#' @return An object of class \code{actiRhythm_multicosinor}: the selected number
#'   of harmonics, the per-harmonic amplitude and acrophase, the MESOR, R-squared,
#'   and the full model-comparison table. Never errors.
#'
#' @references
#' \insertRef{cornelissen2014}{actiRhythm}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' th <- as.numeric(difftime(ts, ts[1], units = "hours"))
#' y  <- 100 + 40 * cos(2 * pi * th / 24) + 20 * cos(2 * pi * 2 * th / 24)
#' cosinor.multicomponent(y, ts)
#'
#' @export
cosinor.multicomponent <- function(counts, timestamps, period = 24,
                                   max_harmonics = 3, criterion = c("AIC", "BIC")) {
  criterion <- match.arg(criterion)
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  y <- suppressWarnings(as.numeric(counts))
  th <- as.numeric(difftime(timestamps, min(timestamps, na.rm = TRUE), units = "hours"))
  ok <- is.finite(y) & is.finite(th); y <- y[ok]; th <- th[ok]
  na_out <- structure(list(n_harmonics = NA_integer_, mesor = NA_real_,
    harmonics = data.frame(), r_squared = NA_real_, comparison = data.frame(),
    criterion = criterion, insufficient = TRUE),
    class = c("actiRhythm_multicosinor", "list"))
  if (length(y) < 4L || stats::sd(y) == 0) return(na_out)

  # Average onto a phase-folded profile and fit with sqrt-count weights, matching
  # cosinor.analysis so R-squared and AIC/BIC are comparable across the package.
  nb  <- max(round(period), 4L)
  bin <- pmin(floor((th %% period) / period * nb) + 1L, nb)
  py  <- tapply(y, bin, mean); pn <- tapply(y, bin, length)
  pt  <- (as.integer(names(py)) - 0.5) * period / nb
  py  <- as.numeric(py); pn <- as.numeric(pn); w <- sqrt(pn)
  np  <- length(py)
  if (np < 3L) return(na_out)
  max_harmonics <- max(1L, min(max_harmonics, (np - 1L) %/% 2L))
  tss <- sum(pn * (py - stats::weighted.mean(py, pn))^2)
  N   <- length(y)
  fit_k <- function(K) {
    Xp <- cbind(1, do.call(cbind, lapply(seq_len(K), function(k)
      cbind(cos(2 * pi * k * pt / period), sin(2 * pi * k * pt / period)))))
    b   <- stats::lm.fit(Xp * w, py * w)$coefficients     # WLS coefficients on the profile
    prss <- sum(pn * (py - Xp %*% b)^2)                    # profile RSS -> R-squared
    Xr  <- cbind(1, do.call(cbind, lapply(seq_len(K), function(k)
      cbind(cos(2 * pi * k * th / period), sin(2 * pi * k * th / period)))))
    rrss <- sum((y - Xr %*% b)^2)                          # raw-epoch RSS -> AIC/BIC selection
    p   <- length(b) + 1
    list(K = K, b = b,
         AIC = N * log(rrss / N) + 2 * p, BIC = N * log(rrss / N) + log(N) * p,
         r2 = 1 - prss / tss)
  }
  fits <- lapply(seq_len(max_harmonics), fit_k)
  comp <- data.frame(harmonics = vapply(fits, `[[`, integer(1), "K"),
                     AIC = vapply(fits, `[[`, numeric(1), "AIC"),
                     BIC = vapply(fits, `[[`, numeric(1), "BIC"),
                     r_squared = vapply(fits, `[[`, numeric(1), "r2"))
  best <- fits[[which.min(comp[[criterion]])]]
  b <- best$b; K <- best$K
  harm <- lapply(seq_len(K), function(k) {
    A <- b[2 * k]; B <- b[2 * k + 1]
    amp <- sqrt(A^2 + B^2)
    acro <- ((atan2(B, A)) %% (2 * pi)) / (2 * pi) * (period / k)
    data.frame(harmonic = k, amplitude = amp, acrophase_h = acro)
  })
  structure(list(n_harmonics = K, mesor = unname(b[1]),
    harmonics = do.call(rbind, harm), r_squared = best$r2,
    comparison = comp, criterion = criterion, insufficient = FALSE),
    class = c("actiRhythm_multicosinor", "list"))
}

#' @export
print.actiRhythm_multicosinor <- function(x, ...) {
  cat("Multicomponent Cosinor\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Selected harmonics: %d (by %s)   MESOR: %.1f   R-squared: %.3f\n\n",
              x$n_harmonics, x$criterion, x$mesor, x$r_squared))
  h <- x$harmonics; h$amplitude <- round(h$amplitude, 2); h$acrophase_h <- round(h$acrophase_h, 2)
  print(h, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' Activity Onset and Offset (Relative-Difference Phase Markers)
#'
#' Finds the daily activity onset and offset by a relative-difference contrast on
#' the averaged 24-hour profile: the onset is where mean activity rises most
#' sharply (the relative difference of the window after versus before is largest)
#' and the offset is where it falls most sharply. These are non-cosinor,
#' non-changepoint phase markers, a normalized-contrast edge detector on the daily
#' profile rather than a published actigraphy algorithm.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param window_hours Half-window, in hours, compared before versus after each
#'   minute (default 6).
#'
#' @return An object of class \code{actiRhythm_aont}: \code{onset_h} and
#'   \code{offset_h} (clock hours) and the relative-difference profile. Never errors.
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
#' h  <- as.numeric(format(ts, "%H"))
#' activity.onset.offset(ifelse(h >= 23 | h < 7, 5, 300), ts)
#'
#' @export
activity.onset.offset <- function(counts, timestamps, window_hours = 6) {
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  mod <- as.POSIXlt(ts)$hour * 60 + as.POSIXlt(ts)$min
  prof <- as.numeric(tapply(x, factor(mod, levels = 0:1439), mean, na.rm = TRUE))
  na_out <- structure(list(onset_h = NA_real_, offset_h = NA_real_,
    rel_diff = rep(NA_real_, 1440), insufficient = TRUE), class = c("actiRhythm_aont", "list"))
  if (all(is.na(prof))) return(na_out)
  prof[is.na(prof)] <- mean(prof, na.rm = TRUE)
  w <- as.integer(window_hours * 60)
  cp <- c(prof, prof, prof)                   # three copies, centre minute i -> 1440 + i
  rd <- vapply(seq_len(1440), function(i) {
    ci <- 1440L + i
    before <- mean(cp[(ci - w):(ci - 1L)]); after <- mean(cp[ci:(ci + w - 1L)])
    if (before + after == 0) 0 else (after - before) / (after + before)
  }, numeric(1))
  structure(list(onset_h = (which.max(rd) - 1) / 60, offset_h = (which.min(rd) - 1) / 60,
    rel_diff = rd, insufficient = FALSE), class = c("actiRhythm_aont", "list"))
}

#' @export
print.actiRhythm_aont <- function(x, ...) {
  cat("Activity Onset / Offset\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  Activity onset:  %05.2f h\n", x$onset_h))
  cat(sprintf("  Activity offset: %05.2f h\n\n", x$offset_h))
  invisible(x)
}


#' Phase Concentration Tests
#'
#' Tests whether a set of daily phase markers (acrophases, onsets, L5/M10 times)
#' are concentrated rather than scattered around the clock. Reports the mean
#' resultant vector, the Rayleigh test of uniformity (Fisher 1993), and the
#' Hermans-Rasson test (Landler et al. 2019), which catches multimodal
#' clustering that Rayleigh misses.
#'
#' @param times_h Numeric vector of clock times (hours, 0-24), one per day.
#' @param period Period the times wrap on, in hours (default 24).
#' @param n_perm Permutations for the Hermans-Rasson p-value (default 2000, a
#'   speed tradeoff; Landler et al. 2019 use 9999 for a finer minimum p-value).
#'
#' @return An object of class \code{actiRhythm_phasetest}: mean direction, mean
#'   resultant length R, and the Rayleigh and Hermans-Rasson statistics and
#'   p-values. Never errors.
#'
#' @references
#' \insertRef{fisher1993}{actiRhythm}
#'
#' \insertRef{landler2019}{actiRhythm}
#'
#' @examples
#' set.seed(1)
#' onsets <- 23 + stats::rnorm(10, 0, 0.5)      # tightly clustered near 23:00
#' phase.concentration(onsets %% 24)
#'
#' @export
phase.concentration <- function(times_h, period = 24, n_perm = 2000) {
  th <- 2 * pi * (suppressWarnings(as.numeric(times_h)) %% period) / period
  th <- th[is.finite(th)]
  n <- length(th)
  na_out <- structure(list(n = n, mean_direction_h = NA_real_, R = NA_real_,
    rayleigh_stat = NA_real_, rayleigh_p = NA_real_, hr_stat = NA_real_,
    hr_p = NA_real_), class = c("actiRhythm_phasetest", "list"))
  if (n < 3L) return(na_out)

  C <- sum(cos(th)); S <- sum(sin(th))
  R <- sqrt(C^2 + S^2) / n
  mean_dir <- (atan2(S, C) %% (2 * pi)) / (2 * pi) * period

  # Rayleigh test (Fisher 1993): Z = nR^2 with small-sample correction.
  Z <- n * R^2
  p_ray <- exp(-Z) * (1 + (2 * Z - Z^2) / (4 * n) -
                      (24 * Z - 132 * Z^2 + 76 * Z^3 - 9 * Z^4) / (288 * n^2))
  p_ray <- min(max(p_ray, 0), 1)

  # Hermans-Rasson test (Landler et al. 2019). The base statistic is the pairwise
  # angular-distance term plus the |sin| term; clustered angles make it small, so
  # the reported V = (n^2 pi - T) / n increases with clustering and its upper-tail
  # permutation p equals the lower tail of T.
  hr_T <- function(a) {
    d <- outer(a, a, "-")
    sum(pi - abs(pi - (d %% (2 * pi)))) + 2.895 * sum(abs(sin(d)))
  }
  hr_V <- function(a) (n^2 * pi - hr_T(a)) / n
  Vobs <- hr_V(th)
  perm <- vapply(seq_len(n_perm), function(i)
    hr_V(stats::runif(n, 0, 2 * pi)), numeric(1))
  p_hr <- (1 + sum(perm >= Vobs)) / (n_perm + 1)

  structure(list(n = n, mean_direction_h = mean_dir, R = R,
    rayleigh_stat = Z, rayleigh_p = p_ray, hr_stat = Vobs, hr_p = p_hr),
    class = c("actiRhythm_phasetest", "list"))
}

#' @export
print.actiRhythm_phasetest <- function(x, ...) {
  cat("Phase Concentration Tests\n\n")
  cat(sprintf("  n days:          %d\n", x$n))
  cat(sprintf("  Mean direction:  %05.2f h    R: %.3f\n", x$mean_direction_h, x$R))
  cat(sprintf("  Rayleigh:        Z = %.2f, p = %.4f\n", x$rayleigh_stat, x$rayleigh_p))
  cat(sprintf("  Hermans-Rasson:  V = %.2f, p = %.4f\n\n", x$hr_stat, x$hr_p))
  invisible(x)
}
