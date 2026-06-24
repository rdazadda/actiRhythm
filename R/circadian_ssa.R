# Diagonal averaging (Hankelization): collapse an L x K matrix back to a
# length-N series by averaging along anti-diagonals.
.ssa_hankelize <- function(Y) {
  ad <- as.vector(row(Y) + col(Y) - 1L)
  as.numeric(rowsum(as.vector(Y), ad)) / as.numeric(tabulate(ad))
}

#' Singular Spectrum Analysis of an Activity Series
#'
#' Decomposes an activity-count series into additive components (trend, a
#' circadian component, ultradian components, and noise) with Basic Singular
#' Spectrum Analysis: embed the series into a Hankel trajectory matrix, take its
#' singular value decomposition, group the resulting elementary series, and
#' reconstruct each group by diagonal averaging. pyActigraphy implements this
#' model-free decomposition for actigraphy.
#'
#' @param counts Numeric activity vector.
#' @param timestamps POSIXct timestamps, one per value.
#' @param window_hours Embedding window length in hours (default 24). The window
#'   in epochs is \code{L = round(window_hours * 3600 / epoch_seconds)}.
#' @param n_components Number of leading elementary components to reconstruct and
#'   keep (default 10).
#' @param groups Optional named list of 1-based component indices, e.g.
#'   \code{list(trend = 1, circadian = 2:3)}. When \code{NULL} the grouping is
#'   chosen automatically (component 1 is the trend; the circadian pair is the
#'   largest component whose period falls in \code{period_range}).
#' @param w_components Number of leading components used for the w-correlation
#'   matrix (default \code{min(n_components, 10)}).
#' @param period_range Two-element period window in hours used to identify the
#'   circadian component (default \code{c(20, 28)}).
#' @param detrend If \code{TRUE}, the trend component is removed from
#'   \code{reconstructed}.
#'
#' @return An object of class \code{actiRhythm_ssa}: the singular values and
#'   partial variances, the reconstructed component series (\code{trend},
#'   \code{circadian}, \code{ultradian}), the w-correlation matrix, the
#'   circadian \code{fundamental_period}, and the share of variance the circadian
#'   component carries. The function never errors; on insufficient data it
#'   returns the same structure with \code{insufficient = TRUE}.
#'
#' @details
#' Singular Spectrum Analysis on a long minute-level series builds a large
#' trajectory matrix; for multi-day recordings, resampling to a coarser epoch
#' (for example 10 to 30 minutes) before calling keeps the decomposition fast.
#'
#' @references
#' Golyandina N, Zhigljavsky A (2013). Singular Spectrum Analysis for Time
#' Series. Springer. \doi{10.1007/978-3-642-34913-3}
#'
#' Vautard R, Yiou P, Ghil M (1992). Singular-spectrum analysis: a toolkit for
#' short, noisy chaotic signals. \emph{Physica D}, 58(1-4):95-126.
#' \doi{10.1016/0167-2789(92)90103-T}
#'
#' Hammad G, Reyt M, Beliy N, et al. (2021). pyActigraphy: open-source python
#' package for actigraphy data visualization and analysis. \emph{PLOS
#' Computational Biology}, 17(10):e1009514. \doi{10.1371/journal.pcbi.1009514}
#'
#' @seealso \code{\link{circadian.period}}, \code{\link{circadian.flm}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 1800, length.out = 7 * 48)
#' h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
#' counts <- 50 + 0.5 * seq_along(ts) + 60 * cos(2 * pi * (h - 14) / 24)
#' ssa <- circadian.ssa(counts, ts)
#' ssa
#'
#' @export
circadian.ssa <- function(counts, timestamps, window_hours = 24,
                          n_components = 10, groups = NULL, w_components = NULL,
                          period_range = c(20, 28), detrend = FALSE) {

  na_result <- function() structure(list(
    window_hours = window_hours, L = NA_integer_, K = NA_integer_, n = 0L,
    epoch_seconds = NA_real_, sigma = numeric(0), d = 0L,
    partial_variances = numeric(0), cumulative_variance = numeric(0),
    components = matrix(numeric(0), 0, 0), wcor = matrix(numeric(0), 0, 0),
    groups = list(), trend = numeric(0), circadian = numeric(0),
    ultradian = numeric(0), fundamental_period = NA_real_,
    circadian_variance = NA_real_, reconstructed = numeric(0),
    residual = numeric(0), span_days = NA_real_, n_used = 0L,
    insufficient = TRUE), class = c("actiRhythm_ssa", "list"))

  if (length(counts) != length(timestamps) || !length(counts)) return(na_result())
  t_sec <- suppressWarnings(as.numeric(timestamps))
  cnt   <- suppressWarnings(as.numeric(counts))
  keep  <- is.finite(cnt) & is.finite(t_sec)
  cnt <- cnt[keep]; t_sec <- t_sec[keep]
  n <- length(cnt)
  if (n < 20L) return(na_result())

  o <- order(t_sec); cnt <- cnt[o]; t_sec <- t_sec[o]
  dt <- diff(t_sec); dt <- dt[dt > 0]
  epoch_seconds <- if (length(dt)) stats::median(dt) else 60
  span_days <- (t_sec[n] - t_sec[1]) / 86400
  t_hours <- (t_sec - t_sec[1]) / 3600

  L <- as.integer(min(max(round(window_hours * 3600 / epoch_seconds), 2L), n - 1L))
  K <- n - L + 1L
  s_sd <- stats::sd(cnt)
  if (n < 2L * L || span_days < 2 || !is.finite(s_sd) || s_sd == 0) return(na_result())
  if (as.double(L) * K > 5e7) return(na_result())

  A <- outer(seq_len(L), seq_len(K), function(i, j) cnt[i + j - 1L])
  nc <- min(n_components, L, K)
  s <- svd(A, nu = nc, nv = nc)
  sigma <- s$d
  total_var <- sum(A * A)
  lambda_s <- sigma^2 / total_var
  d <- sum(sigma > .Machine$double.eps * sigma[1] * max(L, K))
  nc <- max(1L, min(nc, length(sigma), d))

  components <- vapply(seq_len(nc),
    function(r) .ssa_hankelize(s$d[r] * tcrossprod(s$u[, r], s$v[, r])), numeric(n))
  if (is.null(dim(components))) components <- matrix(components, ncol = nc)

  recon_group <- function(idx) {
    idx <- idx[idx >= 1 & idx <= nc]
    if (!length(idx)) rep(0, n) else rowSums(components[, idx, drop = FALSE])
  }

  if (is.null(groups)) {
    circ_cand <- integer(0)
    for (r in seq.int(2L, nc)) {
      lsp <- tryCatch(.lomb_scargle(components[, r], t_hours, period_range[1], period_range[2]),
                      error = function(e) NULL)
      pr <- if (is.null(lsp)) NA_real_ else lsp$peak.at[1]
      if (is.finite(pr) && pr >= period_range[1] && pr <= period_range[2])
        circ_cand <- c(circ_cand, r)
    }
    if (length(circ_cand)) {
      r <- circ_cand[which.max(lambda_s[circ_cand])]
      twin <- NA_integer_
      for (cand in c(r - 1L, r + 1L)) {
        if (cand >= 2L && cand <= nc &&
            abs(lambda_s[cand] - lambda_s[r]) / lambda_s[r] < 0.30) { twin <- cand; break }
      }
      circ <- sort(unique(c(r, if (!is.na(twin)) twin else integer(0))))
    } else circ <- intersect(2:3, seq_len(nc))
    ultra <- intersect(setdiff(4:5, circ), seq_len(nc))
    groups <- list(trend = 1L, circadian = circ, ultradian = ultra)
  }

  wc <- if (is.null(w_components)) min(nc, 10L) else min(w_components, nc)
  wts <- pmin(seq_len(n), min(L, K), n - seq_len(n) + 1L)
  Wm <- components[, seq_len(wc), drop = FALSE] * sqrt(wts)
  inner <- crossprod(Wm)
  dn <- sqrt(diag(inner)); dn[dn == 0] <- NA
  wcor <- inner / outer(dn, dn); diag(wcor) <- 1

  circ_series <- recon_group(groups$circadian)
  lsp <- tryCatch(.lomb_scargle(circ_series, t_hours, period_range[1], period_range[2]),
                  error = function(e) NULL)
  fundamental_period <- if (is.null(lsp)) NA_real_ else lsp$peak.at[1]
  ci <- groups$circadian[groups$circadian >= 1 & groups$circadian <= nc]
  circadian_variance <- if (length(ci)) sum(lambda_s[ci]) else NA_real_

  kept <- sort(unique(unlist(groups)))
  reconstructed <- recon_group(kept)
  if (detrend) reconstructed <- reconstructed - recon_group(groups$trend)
  residual <- cnt - recon_group(kept)

  structure(list(
    window_hours = window_hours, L = L, K = K, n = n, epoch_seconds = epoch_seconds,
    sigma = sigma, d = d, partial_variances = lambda_s,
    cumulative_variance = cumsum(lambda_s), components = components, wcor = wcor,
    groups = groups, trend = recon_group(groups$trend), circadian = circ_series,
    ultradian = recon_group(groups$ultradian), fundamental_period = fundamental_period,
    circadian_variance = circadian_variance, reconstructed = reconstructed,
    residual = residual, span_days = span_days, n_used = n, insufficient = FALSE),
    class = c("actiRhythm_ssa", "list"))
}


#' @export
print.actiRhythm_ssa <- function(x, ...) {
  cat("Singular Spectrum Analysis (Basic SSA)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data for SSA\n\n"); return(invisible(x)) }
  cat(sprintf("  Window length L:    %d epochs (%.1f h)\n", x$L, x$window_hours))
  cat(sprintf("  Series length n:    %d (K = %d, span %.1f days)\n", x$n, x$K, x$span_days))
  cat(sprintf("  Components kept:     %d of rank %d\n", ncol(x$components), x$d))
  cat("\n  Variance explained (leading components):\n")
  for (i in seq_len(min(5L, length(x$partial_variances))))
    cat(sprintf("    ET%-2d  lambda = %.4f  (cumulative %.4f)\n",
                i, x$partial_variances[i], x$cumulative_variance[i]))
  cat("\n  Grouping:\n")
  cat(sprintf("    Trend:      components %s\n", paste(x$groups$trend, collapse = ", ")))
  cat(sprintf("    Circadian:  components %s (%.1f%% of variance)\n",
              paste(x$groups$circadian, collapse = ", "), 100 * x$circadian_variance))
  cat(sprintf("    Fundamental period: %.2f h\n", x$fundamental_period))
  cat("\n  Reference: Golyandina and Zhigljavsky (2013)\n\n")
  invisible(x)
}
