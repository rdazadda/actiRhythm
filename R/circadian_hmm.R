# Scaled forward-backward (Zucchini, MacDonald & Langrock 2016).
.hmm_fb <- function(allprobs, Gamma, delta) {
  n <- nrow(allprobs); K <- ncol(allprobs)
  la <- lb <- matrix(NA_real_, n, K)
  foo <- delta * allprobs[1, ]; s <- sum(foo); ls <- log(s); foo <- foo / s
  la[1, ] <- log(foo) + ls
  for (i in 2:n) {
    foo <- (foo %*% Gamma) * allprobs[i, ]; s <- sum(foo)
    ls <- ls + log(s); foo <- foo / s; la[i, ] <- log(foo) + ls
  }
  llk <- ls
  lb[n, ] <- 0; foo <- rep(1 / K, K); ls <- log(K)
  for (i in (n - 1):1) {
    foo <- Gamma %*% (allprobs[i + 1, ] * foo); lb[i, ] <- log(foo) + ls
    s <- sum(foo); foo <- foo / s; ls <- ls + log(s)
  }
  list(la = la, lb = lb, llk = llk)
}

.hmm_viterbi <- function(allprobs, Gamma, delta) {
  n <- nrow(allprobs); K <- ncol(allprobs)
  xi <- matrix(0, n, K); foo <- delta * allprobs[1, ]; xi[1, ] <- foo / sum(foo)
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * Gamma, 2, max) * allprobs[i, ]; xi[i, ] <- foo / sum(foo)
  }
  iv <- integer(n); iv[n] <- which.max(xi[n, ])
  for (i in (n - 1):1) iv[i] <- which.max(Gamma[, iv[i + 1]] * xi[i, ])
  iv
}

.hmm_em <- function(y, K, max_iter, tol) {
  n <- length(y)
  mu <- as.numeric(stats::quantile(y, probs = (seq_len(K) - 0.5) / K))
  sigma <- rep(stats::sd(y) / K + 1e-6, K)
  Gamma <- matrix(0.1 / (K - 1), K, K); diag(Gamma) <- 0.9
  delta <- rep(1 / K, K)
  old_llk <- -Inf
  for (it in seq_len(max_iter)) {
    allprobs <- sapply(seq_len(K), function(k) stats::dnorm(y, mu[k], sigma[k]))
    allprobs[allprobs < 1e-300] <- 1e-300
    fb <- .hmm_fb(allprobs, Gamma, delta); llk <- fb$llk
    gamma <- exp(fb$la + fb$lb - llk); gamma <- gamma / rowSums(gamma)
    tG <- matrix(0, K, K)
    for (j in seq_len(K)) for (k in seq_len(K))
      tG[j, k] <- Gamma[j, k] * sum(exp(fb$la[1:(n - 1), j] +
        log(allprobs[2:n, k]) + fb$lb[2:n, k] - llk))
    Gamma <- tG / rowSums(tG)
    delta <- gamma[1, ] / sum(gamma[1, ])
    cs <- colSums(gamma)
    mu <- colSums(gamma * y) / cs
    sigma <- sqrt(colSums(gamma * (y - rep(mu, each = n))^2) / cs) + 1e-6
    if (abs(llk - old_llk) < tol) break
    old_llk <- llk
  }
  ord <- order(mu)                                  # label by mean (lowest = rest)
  list(mu = mu[ord], sigma = sigma[ord], Gamma = Gamma[ord, ord, drop = FALSE],
       delta = delta[ord], llk = llk, gamma = gamma[, ord, drop = FALSE], n_iter = it)
}

#' State-Space (Hidden Markov) Rest-Activity Model
#'
#' Fits an unsupervised Gaussian hidden Markov model to the activity counts,
#' inferring latent rest and active states and the rhythm with which the subject
#' moves between them (Huang et al. 2018). A threshold scorer labels each
#' epoch independently; the HMM uses the persistence of states. Its decoded
#' path gives a 24-hour state-occupancy profile, the probability of being at rest
#' across the day.
#'
#' @param counts Numeric activity vector (a coarse epoch is recommended for speed).
#' @param timestamps POSIXct timestamps, one per value.
#' @param states Number of hidden states (default 2 = rest/active; 3 adds a
#'   moderate state).
#' @param transform Variance-stabilizing transform of the counts: \code{"sqrt"}
#'   (default), \code{"log"}, or \code{"none"}.
#' @param max_iter,tol EM iteration cap and log-likelihood tolerance.
#' @param n_starts Random restarts; the highest-likelihood fit is kept (default 5).
#' @param seed Optional seed for the restarts.
#'
#' @return An object of class \code{actiRhythm_hmm}: the per-state emission means
#'   and SDs, the transition matrix, the decoded \code{state_path} and a
#'   \code{sleep_state} ("S"/"W") vector, a 24-hour state-occupancy profile, and
#'   the log-likelihood with AIC/BIC. Never errors.
#'
#' @references
#' \insertRef{huang2018hmm}{actiRhythm}
#'
#' @seealso \code{\link{sleep.changepoints}}, \code{\link{sleep.cole.kripke}}
#'
#' @examples
#' ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
#' h  <- as.numeric(format(ts, "%H"))
#' counts <- ifelse(h >= 23 | h < 7, 2, 250) + pmax(0, stats::rnorm(length(ts), 0, 10))
#' rest.hmm(counts, ts)
#'
#' @export
rest.hmm <- function(counts, timestamps, states = 2L, transform = c("sqrt", "log", "none"),
                     max_iter = 200L, tol = 1e-6, n_starts = 5L, seed = NULL) {
  transform <- match.arg(transform)
  if (length(counts) != length(timestamps))
    stop("counts and timestamps must have same length")
  x <- suppressWarnings(as.numeric(counts))
  ok <- is.finite(x) & !is.na(timestamps); x <- x[ok]; ts <- timestamps[ok]
  o <- order(ts); x <- x[o]; ts <- ts[o]
  n <- length(x)
  na_out <- structure(list(states = states, emission = data.frame(),
    transition = matrix(nrow = 0, ncol = 0), state_path = integer(0),
    sleep_state = character(0), tod_profile = data.frame(), loglik = NA_real_,
    AIC = NA_real_, BIC = NA_real_, insufficient = TRUE),
    class = c("actiRhythm_hmm", "list"))
  if (n < 2L * states + 10L || stats::sd(x) == 0) return(na_out)

  y <- switch(transform, sqrt = sqrt(pmax(x, 0)), log = log(x + 1), none = x)
  if (stats::sd(y) == 0) return(na_out)
  if (!is.null(seed)) set.seed(seed)

  best <- NULL
  for (st in seq_len(n_starts)) {
    yj <- if (st == 1L) y else y + stats::rnorm(n, 0, stats::sd(y) / 20)
    fit <- tryCatch(.hmm_em(yj, states, max_iter, tol), error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$llk) && (is.null(best) || fit$llk > best$llk))
      best <- fit
  }
  if (is.null(best)) return(na_out)

  allprobs <- sapply(seq_len(states), function(k) stats::dnorm(y, best$mu[k], best$sigma[k]))
  allprobs[allprobs < 1e-300] <- 1e-300
  path <- .hmm_viterbi(allprobs, best$Gamma, best$delta)
  sleep_state <- ifelse(path == 1L, "S", "W")        # lowest-mean state = rest/sleep
  p <- states * 2 + states * (states - 1)            # mu, sigma, transition free params
  hour <- as.POSIXlt(ts)$hour
  occ <- tapply(path == 1L, hour, mean)
  tod <- data.frame(hour = as.integer(names(occ)), p_rest = as.numeric(occ))

  structure(list(states = states,
    emission = data.frame(state = seq_len(states),
      mean_transformed = best$mu, sd_transformed = best$sigma,
      label = c("rest", if (states == 3) "moderate", "active")[seq_len(states)]),
    transition = best$Gamma, state_path = path, sleep_state = sleep_state,
    tod_profile = tod, loglik = best$llk,
    AIC = -2 * best$llk + 2 * p, BIC = -2 * best$llk + log(n) * p,
    transform = transform, n_iter = best$n_iter, insufficient = FALSE),
    class = c("actiRhythm_hmm", "list"))
}

#' @export
print.actiRhythm_hmm <- function(x, ...) {
  cat("State-Space Rest-Activity Model (Gaussian HMM)\n\n")
  if (isTRUE(x$insufficient)) { cat("  Insufficient data\n\n"); return(invisible(x)) }
  cat(sprintf("  States: %d   log-likelihood: %.1f   AIC: %.1f\n\n",
              x$states, x$loglik, x$AIC))
  e <- x$emission; e$mean_transformed <- round(e$mean_transformed, 2)
  e$sd_transformed <- round(e$sd_transformed, 2)
  print(e, row.names = FALSE)
  cat(sprintf("\n  Time at rest: %.0f%%   persistence (rest self-transition): %.2f\n\n",
              100 * mean(x$state_path == 1L), x$transition[1, 1]))
  invisible(x)
}
