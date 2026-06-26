# Empty-check + numeric coercion shared by the epoch-level sleep scorers. NA
# epochs are zero-filled for the convolution; how they appear in the output is
# governed by each scorer's na_action argument.
.validate_counts <- function(counts, name = "counts") {
  if (length(counts) == 0) stop(name, " vector is empty")
  counts <- suppressWarnings(as.numeric(counts))
  counts[is.na(counts)] <- 0
  counts
}

#' Cole-Kripke Sleep/Wake Scoring
#'
#' Classifies each epoch as sleep or wake from activity counts with the
#' Cole-Kripke algorithm (Cole et al. 1992), validated for adults on one-minute
#' epochs, with the optional Webster rescoring rules (Webster et al. 1982).
#' Returns the per-epoch \code{sleep_state} that the regularity and
#' locomotor-sleep metrics read (\code{\link{sleep.regularity.index}},
#' \code{\link{lids}}).
#'
#' @param counts Numeric vector of activity counts (vertical axis).
#' @param apply_rescoring Apply Webster's rescoring rules (default \code{TRUE}).
#' @param epoch_length Epoch length in seconds (default 60). The algorithm was
#'   validated on 60-second epochs; other lengths raise a warning.
#' @param na_action How NA-count epochs appear in the output: \code{"na"}
#'   (default) emits \code{NA}, so non-wear gaps are not read as sleep;
#'   \code{"wake"} scores them wake; \code{"zero"} scores them from a zero count.
#'
#' @return Character vector of states, \code{"S"} (sleep) or \code{"W"} (wake),
#'   the same length as \code{counts}.
#'
#' @details
#' The sleep index uses a seven-epoch window (four before, the current epoch, and
#' two after), with counts divided by 100 and capped at 300:
#' \deqn{D = 0.001 (106 P_4 + 54 P_3 + 58 P_2 + 76 P_1 + 230 C + 74 N_1 + 67 N_2)}
#' An epoch is scored sleep when \eqn{D < 1}. Webster's rescoring then re-labels
#' short sleep bouts that follow or are surrounded by sustained wake as wake.
#'
#' @references
#' \insertRef{cole1992}{actiRhythm}
#'
#' \insertRef{webster1982}{actiRhythm}
#'
#' @seealso \code{\link{sleep.sadeh}}, \code{\link{sleep.regularity.index}}
#'
#' @examples
#' \donttest{
#' agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
#' state <- sleep.cole.kripke(agd$axis1)
#' table(state)
#' }
#'
#' @export
sleep.cole.kripke <- function(counts, apply_rescoring = TRUE, epoch_length = 60,
                              na_action = c("na", "wake", "zero")) {
  na_action <- match.arg(na_action)
  was_na <- is.na(suppressWarnings(as.numeric(counts)))
  counts <- .validate_counts(counts)
  if (!is.null(epoch_length) && !is.na(epoch_length) && epoch_length != 60) {
    warning("Cole-Kripke was validated for 60-second epochs. epoch_length = ",
            epoch_length, "s.")
  }

  scaled <- counts / 100
  scaled[scaled > 300] <- 300
  n <- length(scaled)

  # Seven-epoch weighted window (Cole et al. 1992, Table 2): four epochs before,
  # the current epoch, and two after.
  coef <- c(106, 54, 58, 76, 230, 74, 67) / 1000
  padded <- c(rep(0, 4), scaled, rep(0, 2))
  sleep_index <- numeric(n)
  for (i in seq_len(n)) sleep_index[i] <- sum(coef * padded[i:(i + 6)])

  state <- ifelse(sleep_index < 1, "S", "W")
  if (apply_rescoring) state <- .apply.webster.rescoring(state)
  if (na_action != "zero") state[was_na] <- if (na_action == "na") NA_character_ else "W"
  state
}


# Webster's rescoring rules: short sleep bouts following or bracketed by long
# wake spans are re-labelled wake. Run-length encoded for speed.
.apply.webster.rescoring <- function(sleep.state) {
  n <- length(sleep.state)
  if (n == 0) return(sleep.state)
  rescored <- sleep.state

  r <- rle(sleep.state)
  lens <- r$lengths; vals <- r$values; nb <- length(lens)
  ends <- cumsum(lens); starts <- c(1, ends[-nb] + 1)

  # After a wake span of given length, re-score the following sleep.
  for (b in seq_len(nb)) {
    if (vals[b] == "W" && b < nb && vals[b + 1] == "S") {
      wake_len <- lens[b]; sleep_len <- lens[b + 1]; s0 <- starts[b + 1]
      if (wake_len >= 4)  rescored[s0] <- "W"
      if (wake_len >= 10 && sleep_len >= 3) rescored[s0:min(s0 + 2, n)] <- "W"
      if (wake_len >= 15 && sleep_len >= 4) rescored[s0:min(s0 + 3, n)] <- "W"
    }
  }
  # Short sleep bouts surrounded by long wake on both sides.
  for (b in seq_len(nb)) {
    if (vals[b] == "S") {
      sleep_len <- lens[b]; s0 <- starts[b]; s1 <- ends[b]
      wake_before <- if (b > 1 && vals[b - 1] == "W") lens[b - 1] else 0
      wake_after  <- if (b < nb && vals[b + 1] == "W") lens[b + 1] else 0
      if (sleep_len <= 6  && wake_before >= 15 && wake_after >= 15) rescored[s0:s1] <- "W"
      if (sleep_len <= 10 && wake_before >= 20 && wake_after >= 20) rescored[s0:s1] <- "W"
    }
  }
  rescored
}


#' Sadeh Sleep/Wake Scoring
#'
#' Classifies each epoch as sleep or wake from activity counts with the Sadeh
#' algorithm (Sadeh et al. 1994), validated for children and adolescents on
#' one-minute epochs.
#'
#' @param counts Numeric vector of activity counts (vertical axis).
#' @param epoch_length Epoch length in seconds (default 60). The algorithm was
#'   validated on 60-second epochs; other lengths raise a warning.
#' @param na_action How NA-count epochs appear in the output: \code{"na"}
#'   (default) emits \code{NA}, so non-wear gaps are not read as sleep;
#'   \code{"wake"} scores them wake; \code{"zero"} scores them from a zero count.
#'
#' @return Character vector of states, \code{"S"} (sleep) or \code{"W"} (wake),
#'   the same length as \code{counts}.
#'
#' @details
#' The sleep index uses an eleven-epoch window (five before, the current epoch,
#' and five after), with counts capped at 300:
#' \deqn{SI = 7.601 - 0.065 \cdot AVG - 1.08 \cdot NATS - 0.056 \cdot SD - 0.703 \cdot LG}
#' where \eqn{AVG} is the window mean, \eqn{NATS} the number of epochs with counts
#' in [50, 100), \eqn{SD} the standard deviation over the current and five
#' preceding epochs, and \eqn{LG = \log(\mathrm{count} + 1)}. An epoch is scored
#' sleep when \eqn{SI > -4} (the threshold used by validated ActiGraph
#' implementations).
#'
#' @references
#' \insertRef{sadeh1994}{actiRhythm}
#'
#' @seealso \code{\link{sleep.cole.kripke}}, \code{\link{sleep.regularity.index}}
#'
#' @examples
#' \donttest{
#' agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
#' state <- sleep.sadeh(agd$axis1)
#' table(state)
#' }
#'
#' @export
sleep.sadeh <- function(counts, epoch_length = 60,
                        na_action = c("na", "wake", "zero")) {
  na_action <- match.arg(na_action)
  was_na <- is.na(suppressWarnings(as.numeric(counts)))
  counts <- .validate_counts(counts)
  if (!is.null(epoch_length) && !is.na(epoch_length) && epoch_length != 60) {
    warning("Sadeh was validated for 60-second epochs. epoch_length = ",
            epoch_length, "s.")
  }

  capped <- counts
  capped[capped > 300] <- 300
  n <- length(capped)
  padded <- c(rep(0, 5), capped, rep(0, 5))

  AVG <- as.numeric(stats::filter(padded, rep(1 / 11, 11), sides = 2))[6:(5 + n)]

  nat <- as.numeric(padded >= 50 & padded < 100)
  NATS <- as.numeric(stats::filter(nat, rep(1, 11), sides = 2))[6:(5 + n)]

  # Standard deviation over the current epoch and the five preceding it.
  sumsq <- as.numeric(stats::filter(padded^2, rep(1, 6), sides = 1))[6:(5 + n)]
  sums  <- as.numeric(stats::filter(padded,    rep(1, 6), sides = 1))[6:(5 + n)]
  variance <- (sumsq - sums^2 / 6) / 5
  variance[variance < 0] <- 0
  SD <- sqrt(variance); SD[is.na(SD)] <- 0

  LG <- log(capped + 1)

  sleep_index <- 7.601 - 0.065 * AVG - 1.08 * NATS - 0.056 * SD - 0.703 * LG
  state <- ifelse(sleep_index > -4, "S", "W")
  if (na_action != "zero") state[was_na] <- if (na_action == "na") NA_character_ else "W"
  state
}
