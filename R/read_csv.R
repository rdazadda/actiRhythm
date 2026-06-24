# Coerce a time column to POSIXct from POSIXct, character, or numeric (epoch sec).
.coerce_time <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = tz))
  x <- as.character(x)
  for (fmt in c("%Y-%m-%d %H:%M:%S", "%m/%d/%Y %H:%M:%S", "%d/%m/%Y %H:%M:%S",
                "%Y-%m-%dT%H:%M:%S")) {
    t <- as.POSIXct(strptime(x, fmt, tz = tz))
    if (!all(is.na(t))) return(t)
  }
  as.POSIXct(x, tz = tz)
}

#' Read an ActiGraph (ActiLife) Epoch CSV
#'
#' Reads an ActiLife epoch (count) CSV export. The metadata header is parsed for
#' the start time and epoch length, the count columns are recognized by name or
#' position, and a regular timestamp grid is built. The result is a data frame
#' with a POSIXct \code{timestamp} and the available count columns, the same shape
#' as \code{\link{agd.counts}}, so it feeds the analysis functions directly.
#'
#' @param filepath Path to the ActiLife CSV.
#' @param tz Time zone for the timestamps (default \code{"UTC"}; ActiLife stores
#'   local clock time without a zone).
#' @param date_format Date format of the header Start Date (default
#'   \code{"\%m/\%d/\%Y"}, the ActiLife M/d/yyyy default).
#' @param epoch_length Fallback epoch length in seconds if the header lacks an
#'   Epoch Period line (default 60).
#'
#' @return A data frame with \code{timestamp} (POSIXct) and the count columns
#'   present (\code{axis1}, \code{axis2}, \code{axis3}, \code{vm}, \code{steps},
#'   \code{lux}).
#'
#' @seealso \code{\link{counts.from.data.frame}}, \code{\link{read.agd}}
#'
#' @export
read.actigraph.csv <- function(filepath, tz = "UTC", date_format = "%m/%d/%Y",
                               epoch_length = 60) {
  if (!file.exists(filepath)) stop("file not found: ", filepath)
  H <- readLines(filepath, n = 20, warn = FALSE)
  grab <- function(pat) { ln <- grep(pat, H, value = TRUE, ignore.case = TRUE)[1]; ln }
  st <- grab("Start Time");  sd <- grab("Start Date");  ep <- grab("Epoch Period")
  start_time <- if (!is.na(st)) trimws(sub("(?i)^.*Start Time\\s+", "", st, perl = TRUE)) else "00:00:00"
  start_date <- if (!is.na(sd)) trimws(sub("(?i)^.*Start Date\\s+", "", sd, perl = TRUE)) else NA_character_
  if (!is.na(ep)) {
    es <- trimws(sub("(?i)^.*\\)\\s*", "", ep, perl = TRUE))
    p <- as.numeric(strsplit(es, ":")[[1]])
    if (length(p) == 3 && all(is.finite(p))) epoch_length <- p[1] * 3600 + p[2] * 60 + p[3]
  }

  # data starts at the first non-metadata, non-empty line
  meta <- "ActiGraph|Serial Number|Start Time|Start Date|Epoch Period|Download|Memory Address|Battery Voltage|^[- ]*$"
  ds <- which(!grepl(meta, H, ignore.case = TRUE) & nzchar(trimws(H)))[1]
  if (is.na(ds)) ds <- length(H) + 1L
  has_header <- ds <= length(H) &&
    grepl("Axis|Date|Time|Vector|Steps|Lux|Inclinometer", H[ds], ignore.case = TRUE)

  body <- utils::read.csv(filepath, header = has_header, skip = ds - 1L,
                          stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(gsub("[ _]", "", names(body)))

  out <- data.frame(row.names = seq_len(nrow(body)))
  pick <- function(key) { i <- which(nm == key); if (length(i)) suppressWarnings(as.numeric(body[[i[1]]])) else NULL }
  if (has_header) {
    for (k in c("axis1", "axis2", "axis3", "steps", "lux")) {
      v <- pick(k); if (!is.null(v)) out[[k]] <- v
    }
    vm <- pick("vectormagnitude"); if (is.null(vm)) vm <- pick("vm")
    if (!is.null(vm)) out$vm <- vm
  } else {
    axes <- c("axis1", "axis2", "axis3", "steps")
    numcols <- body[, vapply(body, is.numeric, logical(1)), drop = FALSE]
    for (j in seq_len(min(ncol(numcols), length(axes)))) out[[axes[j]]] <- numcols[[j]]
  }
  if (is.null(out$axis1)) stop("no recognizable count column found in ", basename(filepath))
  if (is.null(out$vm) && all(c("axis1", "axis2", "axis3") %in% names(out)))
    out$vm <- round(sqrt(out$axis1^2 + out$axis2^2 + out$axis3^2))

  t0 <- as.POSIXct(strptime(paste(start_date, start_time),
                            paste(date_format, "%H:%M:%S"), tz = tz))
  timestamp <- seq(t0, by = epoch_length, length.out = nrow(out))
  cbind(timestamp = timestamp, out)
}

#' Activity Counts from a Data Frame
#'
#' A device-neutral entry point: pull a count column (and optionally a time
#' column) out of any data frame and return a tidy \code{timestamp}/\code{counts}
#' frame ready for the analysis functions. Use it for non-ActiGraph counts or any
#' pre-extracted series.
#'
#' @param df A data frame.
#' @param count_col Name or index of the count column (default \code{"axis1"}).
#' @param time_col Name or index of the timestamp column; if \code{NULL},
#'   timestamps are synthesized from \code{start} by \code{epoch_length} (with a
#'   warning).
#' @param epoch_length Epoch length in seconds for synthesized timestamps
#'   (default 60).
#' @param tz Time zone (default \code{"UTC"}).
#' @param start Start time for synthesized timestamps (default
#'   \code{"1970-01-01"}).
#'
#' @return A data frame with \code{timestamp} (POSIXct) and \code{counts}.
#'
#' @seealso \code{\link{read.actigraph.csv}}
#'
#' @examples
#' df <- data.frame(activity = c(0, 50, 300), clock = c("2024-01-01 00:00:00",
#'   "2024-01-01 00:01:00", "2024-01-01 00:02:00"))
#' counts.from.data.frame(df, count_col = "activity", time_col = "clock")
#'
#' @export
counts.from.data.frame <- function(df, count_col = "axis1", time_col = NULL,
                                   epoch_length = 60, tz = "UTC",
                                   start = "1970-01-01") {
  if (!is.data.frame(df)) stop("df must be a data.frame")
  resolve <- function(col) {
    if (is.numeric(col)) return(col)
    i <- match(tolower(col), tolower(names(df)))
    if (is.na(i)) stop("column '", col, "' not found; available: ",
                       paste(names(df), collapse = ", "))
    i
  }
  counts <- suppressWarnings(as.numeric(df[[resolve(count_col)]]))
  if (!is.null(time_col)) {
    ts <- .coerce_time(df[[resolve(time_col)]], tz)
  } else {
    ts <- seq(as.POSIXct(start, tz = tz), by = epoch_length, length.out = length(counts))
    warning("time_col is NULL; timestamps synthesized from '", start,
            "' by epoch_length = ", epoch_length, "s")
  }
  data.frame(timestamp = ts, counts = counts)
}
