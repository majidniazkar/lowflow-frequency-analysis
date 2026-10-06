# Reading daily discharge, the D-day mean, and block-minimum extraction.
#
# The D-day mean is rolled over a GAP-FREE daily calendar, and a window that
# contains a missing day is NA rather than an average of whichever days happen
# to be present. Blocks are screened on coverage before contributing a minimum.

suppressPackageStartupMessages(library(zoo))

MONTH_ABBR <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

season_label <- function(months) {
  if (length(months) == 0) return("year")
  months <- as.integer(months)
  contiguous <- length(months) == 1 || all(diff(months) == 1)
  if (contiguous) paste(MONTH_ABBR[months], collapse = "") else paste(months, collapse = "-")
}

#' Read a daily discharge record.
#'
#' @param path .xlsx or .csv with one row per day.
#' @param name station label used in titles and tables.
#' @param D width of the moving-average window in days. 7 is the MAM7 / 7Q10
#'   convention; D = 1 is the single-day minimum, noisier and more sensitive to
#'   gauging error at low stage.
#' @param zero_policy "keep" (default), "drop" or "censor". A recorded zero is
#'   the most extreme low flow a gauge can report, so keeping it is the default;
#'   dropping zeros biases the low tail upward.
#' @return list with `daily` (date, flow, flow_D), `name`, `D` and counts.
read_flow <- function(path, name = "station", D = 7, zero_policy = "keep",
                      date_col = "date", flow_col = "flow") {
  stopifnot(zero_policy %in% c("keep", "drop", "censor"))
  df <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.frame(readxl::read_excel(path))
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE)
  }
  for (cl in c(date_col, flow_col)) {
    if (!cl %in% names(df)) {
      stop(sprintf("column '%s' not found; columns are: %s", cl,
                   paste(names(df), collapse = ", ")))
    }
  }

  # as.Date() keeps the local wall-clock date and discards any UTC offset, so
  # two records exported with different timezone conventions still align.
  dt <- df[[date_col]]
  dates <- if (inherits(dt, "Date")) dt else as.Date(substr(as.character(dt), 1, 10))
  flow <- suppressWarnings(as.numeric(df[[flow_col]]))
  keep <- !is.na(dates)
  dates <- dates[keep]
  flow <- flow[keep]
  ord <- order(dates)
  dates <- dates[ord]
  flow <- flow[ord]

  n_zero <- sum(flow == 0, na.rm = TRUE)
  if (zero_policy == "drop") {
    flow[!is.na(flow) & flow == 0] <- NA_real_
  } else if (zero_policy == "censor") {
    pos <- flow[!is.na(flow) & flow > 0]
    if (length(pos)) flow[!is.na(flow) & flow == 0] <- min(pos) / 2
  }

  # regularise onto a complete daily calendar BEFORE rolling
  full <- seq(min(dates), max(dates), by = "day")
  f <- rep(NA_real_, length(full))
  f[match(dates, full)] <- flow

  flow_D <- if (D <= 1) f else
    as.numeric(zoo::rollapply(f, D, mean, align = "right", fill = NA_real_))

  rec <- list(
    daily = data.frame(date = full, flow = f, flow_D = flow_D),
    name = name, D = D, zero_policy = zero_policy,
    n_days = length(full), n_valid = sum(!is.na(f)),
    n_missing = sum(is.na(f)), n_zero = n_zero
  )
  class(rec) <- "flow_record"
  rec
}

summary.flow_record <- function(object, ...) {
  d <- object$daily
  cat(sprintf("%s\n", object$name))
  cat(sprintf("  span        : %s to %s (%d calendar days)\n",
              min(d$date), max(d$date), object$n_days))
  cat(sprintf("  valid daily : %d   missing: %d   zero flows: %d\n",
              object$n_valid, object$n_missing, object$n_zero))
  cat(sprintf("  D           : %d-day mean, %d days with a complete window\n",
              object$D, sum(!is.na(d$flow_D))))
  cat(sprintf("  zero policy : %s\n", object$zero_policy))
  invisible(object)
}

.block_year <- function(dates, start_month = 1L) {
  y <- as.integer(format(dates, "%Y"))
  m <- as.integer(format(dates, "%m"))
  if (start_month > 1L) y[m < start_month] <- y[m < start_month] - 1L
  y
}

#' Block minima with an explicit coverage screen.
#'
#' A block contributes its minimum only if `min_coverage` of the calendar days
#' it spans carry a valid D-day mean, so a year with 30 usable days cannot stand
#' beside a complete one. Rejected blocks are returned as well, which keeps a
#' gauge outage visible instead of silently absent.
#'
#' @param start_month 1 for calendar years. Where the low-flow season straddles
#'   the new year, set this so no drought is split across two blocks.
block_minima <- function(record, season = NULL, min_coverage = 0.9,
                         start_month = 1L, column = "flow_D") {
  d <- record$daily
  x <- d[[column]]
  yr <- .block_year(d$date, start_month)
  mo <- as.integer(format(d$date, "%m"))
  insea <- if (is.null(season)) rep(TRUE, nrow(d)) else mo %in% as.integer(season)

  rows <- lapply(sort(unique(yr)), function(y) {
    sel <- yr == y & insea
    if (!any(sel)) return(NULL)
    xv <- x[sel]
    dv <- d$date[sel]
    n_possible <- length(xv)
    n_valid <- sum(!is.na(xv))
    cov <- if (n_possible > 0) n_valid / n_possible else 0
    if (n_valid == 0) return(NULL)
    i <- which.min(xv)
    data.frame(year = y, value = xv[i], date = dv[i], n_possible = n_possible,
               n_valid = n_valid, coverage = cov, used = cov >= min_coverage)
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab)) tab <- data.frame(year = integer(), value = numeric(),
                                      date = as.Date(character()), n_possible = integer(),
                                      n_valid = integer(), coverage = numeric(),
                                      used = logical())
  used <- tab[tab$used, , drop = FALSE]
  list(table = tab, used = used, rejected = tab[!tab$used, , drop = FALSE],
       values = used$value, years = used$year, season = season,
       min_coverage = min_coverage, name = record$name, D = record$D, column = column)
}

#' Per block, whether the annual minimum fell inside the window.
season_diagnostic <- function(record, season, min_coverage = 0.9, start_month = 1L,
                              column = "flow_D") {
  ann <- block_minima(record, NULL, min_coverage, start_month, column)$used
  if (!nrow(ann)) return(data.frame())
  m <- as.integer(format(ann$date, "%m"))
  data.frame(year = ann$year, annual_min = ann$value, month_of_min = m,
             in_season = m %in% as.integer(season))
}

#' Rank every contiguous month window by capture rate.
#'
#' `capture_rate` is the fraction of annual minima falling inside the window.
#' Below about 0.8 the window is not describing the annual drought, whatever
#' care goes into the rest of the analysis.
scan_seasons <- function(record, widths = 3:5, min_coverage = 0.9,
                         start_month = 1L, column = "flow_D") {
  ann <- block_minima(record, NULL, min_coverage, start_month, column)$used
  mm <- as.integer(format(ann$date, "%m"))
  out <- list()
  for (w in widths) {
    for (s in 1:12) {
      months <- ((s - 1 + seq_len(w) - 1) %% 12) + 1
      bm <- block_minima(record, months, min_coverage, start_month, column)
      out[[length(out) + 1]] <- data.frame(
        window = season_label(months),
        months = paste(months, collapse = ","),
        width = w,
        capture_rate = if (length(mm)) mean(mm %in% months) else NA_real_,
        n_blocks_used = nrow(bm$used),
        median_min = if (length(bm$values)) stats::median(bm$values) else NA_real_
      )
    }
  }
  tab <- do.call(rbind, out)
  tab[order(-tab$capture_rate, tab$width), ]
}

#' The month window the record itself supports, as an integer vector.
#'
#' Among the candidates whose capture rate is within `tol` of the best, take the
#' NARROWEST window. The tolerance matters because capture rate is monotone in
#' window width -- a wider window is a superset and can only contain more minima
#' -- so ranking on capture rate alone would almost always return the widest
#' candidate offered, which defeats the purpose of a seasonal analysis. A
#' 4-month window capturing 100% of minima describes the drought season better
#' than a 5-month one capturing the same 100%. Use `tol = 0` for the strict
#' highest-capture-rate window.
#'
#' Stops if the best candidate is below `min_capture`: no contiguous season then
#' describes the drought, and annual blocks are the honest alternative.
best_season <- function(record, widths = 3:5, min_coverage = 0.9, start_month = 1L,
                        column = "flow_D", min_capture = 0.8, tol = 0.02) {
  scan <- scan_seasons(record, widths, min_coverage, start_month, column)
  scan <- scan[!is.na(scan$capture_rate), , drop = FALSE]
  if (!nrow(scan)) stop("season scan produced no usable candidate window")
  best <- max(scan$capture_rate)
  near <- scan[scan$capture_rate >= best - tol, , drop = FALSE]
  near <- near[order(near$width, -near$capture_rate), ]
  top <- near[1, ]
  if (top$capture_rate < min_capture) {
    stop(sprintf(paste0("the best window (%s, capture rate %.2f) still misses the annual ",
                        "minimum in %.0f%% of blocks, below min_capture = %.2f. No contiguous ",
                        "season describes the drought on this record -- use annual blocks ",
                        "(season = NULL), or lower min_capture deliberately."),
                 top$window, top$capture_rate, 100 * (1 - top$capture_rate), min_capture))
  }
  as.integer(strsplit(top$months, ",")[[1]])
}
