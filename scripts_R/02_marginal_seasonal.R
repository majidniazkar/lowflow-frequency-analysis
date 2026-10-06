# Univariate low-flow frequency analysis, SEASONAL block minima.
#
# The season is not assumed. Every contiguous 3-to-5-month window is ranked by
# capture rate -- the fraction of annual minima the window actually contains --
# and best_season() takes the narrowest window within a tolerance of the best.
# The full ranking is printed and written to season_scan.csv, so the choice is
# visible and auditable. Pin SEASON to a vector of months only to reproduce an
# earlier analysis. best_season() stops if no window reaches a 0.8 capture rate,
# which means annual blocks (script 01) are the honest answer for that record.

# Run from the repository root:  Rscript scripts_R/02_marginal_seasonal.R
source("R/lowflow.R")

# Point these at your own workbook; the two columns must be date and flow.
PATH <- "data/station_a.xlsx"
STATION <- "Station A"

D <- 7
MIN_COV <- 0.90         # 90% of the days in the window, not of the year
SEASON <- NULL          # NULL = choose from the record; or pin e.g. c(7, 8, 9)

rec <- read_flow(PATH, name = STATION, D = D, zero_policy = "keep")
summary(rec)

scan <- scan_seasons(rec, min_coverage = MIN_COV)
cat("\nSeason scan, best windows by capture rate:\n")
print(head(scan[, setdiff(names(scan), "months")], 6), row.names = FALSE)

season <- if (is.null(SEASON)) best_season(rec, min_coverage = MIN_COV) else SEASON
cat(sprintf("\nSeason in use: %s (%s)\n",
            paste(season, collapse = ", "), season_label(season)))

cfg <- univariate_config(
  path = PATH,
  station = STATION,
  label = sprintf("%s, seasonal minimum %d-day mean, months %s",
                  STATION, D, paste(season, collapse = ", ")),
  outdir = "out_R/02_marginal_seasonal",
  season = season,
  D = D,
  min_coverage = MIN_COV,
  start_month = 1L,
  zero_policy = "keep",
  T = c(2, 5, 10, 20, 50, 100),
  B_gof = 499,
  B_ci = 999,
  scan_season = TRUE
)

res <- run_univariate(cfg)
cat(sprintf("\nDone. %d files in %s\n", length(res$written), res$outdir))
