# Bivariate joint low-flow analysis, two gauges, SEASONAL block minima.
#
# Identical to script 04 except for the block definition, which lives in one
# place and is used for the fit and the figure alike. The season is chosen from
# the record as in script 02, from the gauge named by SEASON_FROM. Both gauges
# must share one season definition: a copula fitted to blocks defined
# differently from its own marginals cannot be reconciled with anything.

# Run from the repository root:  Rscript scripts_R/05_copula_seasonal.R
source("R/lowflow.R")

# Point these at your own workbooks; the two columns must be date and flow.
NAME_X <- "Station A"; PATH_X <- "data/station_a.xlsx"
NAME_Y <- "Station B"; PATH_Y <- "data/station_b.xlsx"

D <- 7
MIN_COV <- 0.90
SEASON <- NULL          # NULL = choose from the record; keep identical to script 02
SEASON_FROM <- "y"      # which gauge the scan is run on: "x" or "y"

season <- SEASON
if (is.null(season)) {
  path <- if (SEASON_FROM == "x") PATH_X else PATH_Y
  name <- if (SEASON_FROM == "x") NAME_X else NAME_Y
  rec <- read_flow(path, name = name, D = D, zero_policy = "keep")
  scan <- scan_seasons(rec, min_coverage = MIN_COV)
  cat(sprintf("Season scan on %s, best windows by capture rate:\n", name))
  print(head(scan[, setdiff(names(scan), "months")], 6), row.names = FALSE)
  season <- best_season(rec, min_coverage = MIN_COV)
}
cat(sprintf("\nSeason in use: %s (%s)   (keep this identical to script 02)\n",
            paste(season, collapse = ", "), season_label(season)))

cfg <- bivariate_config(
  path_x = PATH_X, path_y = PATH_Y,
  name_x = NAME_X, name_y = NAME_Y,
  label = sprintf("%s and %s, seasonal minimum %d-day mean, months %s",
                  NAME_X, NAME_Y, D, paste(season, collapse = ", ")),
  outdir = "out_R/05_copula_seasonal",
  season = season,
  D = D,
  min_coverage = MIN_COV,
  start_month = 1L,
  zero_policy = "keep",
  T = c(2, 5, 10, 20, 50, 100),
  B_gof = 499,
  families = c("Gumbel", "Joe", "Survival Clayton")
)

res <- run_bivariate(cfg)
cat(sprintf("\nDone. %d files in %s\n", length(res$written), res$outdir))
