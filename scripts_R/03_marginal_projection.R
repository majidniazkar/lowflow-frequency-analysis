# Univariate low-flow frequency analysis of a CLIMATE PROJECTION series.
#
# Three things matter here beyond the observational run:
#
# 1. Its own output directory, so projection results never land under
#    observational filenames and titles.
# 2. Non-stationarity is handled, not merely tested. A transient projection to
#    2100 is the case where a stationary GEV is least defensible, so a trend in
#    the GEV location parameter is fitted and compared against the stationary
#    fit by likelihood ratio, with fixed epochs fitted alongside.
# 3. The season is taken from the OBSERVED record at the same gauge rather than
#    scanned on the projection, so the projected and observed analyses describe
#    the same event.
#
# Still the analyst's job, and not something code can supply: bias adjustment of
# the simulated series against the observed record, targeted at the LOW TAIL
# rather than the mean, and a statement of whether the file is one realisation
# or one member of an ensemble. Until both are settled, report the projected
# change as a RELATIVE change, never as an absolute low-flow value.

# Run from the repository root:  Rscript scripts_R/03_marginal_projection.R
source("R/lowflow.R")

# Point these at your own workbooks; the two columns must be date and flow.
PATH <- "data/projection.xlsx"       # the simulated series
OBS_PATH <- "data/station_b.xlsx"    # the observed record at the same gauge
STATION <- "Station B (projection)"

D <- 7
MIN_COV <- 0.90
SEASON <- NULL         # NULL = take the season from OBS_PATH; or pin e.g. c(7, 8, 9)
EPOCHS <- list(c(2015, 2044), c(2045, 2074), c(2075, 2100))

season <- SEASON
if (is.null(season)) {
  obs <- read_flow(OBS_PATH, name = "observed", D = D, zero_policy = "keep")
  season <- best_season(obs, min_coverage = MIN_COV)
}
cat(sprintf("Season in use: %s (%s) -- from the observed record, so the two are comparable\n",
            paste(season, collapse = ", "), season_label(season)))

cfg <- univariate_config(
  path = PATH,
  station = STATION,
  label = sprintf("%s, seasonal minimum %d-day mean, months %s",
                  STATION, D, paste(season, collapse = ", ")),
  outdir = "out_R/03_marginal_projection",
  season = season,
  D = D,
  min_coverage = MIN_COV,
  start_month = 1L,
  zero_policy = "keep",
  T = c(2, 5, 10, 20, 50, 100),
  B_gof = 499,
  B_ci = 999,
  scan_season = TRUE,
  nonstationary = TRUE,
  epochs = EPOCHS
)

res <- run_univariate(cfg)
cat(sprintf("\nDone. %d files in %s\n", length(res$written), res$outdir))
if (!is.null(res$nonstationary)) {
  cat(sprintf("Trend in GEV location: deviance %.2f on 1 df, p = %.4g\n",
              res$nonstationary$deviance, res$nonstationary$p_value))
}
