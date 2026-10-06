# Univariate low-flow frequency analysis, ANNUAL block minima.
#
# Annual blocks are the right starting point: they need no assumption about when
# the basin runs dry. Run script 02 if a seasonal window captures the annual
# minimum often enough to be worth using.

# Run from the repository root:  Rscript scripts_R/01_marginal_annual.R
source("R/lowflow.R")

# Point these at your own workbook; the two columns must be date and flow.
PATH <- "data/station_a.xlsx"
STATION <- "Station A"

cfg <- univariate_config(
  path = PATH,
  station = STATION,
  label = sprintf("%s, annual minimum 7-day mean", STATION),
  outdir = "out_R/01_marginal_annual",
  season = NULL,
  D = 7,                  # 7-day mean (MAM7 / 7Q10); D = 1 is the single-day minimum
  min_coverage = 0.90,    # a block needs 90% of its days to contribute a minimum
  start_month = 1L,       # calendar year; appropriate where minima fall in summer
  zero_policy = "keep",
  T = c(2, 5, 10, 20, 50, 100),
  B_gof = 499,            # raise to 999 for publication
  B_ci = 999,             # raise to 1999 for publication
  scan_season = TRUE      # report the season scan even for an annual run
)

res <- run_univariate(cfg)
cat(sprintf("\nDone. %d files in %s\n", length(res$written), res$outdir))
