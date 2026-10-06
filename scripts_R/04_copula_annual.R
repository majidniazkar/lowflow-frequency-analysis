# Bivariate joint low-flow analysis, two gauges, ANNUAL block minima.
#
# Copulas are fitted by maximum pseudo-likelihood on RANK pseudo-observations,
# so the dependence estimate cannot inherit a defect of the fitted margins.
# Each gauge's marginal is selected on its own evidence. AND, OR and Kendall
# joint return periods are all reported, because they are three different
# definitions of a "T-year joint drought" and give materially different
# discharges; the most likely realisation on each AND contour is reported too.

# Run from the repository root:  Rscript scripts_R/04_copula_annual.R
source("R/lowflow.R")

# Point these at your own workbooks; the two columns must be date and flow.
NAME_X <- "Station A"; PATH_X <- "data/station_a.xlsx"
NAME_Y <- "Station B"; PATH_Y <- "data/station_b.xlsx"

cfg <- bivariate_config(
  path_x = PATH_X, path_y = PATH_Y,
  name_x = NAME_X, name_y = NAME_Y,
  label = sprintf("%s and %s, annual minimum 7-day mean", NAME_X, NAME_Y),
  outdir = "out_R/04_copula_annual",
  season = NULL,
  D = 7,
  min_coverage = 0.90,
  start_month = 1L,
  zero_policy = "keep",
  T = c(2, 5, 10, 20, 50, 100),
  B_gof = 499,
  # Only copulas with upper-tail dependence on the droughtiness scale are
  # admissible: Gaussian, Frank and plain Clayton all have lambda_drought = 0.
  # Add them back as a one-off sensitivity check, never as a candidate set.
  families = c("Gumbel", "Joe", "Survival Clayton")
)

res <- run_bivariate(cfg)
cat(sprintf("\nDone. %d files in %s\n", length(res$written), res$outdir))
