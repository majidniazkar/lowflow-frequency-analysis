# The R implementation

A self-contained R implementation of the same analysis as the Python package,
following the same conventions. Either language can be used on its own; running
both on one record is a useful cross-check, since the two share no code.

## Install and run

```r
install.packages(c("readxl", "writexl", "zoo", "lmomco", "fitdistrplus",
                   "copula", "Kendall", "ggplot2"))
```

From the repository root:

```bash
Rscript scripts_R/00_make_synthetic_data.R   # writes data/{station_a,station_b,projection}.xlsx
Rscript scripts_R/01_marginal_annual.R
Rscript scripts_R/04_copula_annual.R
```

Or interactively:

```r
source("R/lowflow.R")
lowflow_check()          # version, and whether every dependency is present
```

Each script writes to its own directory under `out_R/`: CSV tables, PNG and PDF
figures, and a `run_log.txt` transcribing every decision and diagnostic. Read
the log — it carries the coverage rejections, the season choice and the
tail-dependence check.

## Layout

```
R/
  lowflow.R     loader; source this one file
  dataio.R      reading, gap-free daily index, D-day means, zero policy,
                block minima with coverage screening, season selection
  marginals.R   GEV (multistart), Pearson III (L-moments), Gamma, Weibull;
                AIC; bootstrap Anderson-Darling; bootstrap intervals;
                trend tests; stationary-vs-trend GEV by likelihood ratio
  copulas.R     Gumbel / Joe / survival Clayton by maximum pseudo-likelihood,
                bootstrap Sn, AND / OR / Kendall joint return periods
  figures.R     the figures, generated from the same fitted objects as the tables
  pipeline.R    the single univariate and single bivariate implementation
scripts_R/
  00_make_synthetic_data.R   generate runnable stand-in data
  01_marginal_annual.R       annual minima, one gauge
  02_marginal_seasonal.R     seasonal minima, one gauge
  03_marginal_projection.R   non-stationary GEV on a projection series
  04_copula_annual.R         joint annual drought, two gauges
  05_copula_seasonal.R       joint seasonal drought, two gauges
```

Six scripts, but one univariate and one bivariate implementation; the scripts
differ only by a configuration list built with `univariate_config()` or
`bivariate_config()`. Produce a variant by passing a different configuration,
never by copying the pipeline — the block definition, the output directory and
the figure labels then exist in exactly one place each and cannot drift apart.

## Conventions

These are the same as on the Python side; `docs/CONVENTIONS.md` states them in
full, with the reasoning and the failure modes they guard against.

- **`D = 7`** by default, the 7-day mean behind MAM7 / 7Q10. The mean is rolled
  over a gap-free daily calendar, and a window containing a missing day is `NA`.
- **Calendar-year blocks**, screened at **`min_coverage = 0.9`**. Rejected
  blocks are reported and drawn as open markers, never silently dropped.
- **Zeros are kept.** A recorded zero is the most extreme low flow a gauge can
  report; dropping zeros biases the low tail upward.
- **The season comes from the record.** `scan_seasons()` ranks every contiguous
  window by capture rate and `best_season()` returns the narrowest window
  within a tolerance of the best, because capture rate is monotone in window
  width. It stops rather than returning a window below a 0.8 capture rate.
- **Marginals are selected by AIC** and tested with a parametric-bootstrap
  Anderson-Darling statistic. A Kolmogorov-Smirnov p-value appears in a column
  named `KS_p_invalid` and must not be used for selection: with parameters
  estimated in-sample its rejection rate at alpha = 0.05 is 0.000.
- **Every reported return level carries a confidence interval** (parametric
  bootstrap), and a warning fires if an interval crosses zero discharge.
- **Copulas are fitted on rank pseudo-observations by maximum
  pseudo-likelihood**, so the dependence estimate cannot inherit a defect of
  the fitted margins.
- **Only upper-tail-dependent copulas are candidates** — Gumbel, Joe and
  survival Clayton. Gaussian, Frank and plain Clayton all have
  `lambda_drought = 0`; they are available for a sensitivity check and are not
  a candidate set.
- **AND, OR and Kendall are all reported and named.** They are three different
  definitions of a "T-year joint drought" and give materially different
  discharges.
- **Report to about twice the record length and no further.**

## Orientation

Minima are negated so that block minima become block maxima, which is what
extreme-value theory describes; the Jacobian is 1, so densities and likelihoods
carry over unchanged. Every marginal models discharge **Q** directly: `cdf(q)`
is non-exceedance, `ppf(p)` returns a discharge, and the T-year low flow is
`ppf(1/T)`. The sign flip is confined to `fit_gev_minima()`.

In `copulas.R`, `u` is a **droughtiness level**, `u = 1 - F_Q(q)`, so large `u`
means deep drought. Joint drought is therefore the **upper** set, with
probability `1 - u - v + C(u, v)`. Discharge is recovered from a droughtiness
level only through `discharge_from_u()`, which is `ppf(1 - u)`.

Note that `lmomco` reports the GEV shape as `kappa = -xi` in the Coles
convention used here, and SciPy's `genextreme` uses `c = -xi` as well. The
conversion is applied once, where the L-moment estimate is used as a starting
value, and nowhere else.

## Differences from the Python implementation

Both implementations follow the conventions above and agree on the quantities
they share. Where they differ:

- **Confidence intervals.** R uses a parametric bootstrap for every family.
  The Python package additionally offers profile-likelihood intervals for the
  GEV.
- **Non-stationary GEV.** Both fit their own trend model. The R side seeds the
  trend fit from the stationary optimum, so a negative deviance is impossible
  at a genuine optimum and is reported as an optimiser failure rather than as
  evidence about the data.
- **Survival Clayton.** Fitted here as a plain Clayton on the reflected ranks
  `1 - u`, which is exactly the survival copula in `u`-space and avoids relying
  on rotation support in the fitting routine.
- **Goodness of fit for copulas** uses `copula::gofCopula(method = "Sn",
  simulation = "pb")`; Python implements the same statistic directly.
- **Figures** are ggplot2 here and Matplotlib there, so they are not
  pixel-identical.
