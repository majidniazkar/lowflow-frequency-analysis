# lowflow — extreme low-flow (drought) frequency analysis

Univariate and bivariate frequency analysis of annual and seasonal **minimum**
river discharge: block-minimum extraction with coverage screening, marginal
fitting with AIC selection and bootstrap goodness-of-fit, low-flow return
levels with confidence intervals, non-stationary GEV for climate projections,
and bivariate copulas with AND / OR / Kendall joint return periods.

Two independent implementations are included, one in **Python** and one in
**R**, following the same conventions. Either can be used on its own; running
both on one record is a cross-check, since they share no code.

Nothing here is specific to a basin, a season or a pair of gauges: station
names come from your own workbooks, and the low-flow season is selected from
the record rather than assumed.

Not for flood or high-flow analysis. Everything here is oriented around minima.

## Install

```bash
git clone https://github.com/majidniazkar/lowflow-frequency-analysis.git
cd lowflow-frequency-analysis
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .
```

Requires Python 3.10+, `numpy`, `scipy`, `pandas`, `matplotlib`, `statsmodels`,
`openpyxl` and `lmoments3`. For the notebooks: `pip install -e ".[notebooks]"`.

For the R implementation, see [`R/README.md`](R/README.md):

```r
install.packages(c("readxl", "writexl", "zoo", "lmomco", "fitdistrplus",
                   "copula", "Kendall", "ggplot2"))
source("R/lowflow.R")
```

## Quickstart

No discharge records are committed to this repository (see `data/README.md`).
Generate the synthetic stand-in series first, then run any script:

```bash
python scripts/00_make_synthetic_data.py     # writes data/{station_a,station_b,projection}.xlsx
python scripts/01_marginal_annual.py
python scripts/04_copula_annual.py
```

Each script writes to its own directory under `out/`: CSV tables, PNG and PDF
figures, and a `run_log.txt` transcribing every decision and diagnostic. Read
the log — it carries the coverage rejections, the season warning and the
tail-dependence check.

To use your own data, point `path=` at a daily workbook with `date` and `flow`
columns (override with `date_col=` / `flow_col=`) and set the station name to
whatever you call that gauge. The seasonal scripts need no further input: they
rank every contiguous 3-to-5-month window by how often it captures the annual
minimum and use the best one.

Minimal API use:

```python
from lowflow import read_flow, block_minima, scan_seasons, best_season, fit_all, compare

rec = read_flow("data/my_gauge.xlsx", name="My gauge", D=7, zero_policy="keep")
print(rec.summary())
print(scan_seasons(rec).head(8))            # rank windows by capture rate

season = best_season(rec)                   # the top-ranked window, from the data
bm = block_minima(rec, season=season, min_coverage=0.90)
fits = fit_all(bm.values)
print(compare(fits, bm.values))             # AIC table, bootstrap Anderson-Darling
print(fits["Weibull"].ppf(1 / 10))          # the 10-year low flow
```

## Layout

```
lowflow/
  dataio.py     reading, gap-free daily index, D-day means, zero policy, coverage
  blocks.py     block-minimum extraction, coverage screening, season selection
  marginals.py  GEV (stationary and trend), Pearson III, Gamma, Weibull; AIC;
                bootstrap Anderson-Darling; profile-likelihood intervals
  copulas.py    six bivariate families, MPL fitting, bootstrap Sn, AND/OR/Kendall
  figures.py    publication-grade figures, self-contained styling
  pipeline.py   the single univariate and single bivariate implementation
scripts/
  00_make_synthetic_data.py   generate runnable stand-in data
  01_marginal_annual.py       annual minima, one river
  02_marginal_seasonal.py     seasonal minima, one river
  03_marginal_projection.py   non-stationary GEV on a projection series
  04_copula_annual.py         joint annual drought, two rivers
  05_copula_seasonal.py       joint seasonal drought, two rivers
notebooks/                    the same five analyses, narrated, outputs stripped
R/                            the R implementation (see R/README.md)
scripts_R/                    the same six scripts in R
docs/
  CONVENTIONS.md              conventions, orientation, and failure modes to know about
  CROSS_LANGUAGE.md           what R and Python agree on, and which to work in
  API.md                      every public function and class
```

Six scripts per language, but **one** univariate and **one** bivariate
implementation in each; the scripts differ only by a configuration object.
Never copy the pipeline to make a variant — near-duplicate copies of an
analysis drift in ways that are invisible in review and expensive in results,
and none of those failures is a statistics mistake, which is exactly why they
survive. Pass a different configuration instead.

## Read this before interpreting any output

[`docs/CONVENTIONS.md`](docs/CONVENTIONS.md) is the part of this project that is
not code: the sign and orientation conventions (`u` is a *droughtiness* level,
not a non-exceedance probability), the selection rules, and nine failure modes
that have actually occurred on real records — among them a Pearson III fit that
assigns probability zero to an observed minimum while the Kolmogorov–Smirnov
test ranks it best of four candidates, and a timezone mismatch between two
gauges that silently empties the joint analysis.

Three conventions are worth stating on the front page:

- **`u` is droughtiness**, `u = 1 - F_Q(q)`, so joint drought is the *upper*
  set and the admissible copulas are those with *upper*-tail dependence.
  `families` defaults to `("Gumbel", "Joe", "Survival Clayton")` for that
  reason; Gaussian, Frank and plain Clayton all have zero drought tail
  dependence and are a sensitivity check, not a candidate set.
- **Selection is by AIC with a bootstrap Anderson–Darling test.** The
  Kolmogorov–Smirnov p-value is printed, labelled invalid, purely to show why
  it must not be used: with parameters estimated in-sample its rejection rate
  at alpha = 0.05 is 0.000 and its median p-value 0.92.
- **Name the joint event.** AND, OR and Kendall are three different definitions
  of a "T-year joint drought" and give materially different discharges.

## Design notes

The decisions that shape the numbers, in one place:

- **Dependence is estimated on ranks.** Pseudo-observations are rank based and
  copulas are fitted by maximum pseudo-likelihood, so the copula cannot inherit
  a defect of the fitted margins. A parametric probability transform saturates
  at 0 or 1 whenever a sample extreme sits on a fitted support bound, which both
  breaks the fit and shifts the answer when it does not break it.
- **Each gauge gets its own marginal**, selected on its own evidence rather than
  fixed to one family after testing a single record.
- **Selection is by AIC with a parametric-bootstrap Anderson–Darling test** at
  B = 499, so no family is named best before it has been tested.
- **Every reported return level carries a confidence interval**, drawn as a band
  on the return-level figure, and a warning fires when an interval crosses zero
  discharge.
- **Level curves are solved, not filtered.** `and_curve()` finds the matching
  `v` for each `u` by Brent's method, so every point is on the contour to
  machine tolerance and the curve is monotone and evenly sampled.
- **AND, OR and Kendall return periods are all reported**, plus the most-likely
  realisation on each AND contour, and each figure caption states the same
  inequality as its table.
- **The data layer comes first.** The D-day mean is rolled over a gap-free daily
  calendar so a window containing a missing day is NaN; blocks are screened on
  coverage; zeros are kept; and in the bivariate pipeline the screen counts days
  valid at *both* gauges.
- **The season is chosen from the record**, and a projection is analysed on the
  season taken from the observed record at the same gauge so the two describe
  the same event.
- **Non-stationarity is fitted and tested**, not assumed away, with the trend
  fit seeded from the stationary optimum so that a negative deviance is
  recognised as an optimiser failure rather than reported as a result.
- **The GEV sign flip lives in exactly one place.** Minima are negated so that
  block minima become block maxima (the Jacobian is 1, so densities and
  likelihoods carry over unchanged), and the AND-event algebra is
  `1 - u - v + C(u, v)` — the joint *drought* probability when `u` is a
  droughtiness level.

## Limitations

- **No bias adjustment.** `03_marginal_projection.py` fits a projection series
  as supplied. A projected low flow should not be quoted as an absolute 7Q10
  without a low-tail-targeted adjustment against the observed record
  (quantile mapping or delta change applied to the lower tail, not a mean
  correction), nor without stating whether the file is one realisation or one
  member of an ensemble. Until then, report relative change only.
- **Return levels beyond about half the record length** extrapolate well past
  the data. The pipeline warns when an interval crosses zero discharge, which
  means the fit has left physical plausibility — not that negative discharge
  is possible.
- **Kendall return periods are Monte Carlo estimates** from the fitted copula
  (default 40,000 draws) and carry their own sampling error, unlike the
  analytic AND and OR quantities.
- **Bootstrap counts in the scripts are a compromise** (`B_gof = 499`,
  `B_ci = 999`). Use 999 and 1999 for anything published.
- **Block minima are assumed independent between blocks.** The pipeline tests
  this (lag-1 autocorrelation, Ljung–Box) but does not correct for it; a
  significant result means the effective record is shorter than `n` and the
  intervals are too narrow.

## Citing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23182321.svg)](https://doi.org/10.5281/zenodo.23182321)

If this code supports a publication, please cite the archived release. Machine-
readable metadata is in [`CITATION.cff`](CITATION.cff); GitHub renders a "Cite
this repository" button from it.

```
<AUTHOR>. lowflow: extreme low-flow (drought) frequency analysis.
Version 1.1.0, 2026. DOI: [[<10.5281/zenodo.23182321>](https://doi.org/10.5281/zenodo.23182321)](https://doi.org/10.5281/zenodo.23182321)
https://github.com/majidniazkar/lowflow-frequency-analysis
```

Use the **concept** DOI (it always resolves to the latest version) in a paper;
use the version DOI when you need to pin exactly what was run.

## Licence

MIT — see [LICENSE](LICENSE).
