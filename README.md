# lowflow — extreme low-flow (drought) frequency analysis in Python

Univariate and bivariate frequency analysis of annual and seasonal **minimum**
river discharge: block-minimum extraction with coverage screening, marginal
fitting with AIC selection and bootstrap goodness-of-fit, low-flow return
levels with confidence intervals, non-stationary GEV for climate projections,
and bivariate copulas with AND / OR / Kendall joint return periods.

Written as a port of an earlier R implementation, with the defects found in a
review of that implementation corrected. Nothing in it is specific to a basin,
a season or a pair of gauges: station names come from your own workbooks, and
the low-flow season is selected from the record rather than assumed.

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
docs/
  CONVENTIONS.md              house conventions, orientation, and traps that have bitten
  R_COMPARISON.md             cross-validation against R, and when to prefer which
  API.md                      every public function and class
```

Six scripts, but **one** univariate and **one** bivariate implementation; the
scripts differ only by a configuration object. This is deliberate — four of the
defects in the R predecessor were drift between five near-copies of the same
code rather than statistical error. Never copy the pipeline to make a variant.

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

## What changed from the R original

| Defect in the R scripts | Fix here |
|---|---|
| A missing `library()` call aborts both copula scripts on a clean machine | Root-finding is `scipy.optimize.brentq`; no such dependency |
| Pearson III probability transform saturates at 0/1 and crashes the copula fit (16.5% of synthetic 60-year samples) | `pseudo_obs` uses ranks; copulas are fitted by maximum pseudo-likelihood. Parametric margins are kept only for the back-transform to discharge |
| Ranking by a KS p-value computed with in-sample parameters cannot discriminate | `compare()` ranks by AIC; `ad_test()` is Anderson–Darling with a parametric-bootstrap null |
| The projection script overwrites the observational figures | `outdir` is part of the configuration; titles are built from `label` |
| The seasonal copula script plots annual minima against seasonal curves | The block definition exists once, in `season`, and is used for fit and figure alike |
| The normal-copula figure is labelled with the Gumbel return levels | Figures take one table; the runner-up copula gets its own figure |
| AND curves filtered off a grid: zig-zagged and unevenly sampled | `and_curve()` solves `P_AND(u,v) = 1/T` for `v` at each `u` by Brent's method |
| "Best copula" named before any test; `N = 100`; no AIC, no tau | Candidates ranked by AIC and tested at `B = 499`; Kendall tau, Spearman rho and tail dependence reported |
| Figure subtitle states the opposite inequality | Subtitle is generated from the event definition |
| Only the AND scenario, only at `u = v` | AND, OR and Kendall return periods, plus the most-likely realisation on each AND curve |
| No confidence intervals anywhere | Profile likelihood for the GEV, parametric bootstrap otherwise; drawn as a band |
| Zeros dropped, biasing the low tail upward | `ZeroPolicy`: `keep` (default), `censor`, or `drop`; counts reported either way |
| No per-year coverage requirement | `min_coverage` (default 0.90); rejected blocks listed and drawn as open markers |
| D-day mean rolled across calendar gaps, and dead in two of five scripts | `read_flow` regularises to a gap-free daily index *before* rolling, with `min_periods=D`; `D = 7` is the default and live everywhere |
| Low-flow season assumed; capture rate computed but never used | `scan_seasons` ranks every contiguous window by capture rate; the pipeline warns below 80% |
| Stationary fit to a transient projection | `GEVMinimaNS` fits a trend in the location parameter; `lr_test` compares it to the stationary fit; fixed epochs are fitted alongside |

Two things the R scripts did that were **correct** and are preserved: negating
the minima so that block minima become block maxima (the Jacobian is 1, so
densities and likelihoods carry over unchanged), and the AND-event algebra
`1 - u - v + C(u, v)`, which is the joint *drought* probability when `u` is a
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

If this code supports a publication, please cite the repository. Fill in your
own details:

```
<Majid Niazkar>. lowflow: extreme low-flow frequency analysis in Python. 2026.
https://github.com/majidniazkar/lowflow-frequency-analysis
```

## Licence

MIT — see [LICENSE](LICENSE).
