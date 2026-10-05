# `lowflow` API reference

Import after `lowflow_on_path()` (the kernel sidecar does this for you).

## `lowflow.dataio`

- `read_flow(path, *, name, date_col="date", flow_col="flow", D=1, zero_policy="keep", sheet_name=0) -> FlowRecord`
  Reads `.xlsx`/`.csv`, regularises to a **gap-free daily index before** rolling,
  then computes the D-day right-aligned mean with `min_periods=D` so a window
  containing a gap is NaN rather than an average over fewer days.
- `ZeroPolicy(mode, floor)` — `"keep"` (default), `"censor"`, `"drop"`.
- `FlowRecord.daily` (`flow`, `flow_D`), `.diagnostics`, `.summary()`.

## `lowflow.blocks`

- `block_minima(record, *, season=None, min_coverage=0.9, start_month=1, column="flow_D") -> BlockMinima`
- `BlockMinima.table` (all blocks incl. rejected), `.used`, `.values`, `.years`, `.summary()`
- `season_diagnostic(record, season, ...) -> DataFrame` — per block, whether the
  annual minimum fell inside the window (`in_season`); its mean is the capture rate.
- `best_season(record, *, widths=(3,4,5), min_capture=0.8, tol=0.02, ...) -> tuple`
  — the window the record supports: among candidates within `tol` of the best
  capture rate, the NARROWEST. Raises if the best is below `min_capture`.
- `scan_seasons(record, *, widths=(3,4,5), ...) -> DataFrame` — every contiguous
  window ranked by `capture_rate`.

## `lowflow.marginals`

All classes model discharge Q directly: `cdf(q)` is non-exceedance, `ppf(p)` a
discharge, `return_level(T) == ppf(1/T)`.

- `GEVMinima` — MLE on -Q, Coles convention, 3 starts x 2 optimisers.
  `.params` (loc/scale/shape of -Q), `.lower_bound`, `.return_level_ci(T, method="proflik")`.
  For minima, `shape < 0` means Q has a finite lower bound; `shape > 0` means the
  low tail is unbounded below and long return periods may go negative.
- `PearsonIII` — L-moments (Hosking), equals R's `parpe3`. **Check `.support()`
  and `.notes`**: a support bound above the sample minimum gives zero density
  there and infinite AIC.
- `GammaDist`, `WeibullDist` — MLE with location fixed at 0.
- `GEVMinimaNS.fit(x, covariate)` — linear trend in the location; `loc1 > 0` is
  drying. `.frozen_at(t)`, `.return_level(T, t)`.
- `fit_all(x) -> dict`, `compare(fits, B=499) -> DataFrame` (AIC + bootstrap A^2;
  also prints the invalid KS p for contrast), `ad_test`, `ks_test_naive`, `lr_test`.
- `lmom_ratios(x)`, `theoretical_lmom_curve(family, t3)`, `trend_tests(values, years)`.

## `lowflow.copulas`

`u` is a **droughtiness level** (`u = 1 - F_Q(q)`): large `u` = deep drought.

- `pseudo_obs(*arrays)` — ranks/(n+1). Pass the **negated** minima.
- `discharge_from_u(marginal, u)` — the only sign flip; `q = ppf(1 - u)`.
- `FAMILIES` — Gumbel, Joe, Survival Clayton, Gaussian, Frank, Clayton.
  `lambda_upper` is the **joint-drought** tail dependence; it is 0 for Gaussian,
  Frank and plain Clayton.
- `fit_copula(family, u)` — maximum pseudo-likelihood (ranks only).
- `fit_all_copulas`, `compare_copulas(fits, B=499)` — AIC + bootstrap Cramer-von Mises Sn.
- `dependence_summary(u)`, `pickands_cfg(u, t)`, `lambda_cfg(u)` — model-free.
- `and_probability`, `or_probability`, `and_design_point(model, T, mx, my)`,
  `and_curve(model, T, mx, my, n=300)` (Brent root-finding, exactly on contour),
  `most_likely_design_point`, `kendall_critical_level(model, T, N=40_000)`,
  `kendall_return_period(model, p, N=40_000)`.

## `lowflow.figures`

`apply_style(base=9)` then `fig_block_minima`, `fig_density_fits`, `fig_qq`,
`fig_return_levels`, `fig_lmom_diagram`, `fig_season_scan`, `fig_pseudo_obs`,
`fig_copula_density`, `fig_and_curves`, `plotting_positions(x, a=0.0)`.

## `lowflow.pipeline`

`UnivariateConfig` / `BivariateConfig` dataclasses + `run_univariate` /
`run_bivariate`. One implementation each; variants are configs, never copies.
