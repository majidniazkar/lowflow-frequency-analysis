# Conventions, orientation and known traps

This file is the part of the project that is not code: the conventions the
implementations follow, and the failure modes they exist to prevent. Every
item below was encountered on real records; none of it is hypothetical.

## Conventions

- **One implementation, many configurations.** Never copy the pipeline to make
  a variant; pass a different config object. Near-duplicate copies of an
  analysis drift in ways that are invisible in review and expensive in results:
  an output directory still pointing at another run's folder, divergent block
  definitions, a figure labelled from a different fit's table, an averaging
  window live in one copy and dead in another. None of those is a statistics
  mistake, which is exactly why they survive.
- **`D = 7` by default** -- the 7-day mean behind MAM7 / 7Q10. `D = 1` is
  noisier and more sensitive to gauging error at low stage.
- **Calendar-year blocks** are correct where minima fall in late summer, far
  from the year boundary. Use `start_month=` only where the low-flow season
  straddles the new year.
- **Keep zeros** (`zero_policy="keep"`). A zero is the most extreme low flow a
  gauge can record; dropping it biases the low tail upward.
- **Screen blocks by coverage** (`min_coverage=0.90`). A year with 30 valid days
  must not contribute a minimum on equal footing with a full year. In the
  bivariate pipeline the screen counts days valid at **both** gauges.
- **Choose the season from the record, never by assumption.** `scan_seasons`
  reports, for each contiguous month window, the fraction of annual minima the
  window actually contains; `best_season` returns the window to use. Below about
  0.8 the seasonal analysis is not describing the annual drought.
- **Prefer the narrowest window that captures the drought.** Capture rate is
  monotone in window width -- a wider window is a superset and can only contain
  more minima -- so ranking on capture rate alone always returns the widest
  candidate offered, which defeats the purpose of a seasonal analysis. A 4-month
  window capturing 100% of minima is a better description of the drought season
  than a 5-month one capturing the same 100%. `best_season` takes the narrowest
  window within `tol` (default 0.02) of the best capture rate.
- **Select marginals by AIC; test by parametric-bootstrap Anderson-Darling.**
  Never rank candidates by a Kolmogorov-Smirnov p-value computed with
  in-sample parameters -- with parameters estimated from the same data its
  rejection rate at alpha = 0.05 is 0.000 and its median p-value 0.92, so it
  cannot discriminate between families at all.
- **Report confidence intervals on every return level.** Profile likelihood for
  the GEV, parametric bootstrap otherwise.
- **Restrict the copula set to the upper-tail-dependent families** --
  `("Gumbel", "Joe", "Survival Clayton")`. Gaussian, Frank and plain Clayton
  all have `lambda_drought = 0`, asserting that simultaneous extreme drought
  becomes asymptotically independent; that contradicts the shared meteorology
  of a two-river basin, and AIC will select one of them anyway if offered. The
  full six are a sensitivity check, not a candidate set.
- **Report the tail-dependence range, not the winning family.** Where the three
  admissible copulas agree on `lambda_drought` to within a spread much narrower
  than the bootstrap interval on a return level, that range is a stronger
  statement than defending one AIC win.
- **Name the joint event.** AND, OR and Kendall are three different definitions
  of a "T-year joint drought" and give materially different discharges.
- **Stop at T of about twice the record length.** Twenty blocks do not support a
  T = 100 joint return level. The default `T` tuple runs to 100 because record
  lengths differ, not because 100 is always reportable.

## Orientation, stated once

Minima are negated so that block minima become block maxima, which is what
extreme-value theory describes; the Jacobian is 1, so densities and likelihoods
carry over unchanged.

In `lowflow.marginals` every class models discharge **Q** directly: `cdf(q)` is
non-exceedance, `ppf(p)` returns a discharge, and the T-year low flow is
`ppf(1/T)`. The sign flip is confined to `GEVMinima`.

In `lowflow.copulas`, `u` is a **droughtiness level**, `u = 1 - F_Q(q)`, so
large `u` means deep drought. Joint drought is therefore the **upper** set,
with probability `1 - u - v + C(u, v)`, and the copulas that represent
simultaneous drought are the ones with **upper**-tail dependence. Plain Clayton
on these variables is the wrong tail. Discharge is recovered from a
droughtiness level only through `discharge_from_u(marginal, u)`, which is
`marginal.ppf(1 - u)`.

Three libraries disagree on the GEV shape sign: R's `evd` and `extRemes` use
the Coles convention xi, R's `lmomco` reports kappa = -xi, and SciPy's
`genextreme` takes c = -xi. Never pass parameters between them unconverted.

## Traps that have actually bitten

**Pearson III can assign probability zero to an observed minimum.** Its
L-moment fit has a hard support bound. If that bound lands above the sample
minimum, the density there is 0, the log-likelihood is `-inf` and AIC infinite
-- while the KS test ranks it *best* of four candidates, because KS sees only
the CDF. Check `fit.notes` and `fit.support()`. This is also what makes a
parametric probability transform return exactly 0 or 1 and break copula
maximum likelihood.

**Pearson III fails in the other direction too, with an unbounded low tail.**
Fitted by L-moments to negated annual minima, a negative skew estimate places
the support bound at `+inf` on the negated scale -- so on discharge the fitted
Q has no lower bound at all, and its low tail crosses zero at a droughtiness
level inside the range people routinely report. On one 20-year record tested
here that crossing was at `u = 0.9871`, which is the origin of non-credible
T = 50 and T = 100 low flows an order of magnitude below the observed 20-year
minimum. Either sign of
the skew breaks Pearson III at exactly the end of the record you care about.
Prefer a lower-bounded family, and always evaluate the fitted quantile function
at the `u` you intend to report.

**Fit copulas on rank pseudo-observations by maximum pseudo-likelihood, never
on a parametric probability transform by maximum likelihood.** Pushing data
through fitted marginal CDFs first makes the copula inherit every marginal
defect. On one set of annual blocks tested here that route returned a Gumbel theta implying
Kendall tau = 0.135 against a measured empirical tau = 0.400 -- understating
the dependence roughly threefold, because the Pearson III transform was
saturating at 0 and 1. The same blocks through rank pseudo-observations gave
tau = 0.414. `lowflow.copulas.pseudo_obs` is the rank transform; in R it is
`copula::pobs` with `method="mpl"`.

**Mixed timezone awareness between two records silently empties the join.** A
date column exported with a UTC offset parses to a timezone-aware index; a
plain Excel date column does not. Two records that disagree can never align: an
index union degrades to `dtype=object`, the common-day count is 0, and that
propagates as "no usable blocks" and a full set of NaN results that otherwise
looks like a successful run. `read_flow` drops the timezone on entry, and
`run_bivariate` now raises with a diagnostic breakdown rather than writing NaN
files.

**An optimiser can fail silently on block minima with one dominant outlier.**
R's `extRemes::fevd` with default settings returned an infeasible parameter
vector at its 1e12 penalty, with no warning and with observations outside the
fitted support, whose return levels converge to a spurious bound; a 400-start
search on the same data found a shape of -0.43 against the reported -1.63. If
you cross-check against R, assert `f$results$convergence == 0`,
`f$results$value < 1e6`, and that every observation lies inside the fitted
support. This package uses a multi-start, two-optimiser search instead.

**`scipy.stats.multivariate_normal.cdf` integrates point by point.** A Gaussian
copula CDF evaluated over a large Monte Carlo sample is ruinous -- it once
turned a 7-minute script into 10.4 hours. `lowflow.copulas._bvn_cdf` is a
vectorised Gauss-Legendre replacement that matches SciPy to 1.7e-16.

**A confidence interval crossing zero discharge** means the fit is
extrapolating past physical plausibility, not that negative flow is possible.
Prefer a lower-bounded family, or fit in log space, before quoting it.

**Editing files on disk does not update an already-imported module.** After
replacing `lowflow/`, restart the kernel. `lowflow.check_install()` reports the
version, the directory actually imported, a pass/fail line per known fix, and
any other copies of the package it can find on the machine -- worth running
before trusting a result, because `import lowflow` silently takes whichever
copy is first on `sys.path`.

## Still the analyst's job

Bias adjustment of a simulated series against the observed record, targeted at
the **low tail** rather than the mean, before quoting projected 7Q10 values;
and establishing whether a projection file is one realisation or one member of
an ensemble. Nothing in this package supplies either, and without the first a
projected change is reportable only as a relative change, never as an absolute
7Q10.
