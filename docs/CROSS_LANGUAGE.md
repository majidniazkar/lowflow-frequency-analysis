# R and Python side by side

This repository ships two independent implementations of the same analysis. They
share no code, so running both on one record is a genuine cross-check rather than
a restatement. This file records what they agree on, what has to match before
they can be compared at all, and where each language is the more comfortable
place to work.

## Agreement on the marginal estimators

Both were run on the same 50 annual minima of the 7-day mean, under R 4.4.3
(`lmomco` 2.5.1, `fitdistrplus` 1.2) and Python 3.13 (`scipy` 1.18,
`lmoments3` 1.0). The figures below are the record of that comparison; they are
not reproduced by the current synthetic generator, whose seed and record length
have since changed.

| Quantity | R | Python | Agreement |
|---|---|---|---|
| Pearson III L-moments (mu, sigma, gamma) | 3.0169714286, 1.3656645374, 1.0944052397 | identical | exact |
| Pearson III return levels, T = 2…100 | 2.772781, 1.858468, 1.503243, 1.265851, 1.052542, 0.938466 | identical to 6 d.p. | exact |
| Gamma MLE (shape, scale) | 4.83823449, 0.62356894 | 4.83921179, 0.62344273 | 2e-4 relative |
| Weibull MLE (shape, scale) | 2.41209082, 3.40407090 | 2.41227586, 3.40437162 | 1e-4 relative |
| GEV MLE (multistart) | loc -3.39121, scale 1.41747, shape -0.43082, nllh 83.25071 | loc -3.39122, scale 1.41748, shape -0.43082, nllh 83.25071 | 5 d.p. |

`lmoments3` reproduces Hosking's L-moment estimators to the printed precision of
`lmomco`, which is the result that matters most when Pearson III is among the
candidates.

## Why both implementations use a multistart GEV search

A single-start optimiser with a penalty-based constraint is fragile on block
minima, and the failure is silent. On one 50-minimum sample a default
single-start fit returned `shape = -1.626` with an objective of 1e12 -- an
internal infeasibility penalty -- and that parameter vector puts the sample
maximum outside the fitted support, so the resulting "return levels" converge to
a spurious bound at every return period. There is no warning, no error and no
non-convergence flag; the objective value has to be inspected to notice. Given a
sensible starting value the same routine finds the correct optimum immediately.

The sample that triggers it is unremarkable: 50 annual minima with one dominant
low outlier, which is what a drought record looks like.

Both implementations therefore try several starting values with two optimisers
and accept a candidate only if the objective is finite and every observation
lies inside the fitted support. If you cross-check against any other GEV
routine, assert the same three things:

```r
f <- fevd(-x, type = "GEV")
stopifnot(f$results$convergence == 0, f$results$value < 1e6)
z <- 1 + f$results$par["shape"] * (-x - f$results$par["location"]) / f$results$par["scale"]
stopifnot(all(z > 0))   # every observation inside the fitted support
```

The same discipline applies to the non-stationary fit. The stationary GEV is the
trend model with `mu1 = 0`, so the larger model can never have a worse
likelihood at a genuine optimum; the R implementation seeds the trend fit from
the stationary optimum and treats a negative deviance as proof of an optimiser
failure rather than as evidence about the data. A delegated single-start fit on
the bundled projection series returned a trend-model likelihood *worse* than the
nested stationary one (241.76 against 220.98) with a shape of -4.6, which is
exactly the failure this guard catches.

## Agreement on the copula

Once four things are matched -- the block definition, the averaging window `D`,
the coverage screen, and the copula fitting method -- the two implementations
agree on the Gumbel parameter to 0.002% on two-gauge seasonal blocks
(R 1.8086, Python 1.808631).

Those four settings are also the entire explanation for any disagreement, in
roughly that order of likelihood. The sequence below is from one reconciliation
and shows how much each one moves the answer:

| Configuration | Gumbel theta | implied Kendall tau |
|---|---|---|
| ML on parametric probability values, annual blocks, `D = 1`, no coverage screen | 1.1561 | 0.135 |
| MPL on rank pseudo-observations, same blocks | 1.7057 | 0.414 |
| MPL on ranks, scan-selected seasonal blocks, `D = 7`, 90% coverage screen | 1.8086 | 0.447 |

The empirical Kendall tau on those blocks is 0.400, so the first row understates
the dependence roughly threefold: the parametric probability transform saturates
at 0 and 1 whenever a sample extreme sits on a fitted support bound, and the
copula inherits that. This is the single largest discrepancy anyone comparing two
implementations is likely to meet, and it is a choice of estimator rather than a
numerical difference.

**How to isolate a disagreement.** Compare `u` first -- that is pure dependence.
Then compare the discharge at a matched `u` -- that is pure marginal. Then
compare the block count. On the records used here the marginal term ran about
four times the dependence term and the two partly cancelled, so a small net
difference in a return level can conceal two large offsetting ones.

## A worked consequence: Pearson III ranked best while assigning probability zero

On the 50-minimum record above:

| Distribution | k | nllh | AIC | dAIC | bootstrap A² p | **KS p (invalid)** |
|---|---|---|---|---|---|---|
| Gamma | 2 | 83.11 | 170.23 | 0.00 | 0.362 | 0.778 |
| Weibull | 2 | 83.30 | 170.61 | 0.38 | 0.070 | 0.392 |
| GEV | 3 | 83.25 | 172.50 | 2.27 | 0.046 | 0.432 |
| **Pearson III** | 3 | **inf** | **inf** | **inf** | 0.130 | **0.900** |

The L-moment Pearson III fit had skew 1.094, location 3.017 and scale 1.366,
putting a hard lower bound on discharge at 0.521 m3/s. The smallest annual
minimum in the record was 0.350 m3/s -- 0.171 below that bound. The fitted
density at the observed minimum is therefore exactly 0, the log-likelihood is
`-inf`, and AIC is infinite: *the distribution assigns probability zero to an
event in the record.*

The Kolmogorov-Smirnov test ranked it first of the four, p = 0.90. KS sees only
the CDF, which is a perfectly finite 0 at that point; it has no way to see that
the density is zero. Ranking by that p-value would have selected it, and the
reported return levels (2.77, 1.86, 1.50, 1.27, 1.05, 0.94 m3/s for T = 2…100)
carry no sign of the problem. The same saturation is what breaks copula maximum
likelihood.

The lesson is not "never use Pearson III" -- it is standard and often appropriate
for low flows. It is that a three-parameter distribution with a fitted support
bound must be checked against the sample range, and that the check has to be a
likelihood, not an empirical-distribution distance. See `docs/CONVENTIONS.md`
for the mirror-image failure, in which a negative skew estimate leaves the low
tail unbounded instead.

## Which language to work in

R has the broader extreme-value ecosystem; Python has the broader
everything-else. For an analysis of this shape the deciding factor is usually
neither.

Where R is ahead:

- **Breadth of ready-made EVA.** `extRemes`, `evd`, `ismev`, `POT`, `texmex`,
  `lmomco`, `copula`, `fExtremes`. `lmomco` alone covers around 30 distributions
  with L-moment estimators where `lmoments3` covers about a dozen; `copula`
  offers rotated and vine copulas, several goodness-of-fit statistics and several
  estimation methods as library functions.
- **Published provenance.** `lmomco` implements Hosking's FORTRAN; `copula` is
  the Hofert/Kojadinovic/Maechler/Yan reference implementation. Reviewers in
  hydrology recognise them, which is worth something in a methods section.
- **Profile-likelihood intervals** are available off the shelf.

Where Python is ahead:

- **The data layer.** The gap-aware D-day mean, coverage screening and
  common-day intersection are bookkeeping that `pandas` makes hard to get wrong.
- **Climate projections.** Where the projection series comes from CMIP or
  EURO-CORDEX NetCDF, `xarray` plus `dask` is decisively better for extracting
  it, and the low-tail bias adjustment this analysis still needs has mature
  implementations (`xclim.sdba`, `ibicus`).
- **Packaging.** Keeping one implementation behind a configuration object, with
  tests, is the path of least resistance rather than a discipline to maintain.

Neither language fits a GEV better than the other, as the table at the top
shows. The things that actually decide whether a low-flow analysis is right --
the block definition, the coverage screen, the estimator for the dependence,
whether uncertainty is reported, and whether an optimiser failure would be
noticed -- are language-independent, and both implementations here handle them
the same way.
