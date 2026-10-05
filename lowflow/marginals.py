"""Marginal distributions for annual/seasonal LOW flows.

Design rule, and the reason this module is not a thin wrapper round scipy: every
class models the **discharge** variable Q directly. ``cdf(q)`` is always the
non-exceedance probability P(Q <= q), ``ppf(p)`` is always a discharge, and the
T-year low flow is always ``ppf(1/T)``. The GEV is still *estimated* on -Q (block
minima become block maxima, which is what extreme-value theory describes), but
that sign flip is confined to :class:`GEVMinima` and never escapes into user
code. In the R scripts the flip was spread across six call sites -- ``fevd`` on
``low_flows_neg``, ``dgev(-x, ...)``, ``-qgev(ppoints(n), ...)``,
``ks.test(low_flows_neg, ...)`` -- which is what made the Q-Q plot correct only
by accident and the return levels easy to get wrong.

Implemented here, addressing the R review:

* **A3** -- :func:`compare` ranks by AIC, and :func:`ad_test` is an
  Anderson-Darling test with a *parametric bootstrap* null, because the
  Kolmogorov-Smirnov p-value with estimated parameters is invalid (measured
  rejection rate 0.000 at alpha=0.05, median p 0.919).
* **B6** -- every fit exposes :meth:`return_level_ci`: profile likelihood for
  the GEV, parametric bootstrap for the others.
* **C5** -- :class:`GEVMinimaNS` fits a trend in the location parameter and
  :func:`lr_test` compares it against the stationary fit, for projection series
  where a stationary fit is indefensible.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.special import gammaln

try:  # L-moment estimators (Hosking); a direct port of the FORTRAN lmoments
    import lmoments3.distr as _ld

    _HAVE_LMOM = True
except ImportError:  # pragma: no cover
    _HAVE_LMOM = False

__all__ = [
    "GEVMinima",
    "GEVMinimaNS",
    "PearsonIII",
    "GammaDist",
    "WeibullDist",
    "fit_all",
    "compare",
    "ad_test",
    "ks_test_naive",
    "lr_test",
    "lmom_ratios",
    "trend_tests",
]

CHI2_1_95 = stats.chi2.ppf(0.95, 1)


# ---------------------------------------------------------------------------
# base
# ---------------------------------------------------------------------------
@dataclass(repr=False)
class Marginal:
    """Base class: a fitted distribution for the low-flow variable Q."""

    params: dict
    data: np.ndarray
    name: str = "marginal"
    n_params: int = 0
    notes: list = field(default_factory=list)

    # -- subclasses implement these three in Q-space ------------------------
    def cdf(self, q):  # P(Q <= q)
        raise NotImplementedError

    def pdf(self, q):
        raise NotImplementedError

    def ppf(self, p):  # discharge with non-exceedance probability p
        raise NotImplementedError

    def rvs(self, n, rng=None):
        rng = np.random.default_rng(rng)
        return self.ppf(rng.uniform(size=n))

    @classmethod
    def fit(cls, x):
        raise NotImplementedError

    # -- shared machinery ---------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.data)

    def nllh(self, x=None) -> float:
        x = self.data if x is None else np.asarray(x, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            lp = np.log(self.pdf(x))
        if not np.all(np.isfinite(lp)):
            return np.inf
        return float(-lp.sum())

    @property
    def aic(self) -> float:
        return 2 * self.nllh() + 2 * self.n_params

    @property
    def bic(self) -> float:
        return 2 * self.nllh() + self.n_params * np.log(self.n)

    def return_level(self, T):
        """T-year low flow: the discharge with non-exceedance probability 1/T."""
        return self.ppf(1.0 / np.asarray(T, float))

    def return_level_ci(self, T, *, alpha=0.05, B=999, rng=None, method="bootstrap"):
        """Confidence interval for the T-year low flow.

        Default is a parametric bootstrap: resample from the fitted model,
        refit the *same family* end to end, and take the empirical quantiles of
        the resulting return levels. :class:`GEVMinima` overrides this with
        profile likelihood, which is better behaved for long return periods.
        """
        T = np.atleast_1d(np.asarray(T, float))
        rng = np.random.default_rng(rng)
        est = np.full((B, T.size), np.nan)
        for b in range(B):
            xb = self.rvs(self.n, rng)
            try:
                est[b] = type(self).fit(xb).return_level(T)
            except Exception:
                continue
        lo, hi = np.nanpercentile(est, [100 * alpha / 2, 100 * (1 - alpha / 2)], axis=0)
        return pd.DataFrame(
            {
                "T": T,
                "return_level": self.return_level(T),
                "lower": lo,
                "upper": hi,
                "n_bootstrap_ok": np.isfinite(est).sum(axis=0),
                "method": f"parametric bootstrap (B={B})",
            }
        )

    def __repr__(self) -> str:  # pragma: no cover
        ps = ", ".join(f"{k}={v:.4g}" for k, v in self.params.items())
        return f"<{self.name}: {ps}, nllh={self.nllh():.3f}, AIC={self.aic:.2f}>"


# ---------------------------------------------------------------------------
# GEV for minima
# ---------------------------------------------------------------------------
def _gev_nllh(theta, x):
    """Negative log-likelihood of GEV *maxima* in the Coles parameterisation.

    F(x) = exp(-(1 + xi (x - mu)/sigma)^(-1/xi)).  This is the same convention
    used by ``extRemes::fevd`` and ``evd::qgev``; note that
    ``scipy.stats.genextreme`` uses ``c = -xi`` and ``lmomco`` uses
    ``kappa = -xi``, so parameters must never be passed between them unconverted.
    """
    mu, logsig, xi = theta
    if not np.isfinite(logsig) or logsig < -30 or logsig > 30 or not np.isfinite(mu):
        return 1e12          # optimiser excursion; keep it out of the warning log
    sig = np.exp(logsig)
    y = (x - mu) / sig
    if abs(xi) < 1e-8:
        return len(x) * logsig + np.sum(y) + np.sum(np.exp(-y))
    z = 1 + xi * y
    if np.any(z <= 1e-12):
        return 1e12
    return len(x) * logsig + (1 + 1 / xi) * np.sum(np.log(z)) + np.sum(z ** (-1 / xi))


def _gev_quantile(mu, sig, xi, p_exceed):
    """Quantile of GEV maxima with exceedance probability ``p_exceed``."""
    y = -np.log1p(-np.asarray(p_exceed, float))  # = -log(1 - p)
    if abs(xi) < 1e-8:
        return mu - sig * np.log(y)
    return mu + (sig / xi) * (y ** (-xi) - 1.0)


def _gev_start(x):
    """L-moment starting values, converted into the Coles convention."""
    if _HAVE_LMOM:
        try:
            p = _ld.gev.lmom_fit(x)
            return np.array([p["loc"], np.log(max(p["scale"], 1e-8)), -p["c"]])
        except Exception:
            pass
    s = np.std(x, ddof=1)
    return np.array([np.mean(x) - 0.45 * s, np.log(0.78 * s), 0.05])


def _fit_gev_mle(x):
    x = np.asarray(x, float)
    best, best_f = None, np.inf
    starts = [_gev_start(x)]
    starts.append(np.array([starts[0][0], starts[0][1], 0.0]))
    starts.append(np.array([np.median(x), np.log(np.std(x, ddof=1)), -0.1]))
    for s0 in starts:
        for meth in ("Nelder-Mead", "L-BFGS-B"):
            try:
                r = optimize.minimize(_gev_nllh, s0, args=(x,), method=meth)
            except Exception:
                continue
            if r.fun < best_f and np.isfinite(r.fun):
                best, best_f = r.x, r.fun
    if best is None:
        raise RuntimeError("GEV MLE failed from every starting value")
    return best[0], np.exp(best[1]), best[2], best_f


@dataclass(repr=False)
class GEVMinima(Marginal):
    """GEV fitted to block **minima**, by maximum likelihood on -Q.

    ``params`` holds the parameters of the distribution of X = -Q
    (``loc``, ``scale``, ``shape``), in the Coles convention, matching
    ``extRemes::fevd``. All public methods are in Q-space.

    For minima the sign of the shape parameter reads the opposite way round to
    the familiar maxima case: ``shape < 0`` here means X is bounded above, i.e.
    Q has a finite **lower** bound -- a physically meaningful lowest possible
    discharge -- while ``shape > 0`` means the low tail is unbounded below and
    will eventually predict negative discharge at long return periods. Check
    :attr:`lower_bound` before quoting a 100-year value.
    """

    name: str = "GEV"
    n_params: int = 3

    @classmethod
    def fit(cls, x):
        x = np.asarray(x, float)
        mu, sig, xi, _ = _fit_gev_mle(-x)
        obj = cls(params=dict(loc=mu, scale=sig, shape=xi), data=x)
        if xi > 0:
            obj.notes.append(
                "shape > 0: low tail unbounded below, long return periods may go negative"
            )
        return obj

    # X = -Q lives in the Coles convention; everything below is Q-space
    def _x(self):
        p = self.params
        return p["loc"], p["scale"], p["shape"]

    def cdf(self, q):
        mu, sig, xi = self._x()
        x = -np.asarray(q, float)
        if abs(xi) < 1e-8:
            t = np.exp(-(x - mu) / sig)
        else:
            z = 1 + xi * (x - mu) / sig
            t = np.where(z > 0, np.maximum(z, 1e-300) ** (-1 / xi), np.where(xi > 0, np.inf, 0.0))
        return 1.0 - np.exp(-t)  # = 1 - F_X(-q)

    def pdf(self, q):
        mu, sig, xi = self._x()
        x = -np.asarray(q, float)
        if abs(xi) < 1e-8:
            y = (x - mu) / sig
            return np.exp(-y - np.exp(-y)) / sig
        z = 1 + xi * (x - mu) / sig
        out = np.zeros_like(np.asarray(q, float))
        ok = z > 0
        zz = np.where(ok, z, 1.0)
        out[ok] = (zz[ok] ** (-1 / xi - 1) * np.exp(-(zz[ok] ** (-1 / xi)))) / sig
        return out

    def ppf(self, p):
        mu, sig, xi = self._x()
        p = np.asarray(p, float)
        return -_gev_quantile(mu, sig, xi, p)  # exceedance p in X-space

    @property
    def lower_bound(self):
        """Finite lower bound on Q implied by the fit, or -inf."""
        mu, sig, xi = self._x()
        return -(mu - sig / xi) if xi < 0 else -np.inf

    # -- profile likelihood CI (R finding B6) --------------------------------
    def _profile_nllh(self, q_T, p_exceed):
        """min over (sigma, xi) of the nllh with the return level fixed at q_T."""
        x = -self.data
        x_T = -q_T
        y = -np.log1p(-p_exceed)

        def obj(sp):
            logsig, xi = sp
            sig = np.exp(logsig)
            if abs(xi) < 1e-8:
                mu = x_T + sig * np.log(y)
            else:
                mu = x_T - (sig / xi) * (y ** (-xi) - 1.0)
            return _gev_nllh(np.array([mu, logsig, xi]), x)

        best = np.inf
        for xi0 in (self.params["shape"], 0.0, -0.2, 0.2):
            try:
                r = optimize.minimize(
                    obj, np.array([np.log(self.params["scale"]), xi0]), method="Nelder-Mead"
                )
                best = min(best, r.fun)
            except Exception:
                continue
        return best

    def return_level_ci(self, T, *, alpha=0.05, B=999, rng=None, method="proflik"):
        """Profile-likelihood CI (default) or parametric bootstrap.

        Profile likelihood is strongly preferred at long return periods: the
        return-level likelihood is markedly asymmetric there and a normal/delta
        interval understates the upper limit.
        """
        if method != "proflik":
            return super().return_level_ci(T, alpha=alpha, B=B, rng=rng, method=method)
        T = np.atleast_1d(np.asarray(T, float))
        crit = stats.chi2.ppf(1 - alpha, 1)
        nllh0 = self.nllh()
        rows = []
        for Ti in T:
            p = 1.0 / Ti
            q_hat = float(self.ppf(p))
            dev = lambda qq: 2 * (self._profile_nllh(qq, p) - nllh0) - crit
            scale = max(abs(q_hat), 1e-6)
            lo = hi = np.nan
            # Q is a low flow: search downward for the lower limit, upward for upper
            for step in (0.05, 0.15, 0.4, 0.8, 1.5, 3.0):
                a = q_hat - step * scale
                if a <= 0 and self.lower_bound == -np.inf:
                    a = q_hat * (1 - step) if q_hat > 0 else q_hat - step * scale
                try:
                    if dev(a) > 0:
                        lo = optimize.brentq(dev, a, q_hat, xtol=1e-8 * scale)
                        break
                except Exception:
                    continue
            for step in (0.05, 0.15, 0.4, 0.8, 1.5, 3.0):
                b = q_hat + step * scale
                try:
                    if dev(b) > 0:
                        hi = optimize.brentq(dev, q_hat, b, xtol=1e-8 * scale)
                        break
                except Exception:
                    continue
            rows.append(dict(T=Ti, return_level=q_hat, lower=lo, upper=hi,
                             method=f"profile likelihood ({100*(1-alpha):.0f}%)"))
        return pd.DataFrame(rows)


@dataclass(repr=False)
class GEVMinimaNS(Marginal):
    """Non-stationary GEV for minima with a linear trend in the location.

    ``mu(t) = mu0 + mu1 * s(t)`` where ``s`` is the covariate standardised to
    zero mean and unit standard deviation, so ``mu1`` is the change in location
    (in units of -Q) per standard deviation of the covariate. Fitted on X = -Q,
    so **``mu1 > 0`` means drying** (the negated minima are trending up).

    Use :func:`lr_test` against the stationary :class:`GEVMinima` to decide
    whether the trend is supported. This is the fix for R finding C5: the
    projection script ran Mann-Kendall and then fitted a stationary model
    regardless.
    """

    name: str = "GEV (non-stationary location)"
    n_params: int = 4
    covariate: np.ndarray | None = None
    _cov_mean: float = 0.0
    _cov_sd: float = 1.0

    @classmethod
    def fit(cls, x, covariate=None):
        x = np.asarray(x, float)
        t = np.arange(len(x), dtype=float) if covariate is None else np.asarray(covariate, float)
        tm, ts = float(t.mean()), float(t.std(ddof=0)) or 1.0
        s = (t - tm) / ts
        xx = -x

        def nllh(th):
            mu0, mu1, logsig, xi = th
            mu = mu0 + mu1 * s
            sig = np.exp(logsig)
            y = (xx - mu) / sig
            if abs(xi) < 1e-8:
                return len(xx) * logsig + np.sum(y) + np.sum(np.exp(-y))
            z = 1 + xi * y
            if np.any(z <= 1e-12):
                return 1e12
            return len(xx) * logsig + (1 + 1 / xi) * np.sum(np.log(z)) + np.sum(z ** (-1 / xi))

        m0, s0, x0, _ = _fit_gev_mle(xx)
        best, bf = None, np.inf
        for mu1_0 in (0.0, 0.1 * s0, -0.1 * s0):
            r = optimize.minimize(
                nllh, np.array([m0, mu1_0, np.log(s0), x0]), method="Nelder-Mead",
                options=dict(maxiter=20000, maxfev=20000),
            )
            if r.fun < bf:
                best, bf = r.x, r.fun
        obj = cls(
            params=dict(loc0=best[0], loc1=best[1], scale=np.exp(best[2]), shape=best[3]),
            data=x,
            covariate=t,
            _cov_mean=tm,
            _cov_sd=ts,
        )
        obj._nllh_cached = bf
        return obj

    def _mu(self, t):
        s = (np.asarray(t, float) - self._cov_mean) / self._cov_sd
        return self.params["loc0"] + self.params["loc1"] * s

    def nllh(self, x=None) -> float:
        return float(getattr(self, "_nllh_cached", np.inf))

    def frozen_at(self, t) -> GEVMinima:
        """The stationary GEV implied at covariate value ``t`` (an epoch slice)."""
        return GEVMinima(
            params=dict(loc=float(self._mu(t)), scale=self.params["scale"], shape=self.params["shape"]),
            data=self.data,
        )

    def return_level(self, T, t=None):
        t = self.covariate[-1] if t is None else t
        return self.frozen_at(t).return_level(T)

    def cdf(self, q):
        return self.frozen_at(self.covariate[-1]).cdf(q)

    def pdf(self, q):
        return self.frozen_at(self.covariate[-1]).pdf(q)

    def ppf(self, p):
        return self.frozen_at(self.covariate[-1]).ppf(p)


# ---------------------------------------------------------------------------
# lower-bounded families fitted in Q-space
# ---------------------------------------------------------------------------
@dataclass(repr=False)
class PearsonIII(Marginal):
    """Pearson type III fitted by L-moments (Hosking), i.e. R's ``parpe3``.

    ``scipy.stats.pearson3`` takes ``(skew, loc, scale)`` which is exactly
    lmomco's ``(gamma, mu, sigma)``, so the two agree parameter for parameter.
    Note the hard support bound at ``loc - 2*scale/skew`` for positive skew: an
    observation outside it makes ``cdf`` return exactly 0 or 1, which is what
    broke ``fitCopula`` in the R scripts (16.5% of 60-year samples). The bound
    is recorded in :attr:`notes` when it bites.
    """

    name: str = "Pearson III"
    n_params: int = 3

    @classmethod
    def fit(cls, x):
        x = np.asarray(x, float)
        if not _HAVE_LMOM:
            raise ImportError("lmoments3 is required for the L-moment Pearson III fit")
        p = _ld.pe3.lmom_fit(x)
        obj = cls(
            params=dict(skew=float(p["skew"]), loc=float(p["loc"]), scale=float(p["scale"])),
            data=x,
        )
        lo, hi = obj.support()
        if np.any(x < lo) or np.any(x > hi):
            obj.notes.append("sample value outside fitted support: cdf saturates at 0 or 1")
        if lo > 0:
            obj.notes.append(f"implied lower bound on Q = {lo:.4g}")
        return obj

    def _frozen(self):
        p = self.params
        return stats.pearson3(p["skew"], loc=p["loc"], scale=p["scale"])

    def support(self):
        p = self.params
        if abs(p["skew"]) < 1e-9:
            return -np.inf, np.inf
        b = p["loc"] - 2 * p["scale"] / p["skew"]
        return (b, np.inf) if p["skew"] > 0 else (-np.inf, b)

    def cdf(self, q):
        return self._frozen().cdf(q)

    def pdf(self, q):
        return self._frozen().pdf(q)

    def ppf(self, p):
        return self._frozen().ppf(p)


@dataclass(repr=False)
class GammaDist(Marginal):
    """Two-parameter gamma, MLE with the location fixed at zero."""

    name: str = "Gamma"
    n_params: int = 2

    @classmethod
    def fit(cls, x):
        x = np.asarray(x, float)
        a, loc, sc = stats.gamma.fit(x, floc=0)
        return cls(params=dict(shape=a, scale=sc), data=x)

    def _frozen(self):
        return stats.gamma(self.params["shape"], loc=0, scale=self.params["scale"])

    def cdf(self, q):
        return self._frozen().cdf(q)

    def pdf(self, q):
        return self._frozen().pdf(q)

    def ppf(self, p):
        return self._frozen().ppf(p)


@dataclass(repr=False)
class WeibullDist(Marginal):
    """Two-parameter Weibull, MLE with the location fixed at zero.

    The two-parameter Weibull is the extreme-value limit for minima of a
    bounded-below variable with a lower bound at zero, which is why it is a
    standard low-flow choice rather than merely a convenient one.
    """

    name: str = "Weibull"
    n_params: int = 2

    @classmethod
    def fit(cls, x):
        x = np.asarray(x, float)
        c, loc, sc = stats.weibull_min.fit(x, floc=0)
        return cls(params=dict(shape=c, scale=sc), data=x)

    def _frozen(self):
        return stats.weibull_min(self.params["shape"], loc=0, scale=self.params["scale"])

    def cdf(self, q):
        return self._frozen().cdf(q)

    def pdf(self, q):
        return self._frozen().pdf(q)

    def ppf(self, p):
        return self._frozen().ppf(p)


FAMILIES = {
    "GEV": GEVMinima,
    "Pearson III": PearsonIII,
    "Gamma": GammaDist,
    "Weibull": WeibullDist,
}


def fit_all(x, families=None) -> dict:
    """Fit every family, skipping (with a warning note) any that fails."""
    families = FAMILIES if families is None else families
    out = {}
    for nm, cls in families.items():
        try:
            out[nm] = cls.fit(x)
        except Exception as exc:  # keep going; report at the end
            out[nm] = exc
    return out


# ---------------------------------------------------------------------------
# goodness of fit
# ---------------------------------------------------------------------------
def _ad_statistic(fit, x) -> float:
    """Anderson-Darling A^2 against the fitted cdf."""
    z = np.sort(np.clip(np.asarray(fit.cdf(x), float), 1e-12, 1 - 1e-12))
    n = len(z)
    i = np.arange(1, n + 1)
    return float(-n - np.sum((2 * i - 1) * (np.log(z) + np.log1p(-z[::-1]))) / n)


def ad_test(fit, *, B: int = 499, rng=None) -> dict:
    """Anderson-Darling test with a **parametric bootstrap** null distribution.

    The null is generated by simulating from the fitted model and *refitting* on
    each simulated sample, which is what makes the p-value valid when the
    parameters were estimated from the data. Compare with
    :func:`ks_test_naive`, which does not do this and consequently never
    rejects.
    """
    rng = np.random.default_rng(rng)
    cls, x, n = type(fit), fit.data, fit.n
    a_obs = _ad_statistic(fit, x)
    null = []
    for _ in range(B):
        xb = fit.rvs(n, rng)
        try:
            fb = cls.fit(xb)
            null.append(_ad_statistic(fb, xb))
        except Exception:
            continue
    null = np.asarray(null, float)
    p = (1 + np.sum(null >= a_obs)) / (len(null) + 1) if null.size else np.nan
    return dict(
        distribution=fit.name,
        A2=a_obs,
        p_value=float(p),
        B_effective=int(null.size),
        A2_null_95=float(np.percentile(null, 95)) if null.size else np.nan,
    )


def ks_test_naive(fit) -> dict:
    """Kolmogorov-Smirnov against the fitted cdf, parameters estimated in-sample.

    Provided only for comparison with the R scripts. The p-value is **not
    valid** -- the null distribution of D is not the Kolmogorov one once the
    parameters have been estimated from the same sample. Simulation with data
    drawn from a true gamma and Pearson III fitted to it gives a rejection rate
    of 0.000 at alpha=0.05 and a median p-value of 0.92, so ranking candidate
    distributions by this number ranks noise. Use :func:`ad_test`.
    """
    d = stats.kstest(fit.data, fit.cdf)
    return dict(distribution=fit.name, KS=float(d.statistic), p_value_invalid=float(d.pvalue))


def compare(fits: dict, *, B: int = 499, rng=None, include_ks=True) -> pd.DataFrame:
    """AIC table with bootstrap Anderson-Darling, sorted by AIC."""
    rng = np.random.default_rng(rng)
    rows = []
    for nm, f in fits.items():
        if isinstance(f, Exception):
            rows.append(dict(distribution=nm, notes=f"FIT FAILED: {type(f).__name__}: {f}",
                             error=str(f)))
            continue
        r = dict(
            distribution=nm,
            k=f.n_params,
            nllh=f.nllh(),
            AIC=f.aic,
            BIC=f.bic,
            **{k: v for k, v in f.params.items()},
        )
        r.update({k: v for k, v in ad_test(f, B=B, rng=rng).items() if k != "distribution"})
        if include_ks:
            r["KS_p_invalid"] = ks_test_naive(f)["p_value_invalid"]
        r["notes"] = "; ".join(f.notes)
        rows.append(r)
    df = pd.DataFrame(rows)
    if "AIC" in df:
        df = df.sort_values("AIC").reset_index(drop=True)
        df.insert(df.columns.get_loc("AIC") + 1, "dAIC", df["AIC"] - df["AIC"].min())
    return df


def lr_test(fit_small, fit_big) -> dict:
    """Likelihood-ratio test of a nested pair (e.g. stationary vs trend GEV)."""
    d = 2 * (fit_small.nllh() - fit_big.nllh())
    df = fit_big.n_params - fit_small.n_params
    return dict(
        deviance=float(d),
        df=int(df),
        p_value=float(stats.chi2.sf(d, df)),
        verdict="trend supported" if stats.chi2.sf(d, df) < 0.05 else "trend not supported",
    )


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------
def lmom_ratios(x) -> dict:
    """Sample L-moments and L-moment ratios (L-CV, L-skew, L-kurtosis)."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    # Hosking's probability-weighted-moment estimator
    b = []
    for r in range(4):
        j = np.arange(1, n + 1)
        w = np.ones(n, float)
        for k in range(r):
            w *= (j - 1 - k) / (n - 1 - k)
        b.append(np.sum(w * x) / n)
    l1 = b[0]
    l2 = 2 * b[1] - b[0]
    l3 = 6 * b[2] - 6 * b[1] + b[0]
    l4 = 20 * b[3] - 30 * b[2] + 12 * b[1] - b[0]
    return dict(l1=l1, l2=l2, t=l2 / l1, t3=l3 / l2, t4=l4 / l2)


def theoretical_lmom_curve(family: str, t3=np.linspace(-0.4, 0.6, 121)):
    """(t3, t4) locus of a two/three-parameter family, for an L-moment ratio diagram."""
    t3 = np.asarray(t3, float)
    if family == "GEV":
        # Hosking & Wallis (1997) polynomial approximation
        c = [0.10701, 0.11090, 0.84838, -0.06669, 0.00567, -0.04208, 0.03763]
        return c[0] + c[1] * t3 + c[2] * t3**2 + c[3] * t3**3 + c[4] * t3**4 \
            + c[5] * t3**5 + c[6] * t3**6
    if family == "Pearson III":
        c = [0.1224, 0.30115, 0.95812, -0.57488, 0.19383]
        return c[0] + c[1] * t3**2 + c[2] * t3**4 + c[3] * t3**6 + c[4] * t3**8
    if family == "Generalised logistic":
        c = [0.16667, 0.0, 0.83333]
        return c[0] + c[2] * t3**2
    if family == "Lognormal":
        c = [0.12282, 0.77518, 0.12279, -0.13638, 0.11368]
        return c[0] + c[1] * t3**2 + c[2] * t3**4 + c[3] * t3**6 + c[4] * t3**8
    raise ValueError(family)


def trend_tests(values, years=None) -> pd.DataFrame:
    """Mann-Kendall trend, Theil-Sen slope, lag-1 autocorrelation, Ljung-Box.

    Replaces the R scripts' ``MannKendall`` + ``Box.test(lag = 10)`` pair. The
    Ljung-Box lag is reduced to a data-driven ``min(10, n//5)`` because lag 10
    on a 55-point series has almost no power.
    """
    v = np.asarray(values, float)
    n = len(v)
    t = np.arange(n, dtype=float) if years is None else np.asarray(years, float)
    mk = stats.kendalltau(t, v)
    ts = stats.theilslopes(v, t, 0.95)
    r1 = float(np.corrcoef(v[:-1], v[1:])[0, 1]) if n > 2 else np.nan
    from statsmodels.stats.diagnostic import acorr_ljungbox

    lag = max(1, min(10, n // 5))
    lb = acorr_ljungbox(v - v.mean(), lags=[lag], return_df=True)
    return pd.DataFrame(
        [
            dict(test="Mann-Kendall (tau)", statistic=mk.statistic, p_value=mk.pvalue,
                 detail="monotonic trend in block minima"),
            dict(test="Theil-Sen slope", statistic=ts[0], p_value=np.nan,
                 detail=f"units/yr, 95% CI [{ts[2]:.4g}, {ts[3]:.4g}]"),
            dict(test="lag-1 autocorrelation", statistic=r1, p_value=np.nan,
                 detail="serial dependence between consecutive block minima"),
            dict(test=f"Ljung-Box (lag {lag})", statistic=float(lb['lb_stat'].iloc[0]),
                 p_value=float(lb['lb_pvalue'].iloc[0]), detail="independence of block minima"),
        ]
    )
