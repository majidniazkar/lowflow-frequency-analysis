"""Publication-grade figures for low-flow frequency analysis.

Self-contained: no dependency on a notebook style plugin, so the scripts
reproduce identical output from a cold start. Every function takes fitted
objects and returns a ``matplotlib`` Figure; nothing is written to disk here.

Conventions applied throughout: a three-step role-mapped font ladder, outward
ticks, no top/right spines, frameless legends or direct end-of-line labels,
and 300 dpi on save. Return-period axes are logarithmic with human-readable
ticks. Rejected blocks are drawn as open markers and never enter a fitted
summary.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .marginals import lmom_ratios, theoretical_lmom_curve

__all__ = [
    "apply_style",
    "PALETTE",
    "fig_block_minima",
    "fig_density_fits",
    "fig_qq",
    "fig_return_levels",
    "fig_lmom_diagram",
    "fig_season_scan",
    "fig_pseudo_obs",
    "fig_copula_density",
    "fig_and_curves",
    "plotting_positions",
]

# distribution -> colour, threaded across every figure in the set (§4.1)
PALETTE = {
    "GEV": "#5B2C8D",
    "Pearson III": "#1B7F4B",
    "Gamma": "#C2410C",
    "Weibull": "#1F5FA8",
    "empirical": "#111111",
    "band": "#9A7FC0",
}
COPULA_PALETTE = {
    "Gumbel": "#5B2C8D",
    "Joe": "#8E44AD",
    "Survival Clayton": "#1B7F4B",
    "Gaussian": "#C2410C",
    "Frank": "#1F5FA8",
    "Clayton": "#7A7A7A",
}
GREY = "#6E6E6E"


def apply_style(base: int = 9) -> None:
    """Role-mapped font ladder and clean axes. Call once before plotting."""
    s1, s2, s3 = base, base - 1, base - 2
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": s1,
            "axes.titlesize": s1,
            "axes.labelsize": s1,
            "legend.fontsize": s2,
            "xtick.labelsize": s3,
            "ytick.labelsize": s3,
            "axes.titlelocation": "left",
            "axes.titleweight": "regular",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.frameon": False,
            "lines.linewidth": 1.6,
            "figure.constrained_layout.use": True,
        }
    )


def plotting_positions(x, a: float = 0.0) -> pd.DataFrame:
    """Weibull (a=0) or Gringorten (a=0.44) positions for block **minima**.

    The i-th smallest value has non-exceedance probability ``p_i`` and therefore
    low-flow return period ``T_i = 1 / p_i``.
    """
    x = np.sort(np.asarray(x, float))
    n = len(x)
    i = np.arange(1, n + 1)
    p = (i - a) / (n + 1 - 2 * a)
    return pd.DataFrame(dict(rank=i, value=x, p_nonexceed=p, T=1.0 / p))


def _finish(ax, *, xlabel=None, ylabel=None, title=None):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.margins(0.04)
    return ax


# ---------------------------------------------------------------------------
# univariate
# ---------------------------------------------------------------------------
def fig_block_minima(bm, trend: pd.DataFrame | None = None, *, unit="m$^3$/s"):
    """Block-minimum series, with rejected blocks shown as open markers."""
    t = bm.table
    used, rej = t[t.used], t[~t.used]
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    ax.plot(used.year, used.min_flow, "-o", color=PALETTE["empirical"], ms=3.4, lw=1.1,
            label=f"used (n = {len(used)})")
    if len(rej):
        ax.plot(rej.year, rej.min_flow, "o", mfc="none", mec="#B00020", ms=5, mew=1.1,
                label=f"rejected, coverage < threshold (n = {len(rej)})")
    if trend is not None:
        sl = trend.loc[trend.test == "Theil-Sen slope", "statistic"]
        mk = trend.loc[trend.test == "Mann-Kendall (tau)"]
        if len(sl):
            slope = float(sl.iloc[0])
            yr = used.year.to_numpy(float)
            med_y, med_x = np.median(used.min_flow), np.median(yr)
            ax.plot(yr, med_y + slope * (yr - med_x), "--", color="#C2410C", lw=1.3,
                    label=f"Theil-Sen {slope:+.3g} {unit}/yr"
                          f" (Mann-Kendall p = {float(mk.p_value.iloc[0]):.3f})")
    ax.legend(loc="upper right", ncols=1, handlelength=1.6)
    return _finish(fig.axes[0], xlabel="Year", ylabel=f"Minimum discharge ({unit})",
                   title=bm.label).figure


def fig_density_fits(x, fits: dict, *, unit="m$^3$/s", title=None):
    """Histogram of block minima with each fitted density overlaid."""
    x = np.asarray(x, float)
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.hist(x, bins="fd", density=True, color="#DCDCDC", edgecolor="white", lw=0.6,
            label=f"observed (n = {len(x)})")
    lo = max(0.0, x.min() - 0.25 * np.ptp(x))
    grid = np.linspace(lo, x.max() + 0.25 * np.ptp(x), 600)
    for nm, f in fits.items():
        if isinstance(f, Exception):
            continue
        with np.errstate(all="ignore"):
            y = f.pdf(grid)
        ax.plot(grid, np.where(np.isfinite(y), y, np.nan), color=PALETTE.get(nm, GREY),
                label=nm, lw=1.5)
    ax.legend(loc="upper right", handlelength=1.5)
    ax.set_ylim(bottom=0)
    return _finish(ax, xlabel=f"Block-minimum discharge ({unit})", ylabel="Probability density",
                   title=title or "Fitted distributions").figure


def fig_qq(x, fits: dict, *, unit="m$^3$/s", title=None):
    """Q-Q panel, one per family, all in discharge space.

    The sign handling that made the R version correct only by accident is gone:
    ``ppf`` is always a discharge, so the theoretical quantile for the i-th
    smallest observation is just ``ppf(p_i)``.
    """
    x = np.sort(np.asarray(x, float))
    pp = plotting_positions(x)
    items = [(nm, f) for nm, f in fits.items() if not isinstance(f, Exception)]
    ncol = 2
    nrow = int(np.ceil(len(items) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6, 2.7 * nrow), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, (nm, f) in zip(axes, items):
        q = f.ppf(pp.p_nonexceed.to_numpy())
        ax.plot(q, pp.value, "o", ms=3.2, color=PALETTE.get(nm, GREY), mec="white", mew=0.3)
        lim = [min(np.nanmin(q), x.min()), max(np.nanmax(q), x.max())]
        ax.plot(lim, lim, "-", color=GREY, lw=0.9)
        ax.set_title(nm)
        ax.margins(0.05)
    for ax in axes[len(items):]:
        ax.set_visible(False)
    fig.supxlabel(f"Theoretical quantile ({unit})", fontsize=mpl.rcParams["axes.labelsize"])
    fig.supylabel(f"Observed quantile ({unit})", fontsize=mpl.rcParams["axes.labelsize"])
    if title:
        fig.suptitle(title, x=0.01, ha="left", fontsize=mpl.rcParams["axes.titlesize"])
    return fig


def fig_return_levels(fits: dict, ci: pd.DataFrame | None, x, *, best: str | None = None,
                      unit="m$^3$/s", title=None, T=(2, 5, 10, 20, 50, 100)):
    """Low-flow return-level curves with a confidence band for the selected model.

    The band is the point of this figure. With ~50 block minima the 100-year low
    flow is estimated to within roughly a factor of two, and a figure that shows
    only the point estimates invites the reader to over-read them.
    """
    Tgrid = np.logspace(np.log10(1.02), np.log10(max(T) * 1.3), 200)
    fig, ax = plt.subplots(figsize=(6.0, 3.6))

    if ci is not None and len(ci):
        c = ci.dropna(subset=["lower", "upper"]).sort_values("T")
        if len(c) > 1:
            ax.fill_between(c["T"], c["lower"], c["upper"], color=PALETTE["band"],
                            alpha=0.25, lw=0,
                            label=f"95% interval, {best} \u2014 {c['method'].iloc[0]}")
    for nm, f in fits.items():
        if isinstance(f, Exception):
            continue
        focal = nm == best
        with np.errstate(all="ignore"):
            y = f.return_level(Tgrid)
        ax.plot(Tgrid, y, color=PALETTE.get(nm, GREY), lw=2.0 if focal else 1.1,
                alpha=1.0 if focal else 0.65, label=nm + (" (selected)" if focal else ""),
                zorder=3 if focal else 2)
    pp = plotting_positions(x)
    ax.plot(pp["T"], pp["value"], "o", ms=3.6, color=PALETTE["empirical"], mec="white",
            mew=0.4, label=f"observed (n = {len(pp)})", zorder=4)

    ax.set_xscale("log")
    ax.set_xticks([2, 5, 10, 20, 50, 100])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlim(1.6, max(T) * 1.35)
    ax.axhline(0, color="#B00020", lw=0.8, ls=":", zorder=1)
    ax.text(0.015, 0.045, "lower = more severe drought", transform=ax.transAxes,
            fontsize=mpl.rcParams["legend.fontsize"], color=GREY)
    ax.legend(loc="upper right", handlelength=1.6)
    return _finish(ax, xlabel="Return period T (years)",
                   ylabel=f"T-year low flow ({unit})",
                   title=title or "Low-flow return levels").figure


def fig_lmom_diagram(x, *, title=None):
    """L-moment ratio diagram: which three-parameter family the sample supports.

    More discriminating than any single-sample EDF test at n ~ 50, and it costs
    nothing: the sample is one point, each candidate family is one curve.
    """
    lm = lmom_ratios(x)
    t3 = np.linspace(-0.2, 0.7, 200)
    fig, ax = plt.subplots(figsize=(4.4, 3.6))
    for fam, col in [("GEV", PALETTE["GEV"]), ("Pearson III", PALETTE["Pearson III"]),
                     ("Lognormal", "#0F766E"), ("Generalised logistic", "#B45309")]:
        ax.plot(t3, theoretical_lmom_curve(fam, t3), color=col, lw=1.4, label=fam)
    ax.plot(lm["t3"], lm["t4"], "*", ms=15, color=PALETTE["empirical"],
            label=f"sample (n = {len(x)})", zorder=5)
    ax.set_xlim(-0.05, 0.65)
    ax.set_ylim(-0.05, 0.45)
    ax.legend(loc="upper left", handlelength=1.5)
    return _finish(ax, xlabel="L-skewness $\\tau_3$", ylabel="L-kurtosis $\\tau_4$",
                   title=title or "L-moment ratio diagram").figure


def fig_season_scan(scan: pd.DataFrame, chosen=None, *, top: int = 12, title=None):
    """Candidate low-flow windows ranked by how often they capture the annual minimum."""
    d = scan.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.0, 0.26 * len(d) + 1.4))
    colors = [
        "#B00020" if (chosen and w == chosen) else ("#1B7F4B" if r >= 0.8 else GREY)
        for w, r in zip(d.window, d.capture_rate)
    ]
    ax.barh(d.window, d.capture_rate, color=colors, height=0.72)
    for w, r in zip(d.window, d.capture_rate):
        ax.text(r + 0.012, w, f"{r:.0%}", va="center",
                fontsize=mpl.rcParams["legend.fontsize"])
    ax.axvline(0.8, color="#111111", lw=0.9, ls="--")
    ax.text(0.805, -0.75, "80%", fontsize=mpl.rcParams["legend.fontsize"], color=GREY)
    ax.set_xlim(0, 1.09)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.get_xaxis().set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    return _finish(ax, xlabel="Annual minima falling inside the window",
                   title=title or "Does the season contain the drought?").figure


# ---------------------------------------------------------------------------
# bivariate
# ---------------------------------------------------------------------------
def fig_pseudo_obs(u, names=("X", "Y"), *, dep: pd.DataFrame | None = None, title=None):
    """Pseudo-observations of droughtiness, with the dependence measures on it."""
    fig, ax = plt.subplots(figsize=(3.9, 3.7))
    ax.plot(u[:, 0], u[:, 1], "o", ms=4.2, color=PALETTE["empirical"], mec="white", mew=0.4)
    ax.plot([0, 1], [0, 1], color=GREY, lw=0.8, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if dep is not None:
        tau = dep.loc[dep.measure == "Kendall tau"]
        lam = dep.loc[dep.measure == "lambda_drought (CFG)"]
        bits = []
        if len(tau):
            bits.append(f"$\\tau$ = {float(tau.value.iloc[0]):.2f}"
                        f" (p = {float(tau.p_value.iloc[0]):.3f})")
        if len(lam):
            bits.append(f"$\\lambda_{{drought}}$ = {float(lam.value.iloc[0]):.2f}")
        ax.text(0.04, 0.96, "\n".join(bits), transform=ax.transAxes, va="top",
                fontsize=mpl.rcParams["legend.fontsize"])
    ax.text(0.97, 0.03, "upper right = both rivers in drought", transform=ax.transAxes,
            ha="right", fontsize=mpl.rcParams["legend.fontsize"], color=GREY)
    return _finish(ax, xlabel=f"{names[0]} droughtiness $u$", ylabel=f"{names[1]} droughtiness $v$",
                   title=title or "Pseudo-observations").figure


def fig_copula_density(model, *, title=None):
    """Contours of the fitted copula density on the unit square."""
    g = np.linspace(0.01, 0.99, 220)
    U, V = np.meshgrid(g, g, indexing="ij")
    with np.errstate(all="ignore"):
        Z = model.pdf(U, V)
    Z = np.where(np.isfinite(Z), Z, np.nan)
    fig, ax = plt.subplots(figsize=(4.1, 3.7))
    lv = np.nanpercentile(Z, [50, 70, 82, 90, 95, 98, 99.5])
    cs = ax.contour(U, V, Z, levels=np.unique(lv), colors=COPULA_PALETTE.get(model.name, GREY),
                    linewidths=1.1)
    ax.clabel(cs, inline=True, fontsize=mpl.rcParams["xtick.labelsize"], fmt="%.1f")
    ax.plot(model.u[:, 0], model.u[:, 1], "o", ms=3.0, color=PALETTE["empirical"],
            mec="white", mew=0.3, alpha=0.85)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return _finish(ax, xlabel="$u$", ylabel="$v$",
                   title=title or f"{model.name} copula density "
                                  f"($\\theta$ = {model.theta:.2f}, "
                                  f"$\\lambda_{{drought}}$ = {model.lambda_drought:.2f})").figure


def fig_and_curves(curves: pd.DataFrame, design: pd.DataFrame, obs_x, obs_y,
                   *, names=("X", "Y"), unit="m$^3$/s", mld: list | None = None,
                   title=None, subtitle=None):
    """Joint AND return-level curves in discharge space.

    Each curve is the locus of threshold pairs whose *joint drought*
    probability -- both rivers simultaneously below their threshold -- equals
    1/T. Solved by root-finding, so the curves are exactly on the contour and
    monotone; the R version filtered a grid and zig-zagged.
    """
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    Ts = sorted(curves["T"].unique())
    cmap = plt.get_cmap("viridis")
    cols = {T: cmap(i / max(len(Ts) - 1, 1)) for i, T in enumerate(Ts)}
    for T in Ts:
        c = curves[curves["T"] == T].sort_values("u")
        ax.plot(c["x"], c["y"], color=cols[T], lw=1.7, zorder=3)
        # label at the knee (closest approach to u = v), where the six curves are
        # well separated; at the curve ends they all pile onto the asymptotes
        i = int(np.argmin(np.abs(c["u"].to_numpy() - c["v"].to_numpy())))
        ax.annotate(f"{int(T)} yr", (c["x"].iloc[i], c["y"].iloc[i]),
                    textcoords="offset points", xytext=(7, 7), color=cols[T],
                    fontsize=mpl.rcParams["legend.fontsize"], zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    ax.plot(obs_x, obs_y, "o", ms=4.0, color=PALETTE["empirical"], mec="white", mew=0.4,
            label=f"observed block minima (n = {len(obs_x)})", zorder=4)
    if design is not None and len(design):
        d = design.dropna(subset=["x", "y"])
        ax.plot(d["x"], d["y"], "s", ms=5.0, mfc="white", mec="#B00020", mew=1.2,
                label="equal-droughtiness design point ($u = v$)", zorder=6)
    if mld:
        m = pd.DataFrame(mld).dropna(subset=["x", "y"])
        ax.plot(m["x"], m["y"], "D", ms=4.6, mfc="#B00020", mec="white", mew=0.6,
                label="most-likely realisation on the curve", zorder=7)
    ax.legend(loc="lower right", handlelength=1.4)
    ax.text(0.015, 0.02, "toward the origin = both rivers drier", transform=ax.transAxes,
            va="bottom", fontsize=mpl.rcParams["legend.fontsize"], color=GREY)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ttl = title or "Joint low-flow return levels"
    if subtitle:
        ttl = f"{ttl}\n{subtitle}"
    return _finish(ax, xlabel=f"{names[0]} discharge ({unit})",
                   ylabel=f"{names[1]} discharge ({unit})", title=ttl).figure
