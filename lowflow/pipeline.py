"""End-to-end pipelines, driven by configuration rather than by copy-paste.

The R original was five scripts that were near-copies of one another, and four
of the review findings were caused by drift between those copies: the
projection script still wrote into the observed-record output folder (A4), the
seasonal copula script recomputed its observation cloud with the wrong block
definition (A5), the normal-copula figure kept the Gumbel label table (B1), and
the ``D``-day switch was live in one copy, dead in another and absent from the
third (C3). Those are not statistics mistakes; they are duplication mistakes.

So there is exactly one univariate implementation and one bivariate
implementation here. The five scripts in ``scripts/`` differ only by a
configuration object, and each one writes to its own output directory because
the directory is part of that configuration.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import copulas as cop
from . import figures as figs
from .blocks import block_minima, scan_seasons, season_diagnostic, season_label
from .dataio import FlowRecord, read_flow
from .marginals import (
    GEVMinima,
    GEVMinimaNS,
    compare,
    fit_all,
    lr_test,
    trend_tests,
)

__all__ = ["UnivariateConfig", "BivariateConfig", "run_univariate", "run_bivariate"]

T_DEFAULT = (2, 5, 10, 20, 50, 100)


# ---------------------------------------------------------------------------
@dataclass
class UnivariateConfig:
    """Everything that distinguishes one univariate run from another."""

    path: str
    station: str
    outdir: str
    label: str                      # appears in figure titles and the log
    season: tuple | None = None     # None = annual blocks
    D: int = 7                      # 7-day mean (MAM7 / 7Q10 convention)
    min_coverage: float = 0.9
    start_month: int = 1
    zero_policy: str = "keep"
    unit: str = "m$^3$/s"
    T: tuple = T_DEFAULT
    B_gof: int = 499                # bootstrap replicates for Anderson-Darling
    B_ci: int = 999                 # bootstrap replicates for non-GEV CIs
    nonstationary: bool = False     # fit a trend in the GEV location and test it
    epochs: tuple = ()              # e.g. ((2011, 2040), (2041, 2070), (2071, 2100))
    scan_season: bool = True        # rank candidate low-flow windows
    seed: int = 20240617
    date_col: str = "date"
    flow_col: str = "flow"


@dataclass
class BivariateConfig:
    """Configuration for a two-river joint low-flow analysis."""

    path_x: str
    path_y: str
    name_x: str
    name_y: str
    outdir: str
    label: str
    season: tuple | None = None
    D: int = 7
    min_coverage: float = 0.9
    start_month: int = 1
    zero_policy: str = "keep"
    unit: str = "m$^3$/s"
    T: tuple = T_DEFAULT
    B_gof: int = 499
    B_ci: int = 999
    families: tuple = ("Gumbel", "Joe", "Survival Clayton")
    seed: int = 20240617
    date_col: str = "date"
    flow_col: str = "flow"


class _Log:
    """Collects a human-readable run log and writes it next to the outputs."""

    def __init__(self):
        self.lines: list[str] = []

    def __call__(self, *parts):
        s = " ".join(str(p) for p in parts)
        self.lines.append(s)
        print(s)

    def rule(self, title=""):
        self("\n" + "=" * 78)
        if title:
            self(title)
            self("=" * 78)

    def table(self, df: pd.DataFrame, floatfmt="%.5g"):
        with pd.option_context("display.width", 200, "display.max_columns", 50,
                               "display.float_format", lambda v: floatfmt % v):
            self(df.to_string(index=False))

    def write(self, path):
        Path(path).write_text("\n".join(self.lines), encoding="utf-8")


def _save(fig, outdir: Path, stem: str, written: list):
    for ext in ("png", "pdf"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p)
        written.append(str(p))
    plt.close(fig)


def _csv(df: pd.DataFrame, outdir: Path, stem: str, written: list):
    p = outdir / f"{stem}.csv"
    df.to_csv(p, index=False)
    written.append(str(p))


# ---------------------------------------------------------------------------
def run_univariate(cfg: UnivariateConfig) -> dict:
    """Marginal low-flow frequency analysis for one station."""
    rng = np.random.default_rng(cfg.seed)
    out = Path(cfg.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figs.apply_style()
    log, written = _Log(), []
    log.rule(f"UNIVARIATE LOW-FLOW ANALYSIS -- {cfg.label}")

    # 1. data ---------------------------------------------------------------
    rec = read_flow(cfg.path, name=cfg.station, D=cfg.D, zero_policy=cfg.zero_policy,
                    date_col=cfg.date_col, flow_col=cfg.flow_col)
    log("\n[1] Record screening")
    log(rec.summary())

    # 2. is the season the right one? --------------------------------------
    scan = None
    if cfg.scan_season:
        scan = scan_seasons(rec, min_coverage=cfg.min_coverage, start_month=cfg.start_month)
        log("\n[2] Candidate low-flow windows, ranked by capture of the annual minimum")
        log.table(scan.head(8).drop(columns=["months"]))
        _csv(scan, out, "season_scan", written)
        if cfg.season:
            chosen = season_label(cfg.season)
            row = scan.loc[scan.window == chosen]
            if len(row):
                cr = float(row.capture_rate.iloc[0])
                log(f"\n    chosen window {chosen}: capture rate {cr:.1%}")
                if cr < 0.8:
                    log(f"    WARNING: {chosen} contains the annual minimum in only {cr:.0%} of "
                        f"years. The best window is {scan.window.iloc[0]} at "
                        f"{scan.capture_rate.iloc[0]:.0%}. A seasonal analysis on {chosen} is "
                        "not describing the annual drought.")
            _save(figs.fig_season_scan(scan, chosen=chosen,
                                       title=f"{cfg.station}: does the window contain the drought?"),
                  out, "fig_season_scan", written)

    if cfg.season:
        diag = season_diagnostic(rec, cfg.season, min_coverage=cfg.min_coverage,
                                 start_month=cfg.start_month)
        _csv(diag, out, "season_vs_annual_minima", written)
        log(f"    annual minimum falls inside the window in "
            f"{diag['in_season'].mean():.1%} of blocks "
            f"({int(diag['in_season'].sum())}/{len(diag)})")

    # 3. block minima -------------------------------------------------------
    bm = block_minima(rec, season=cfg.season, min_coverage=cfg.min_coverage,
                      start_month=cfg.start_month)
    log("\n[3] Block minima")
    log(bm.summary())
    _csv(bm.table, out, "block_minima", written)
    x = bm.values
    if len(x) < 15:
        log(f"    WARNING: only {len(x)} usable blocks; three-parameter fits are unstable.")

    # 4. stationarity and independence -------------------------------------
    tr = trend_tests(x, bm.years)
    log("\n[4] Trend and independence of the block minima")
    log.table(tr)
    _csv(tr, out, "trend_tests", written)
    _save(figs.fig_block_minima(bm, tr, unit=cfg.unit), out, "fig_block_minima", written)

    # 5. marginal fits ------------------------------------------------------
    fits = fit_all(x)
    log("\n[5] Candidate distributions, ranked by AIC")
    cmp = compare(fits, B=cfg.B_gof, rng=rng)
    log.table(cmp)
    _csv(cmp, out, "distribution_comparison", written)
    log("\n    The KS column is reported only to show why it cannot be used for selection:")
    log("    with parameters estimated in-sample its rejection rate at alpha=0.05 is 0.000")
    log("    and its median p-value 0.92. Selection is by AIC; the test is bootstrap A^2.")

    ok = {k: v for k, v in fits.items() if not isinstance(v, Exception)}
    for k, v in fits.items():
        if isinstance(v, Exception):
            log(f"    {k} FAILED: {v}")
    best = cmp.loc[0, "distribution"] if "distribution" in cmp else None
    log(f"\n    selected by AIC: {best}")
    for nm, f in ok.items():
        if f.notes:
            log(f"    note ({nm}): {'; '.join(f.notes)}")

    _save(figs.fig_density_fits(x, ok, unit=cfg.unit, title=f"{cfg.label}: fitted densities"),
          out, "fig_density_fits", written)
    _save(figs.fig_qq(x, ok, unit=cfg.unit, title=f"{cfg.label}: Q-Q"), out, "fig_qq", written)
    _save(figs.fig_lmom_diagram(x, title=f"{cfg.label}: L-moment ratio diagram"),
          out, "fig_lmom_diagram", written)

    # 6. return levels with uncertainty ------------------------------------
    log("\n[6] Low-flow return levels")
    rl_all = pd.DataFrame({"T": list(cfg.T)})
    for nm, f in ok.items():
        rl_all[nm] = f.return_level(np.array(cfg.T, float))
    log.table(rl_all)
    _csv(rl_all, out, "return_levels_all_distributions", written)

    ci = None
    if best in ok:
        bf = ok[best]
        ci = bf.return_level_ci(np.array(cfg.T, float), B=cfg.B_ci, rng=rng)
        log(f"\n    95% confidence interval for the selected model ({best}):")
        log.table(ci)
        _csv(ci, out, "return_levels_ci_selected", written)
        if "lower" in ci and (ci["lower"] < 0).any():
            log("    WARNING: the interval extends below zero discharge at long return periods.")
            log("    The fit is extrapolating past physical plausibility; prefer a "
                "lower-bounded family or fit in log space before quoting those values.")
        if isinstance(bf, GEVMinima):
            log(f"    GEV shape = {bf.params['shape']:.4f}; implied lower bound on Q = "
                f"{bf.lower_bound:.4g}")
        # denser grid so the band reads as a band rather than a few blocks
        ci_band = bf.return_level_ci(
            np.logspace(np.log10(min(cfg.T)), np.log10(max(cfg.T)), 14),
            B=cfg.B_ci, rng=np.random.default_rng(cfg.seed + 1))
    else:
        ci_band = None
    _save(figs.fig_return_levels(ok, ci_band, x, best=best, unit=cfg.unit,
                                 title=f"{cfg.label}: low-flow return levels", T=cfg.T),
          out, "fig_return_levels", written)

    # 7. non-stationarity ---------------------------------------------------
    ns_res = None
    if cfg.nonstationary:
        log("\n[7] Non-stationary GEV (linear trend in the location parameter)")
        st = ok.get("GEV") or GEVMinima.fit(x)
        ns = GEVMinimaNS.fit(x, covariate=bm.years.astype(float))
        test = lr_test(st, ns)
        log(f"    stationary   nllh = {st.nllh():.4f}  AIC = {st.aic:.2f}")
        log(f"    with trend   nllh = {ns.nllh():.4f}  AIC = {ns.aic:.2f}")
        log(f"    loc1 = {ns.params['loc1']:+.4f} per SD of year "
            f"({'drying' if ns.params['loc1'] > 0 else 'wetting'})")
        log(f"    likelihood-ratio deviance = {test['deviance']:.3f} on {test['df']} df, "
            f"p = {test['p_value']:.4g} -> {test['verdict']}")
        rows = []
        for yr in (bm.years.min(), int(np.median(bm.years)), bm.years.max()):
            rl = ns.return_level(np.array(cfg.T, float), t=float(yr))
            rows.append(pd.DataFrame({"year": yr, "T": list(cfg.T), "return_level": rl}))
        ns_tab = pd.concat(rows, ignore_index=True)
        log("\n    effective return levels by epoch:")
        log.table(ns_tab.pivot(index="T", columns="year", values="return_level").reset_index())
        _csv(ns_tab, out, "return_levels_nonstationary", written)
        ns_res = dict(test=test, params=ns.params, table=ns_tab)

    if cfg.epochs:
        log("\n[7b] Fixed-epoch stationary fits (an alternative to the trend model)")
        rows = []
        for lo, hi in cfg.epochs:
            sel = (bm.years >= lo) & (bm.years <= hi)
            if sel.sum() < 15:
                log(f"    {lo}-{hi}: only {int(sel.sum())} blocks, skipped")
                continue
            ef = fit_all(x[sel])
            ec = compare(ef, B=0 if cfg.B_gof == 0 else 199, rng=rng)
            ebest = ec.loc[0, "distribution"]
            rl = ef[ebest].return_level(np.array(cfg.T, float))
            rows.append(pd.DataFrame({"epoch": f"{lo}-{hi}", "n": int(sel.sum()),
                                      "selected": ebest, "T": list(cfg.T), "return_level": rl}))
        if rows:
            ep = pd.concat(rows, ignore_index=True)
            log.table(ep.pivot(index="T", columns="epoch", values="return_level").reset_index())
            _csv(ep, out, "return_levels_by_epoch", written)

    log.rule("OUTPUTS")
    for w in written:
        log("   ", w)
    log.write(out / "run_log.txt")
    return dict(record=rec, blocks=bm, fits=fits, comparison=cmp, return_levels=rl_all,
                ci=ci, trend=tr, season_scan=scan, nonstationary=ns_res,
                outdir=str(out), written=written, log="\n".join(log.lines))


# ---------------------------------------------------------------------------
def _joint_records(ra: FlowRecord, rb: FlowRecord) -> tuple[FlowRecord, FlowRecord]:
    """Restrict two records to the days on which BOTH carry a valid D-day mean.

    The R scripts did this implicitly with ``inner_join`` on date and never
    reported it. Doing it explicitly means the per-block coverage screen counts
    *common* days, so a year in which the two gauges have disjoint gaps is
    rejected rather than contributing a spurious joint minimum.
    """
    ia, ib = pd.DatetimeIndex(ra.daily.index), pd.DatetimeIndex(rb.daily.index)
    if (ia.tz is None) != (ib.tz is None):
        raise ValueError(
            f"{ra.name} and {rb.name} disagree on timezone awareness "
            f"({ia.dtype} vs {ib.dtype}). Their indexes can never align, so every "
            "joint block would be empty. read_flow normalises this; if you built "
            "these records by hand, strip the timezone from both first."
        )
    if ia.tz is not None:
        ia, ib = ia.tz_localize(None), ib.tz_localize(None)
    a = ra.daily.set_axis(ia).reindex(ia.union(ib))
    b = rb.daily.set_axis(ib).reindex(ia.union(ib))
    full = a.index
    both = a["flow_D"].notna() & b["flow_D"].notna()
    outs = []
    for rec, d in ((ra, a), (rb, b)):
        d = d.copy()
        # Mask, do NOT drop: the block still spans the whole calendar period, so
        # n_possible stays the calendar day count and the coverage screen in
        # block_minima remains meaningful. Reindexing onto the common days alone
        # would make coverage identically 1.0 and silently disable the screen.
        d.loc[~both, ["flow", "flow_D"]] = np.nan
        new = FlowRecord(daily=d, name=rec.name, D=rec.D, diagnostics=dict(rec.diagnostics))
        new.diagnostics["n_common_days"] = int(both.sum())
        outs.append(new)
    return outs[0], outs[1]


def run_bivariate(cfg: BivariateConfig) -> dict:
    """Joint low-flow frequency analysis for two rivers."""
    rng = np.random.default_rng(cfg.seed)
    out = Path(cfg.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figs.apply_style()
    log, written = _Log(), []
    log.rule(f"BIVARIATE JOINT LOW-FLOW ANALYSIS -- {cfg.label}")

    # 1. data, restricted to common days -----------------------------------
    ra = read_flow(cfg.path_x, name=cfg.name_x, D=cfg.D, zero_policy=cfg.zero_policy,
                   date_col=cfg.date_col, flow_col=cfg.flow_col)
    rb = read_flow(cfg.path_y, name=cfg.name_y, D=cfg.D, zero_policy=cfg.zero_policy,
                   date_col=cfg.date_col, flow_col=cfg.flow_col)
    log("\n[1] Records")
    log(ra.summary())
    log("")
    log(rb.summary())
    ja, jb = _joint_records(ra, rb)
    log(f"\n    days with a valid {cfg.D}-day mean at BOTH gauges: "
        f"{ja.diagnostics['n_common_days']}")

    bma = block_minima(ja, season=cfg.season, min_coverage=cfg.min_coverage,
                       start_month=cfg.start_month)
    bmb = block_minima(jb, season=cfg.season, min_coverage=cfg.min_coverage,
                       start_month=cfg.start_month)
    keep = sorted(set(bma.used.year) & set(bmb.used.year))
    ua = bma.used.set_index("year").loc[keep]
    ub = bmb.used.set_index("year").loc[keep]
    log("\n[2] Block minima on the common record")
    log(bma.summary())
    log("")
    log(bmb.summary())
    if not keep:
        raise ValueError(
            "No block is usable for BOTH rivers, so there is nothing to fit.\n"
            f"  common days with a valid {cfg.D}-day mean : "
            f"{ja.diagnostics['n_common_days']}\n"
            f"  {cfg.name_x}: {len(bma.used)} of {len(bma.table)} blocks passed coverage\n"
            f"  {cfg.name_y}: {len(bmb.used)} of {len(bmb.table)} blocks passed coverage\n"
            f"  overlapping years: {sorted(set(bma.used.year) & set(bmb.used.year))}\n"
            "Usual causes, in order of likelihood: (a) the two records share no "
            "calendar overlap; (b) min_coverage is too strict for a gappy record; "
            "(c) the chosen season selects too few days per block."
        )
    log(f"\n    blocks usable for BOTH rivers: {len(keep)} ({min(keep)}-{max(keep)})")
    if len(keep) < 20:
        log(f"    WARNING: {len(keep)} paired blocks is a thin sample. Three-parameter "
            "marginals and long return periods are very poorly determined here; treat "
            "T > ~2x the record length as indicative only.")
    pair = pd.DataFrame(
        {
            "year": keep,
            f"min_{cfg.name_x}": ua["min_flow"].to_numpy(),
            f"date_{cfg.name_x}": ua["date"].to_numpy(),
            f"min_{cfg.name_y}": ub["min_flow"].to_numpy(),
            f"date_{cfg.name_y}": ub["date"].to_numpy(),
        }
    )
    pair["lag_days"] = (pd.to_datetime(pair[f"date_{cfg.name_y}"])
                        - pd.to_datetime(pair[f"date_{cfg.name_x}"])).dt.days
    _csv(pair, out, "paired_block_minima", written)
    log(f"    median |lag| between the two rivers' minima: "
        f"{pair['lag_days'].abs().median():.0f} days")
    qx = pair[f"min_{cfg.name_x}"].to_numpy(float)
    qy = pair[f"min_{cfg.name_y}"].to_numpy(float)

    # 3. marginals ----------------------------------------------------------
    log("\n[3] Marginal distributions (both rivers, ranked by AIC)")
    marg, best_marg = {}, {}
    for nm, q in ((cfg.name_x, qx), (cfg.name_y, qy)):
        f = fit_all(q)
        c = compare(f, B=cfg.B_gof, rng=rng)
        log(f"\n  {nm}:")
        log.table(c)
        _csv(c, out, f"marginal_comparison_{nm}", written)
        b = c.loc[0, "distribution"]
        marg[nm], best_marg[nm] = f, b
        log(f"    selected: {b}")
        for k, v in f.items():
            if not isinstance(v, Exception) and v.notes:
                log(f"    note ({k}): {'; '.join(v.notes)}")
    mx = marg[cfg.name_x][best_marg[cfg.name_x]]
    my = marg[cfg.name_y][best_marg[cfg.name_y]]
    log(f"\n    NOTE: the R scripts fixed Pearson III for both rivers after testing")
    log(f"    only one of them. Here each river's marginal is selected on its own evidence:")
    log(f"    {cfg.name_x} -> {best_marg[cfg.name_x]}, {cfg.name_y} -> {best_marg[cfg.name_y]}.")

    # 4. dependence ---------------------------------------------------------
    u = cop.pseudo_obs(-qx, -qy)   # large u = deep drought
    dep = cop.dependence_summary(u)
    log("\n[4] Dependence of the two drought series (model-free)")
    log.table(dep)
    _csv(dep, out, "dependence_summary", written)
    _save(figs.fig_pseudo_obs(u, names=(cfg.name_x, cfg.name_y), dep=dep,
                              title=f"{cfg.label}: droughtiness pseudo-observations"),
          out, "fig_pseudo_obs", written)

    # 5. copulas ------------------------------------------------------------
    fits = cop.fit_all_copulas(u, families=cfg.families)
    cmp = cop.compare_copulas(fits, B=cfg.B_gof, rng=rng)
    log("\n[5] Copulas, ranked by AIC (Sn = Cramer-von Mises, parametric bootstrap)")
    log.table(cmp)
    _csv(cmp, out, "copula_comparison", written)
    best = cmp.loc[0, "copula"]
    model = fits[best]
    log(f"\n    selected by AIC: {best} (theta = {model.theta:.4f}, tau = {model.tau:.3f}, "
        f"lambda_drought = {model.lambda_drought:.3f})")
    lam_emp = float(dep.loc[dep.measure == "lambda_drought (CFG)", "value"].iloc[0])
    log(f"    non-parametric lambda_drought = {lam_emp:.3f} for comparison")
    if model.lambda_drought < 0.05 and lam_emp > 0.25:
        log(f"    CHECK: the selected copula implies lambda_drought = "
            f"{model.lambda_drought:.2f} (asymptotic independence in joint drought) while the")
        log(f"    non-parametric estimate is {lam_emp:.2f}. At n = {len(u)} the CFG estimator is")
        log("    biased upward, so this is not proof of a wrong choice, but the joint design")
        log("    values at long return periods depend on it. Compare the runner-up figure, and")
        log("    prefer a tail-dependent family if the two disagree materially.")
    log("    Families with lambda_drought = 0 (Gaussian, Frank, plain Clayton on these")
    log("    variables) assert that simultaneous extreme drought becomes independent in")
    log("    the limit. That is a modelling assertion, not a neutral default.")
    _save(figs.fig_copula_density(model, title=f"{cfg.label}: {best} copula"),
          out, "fig_copula_density", written)

    # 6. joint return periods ----------------------------------------------
    log("\n[6] Joint return levels")
    design = cop.and_design_point(model, np.array(cfg.T, float), mx, my)
    design = design.rename(columns={"x": cfg.name_x, "y": cfg.name_y})
    log("\n  AND event: both rivers simultaneously below their thresholds,")
    log("  at the equal-droughtiness design point u = v:")
    log.table(design)
    _csv(design, out, "and_design_points", written)

    mld = [cop.most_likely_design_point(model, T, mx, my) for T in cfg.T]
    mld = [m for m in mld if m]
    if mld:
        m = pd.DataFrame(mld).rename(columns={"x": cfg.name_x, "y": cfg.name_y})
        log("\n  Most-likely realisation on each AND curve:")
        log.table(m.drop(columns=["density"]))
        _csv(m, out, "and_most_likely_realisations", written)

    kc = cop.kendall_critical_level(model, np.array(cfg.T, float), rng=rng)
    log("\n  Kendall (survival) critical levels -- the whole region at least as rare")
    log("  as the critical layer, which is usually what a design criterion means:")
    log.table(kc)
    _csv(kc, out, "kendall_critical_levels", written)

    or_tab = pd.DataFrame({"T": list(cfg.T)})
    or_u = []
    for T in cfg.T:
        from scipy import optimize as _o
        try:
            w = _o.brentq(lambda z: float(cop.or_probability(model, z, z)) - 1.0 / T,
                          1e-10, 1 - 1e-10)
        except ValueError:
            w = np.nan
        or_u.append(w)
    or_tab["u"] = or_u
    or_tab[cfg.name_x] = [float(cop.discharge_from_u(mx, w)) if w == w else np.nan for w in or_u]
    or_tab[cfg.name_y] = [float(cop.discharge_from_u(my, w)) if w == w else np.nan for w in or_u]
    log("\n  OR event: at least one river below its threshold (u = v):")
    log.table(or_tab)
    _csv(or_tab, out, "or_design_points", written)
    log("\n    AND, OR and Kendall are three different definitions of a 'T-year joint")
    log("    drought' and give different discharges. Whichever is quoted must be named;")
    log("    the R scripts reported AND only, and labelled the figure with the")
    log("    opposite inequality.")

    curves = pd.concat([cop.and_curve(model, T, mx, my, n=400) for T in cfg.T],
                       ignore_index=True)
    _csv(curves, out, "and_curves", written)
    _save(
        figs.fig_and_curves(
            curves, design.rename(columns={cfg.name_x: "x", cfg.name_y: "y"}), qx, qy,
            names=(cfg.name_x, cfg.name_y), unit=cfg.unit,
            mld=[{**m, "x": m[cfg.name_x] if cfg.name_x in m else m["x"],
                  "y": m[cfg.name_y] if cfg.name_y in m else m["y"]} for m in mld] or None,
            title=f"{cfg.label}: joint AND return levels ({best} copula)",
            subtitle=r"$P(Q_{%s} < q_{%s}\ \mathrm{and}\ Q_{%s} < q_{%s}) = 1/T$"
                     % (cfg.name_x, cfg.name_x, cfg.name_y, cfg.name_y),
        ),
        out, "fig_and_curves", written,
    )

    # also draw the runner-up copula so the sensitivity is visible
    second = cmp.loc[1, "copula"] if len(cmp) > 1 else None
    if second and not isinstance(fits[second], Exception):
        m2 = fits[second]
        c2 = pd.concat([cop.and_curve(m2, T, mx, my, n=400) for T in cfg.T], ignore_index=True)
        d2 = cop.and_design_point(m2, np.array(cfg.T, float), mx, my)
        _save(
            figs.fig_and_curves(
                c2, d2, qx, qy, names=(cfg.name_x, cfg.name_y), unit=cfg.unit,
                title=f"{cfg.label}: joint AND return levels ({second} copula)",
                subtitle="runner-up by AIC, shown for sensitivity",
            ),
            out, "fig_and_curves_runner_up", written,
        )
        comp = design[["T", cfg.name_x, cfg.name_y]].merge(
            d2.rename(columns={"x": cfg.name_x, "y": cfg.name_y})[["T", cfg.name_x, cfg.name_y]],
            on="T", suffixes=(f"_{best}", f"_{second}"))
        log(f"\n  Sensitivity of the AND design point to the copula choice "
            f"({best} vs {second}):")
        log.table(comp)
        _csv(comp, out, "and_design_point_copula_sensitivity", written)

    (out / "config.json").write_text(json.dumps(asdict(cfg), indent=2, default=str),
                                     encoding="utf-8")
    written.append(str(out / "config.json"))
    log.rule("OUTPUTS")
    for w in written:
        log("   ", w)
    log.write(out / "run_log.txt")
    return dict(pairs=pair, marginals=marg, best_marginals=best_marg, u=u, dependence=dep,
                copulas=fits, copula_comparison=cmp, design=design, curves=curves,
                kendall=kc, outdir=str(out), written=written, log="\n".join(log.lines))
