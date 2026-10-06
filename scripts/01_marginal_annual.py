"""Univariate low-flow frequency analysis, ANNUAL block minima.

Fits GEV / Pearson III / Gamma / Weibull to the annual minimum D-day mean
discharge, selects by AIC, tests fit by bootstrap Anderson-Darling, and
reports return levels with confidence intervals.

Annual blocks need no assumption about when the basin runs dry, which makes
this the right starting point; move to script 02 only if a seasonal window
captures the annual minimum often enough to be worth using.

Run from the repository root:  python scripts/01_marginal_annual.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lowflow.pipeline import UnivariateConfig, run_univariate

# --- configuration ---------------------------------------------------------
# D = 7 gives the 7-day mean behind the conventional MAM7 / 7Q10 low-flow
# index. D = 1 is the single-day minimum: noisier, and more sensitive to
# gauging error at low stage.
CFG = UnivariateConfig(
    # Point these at your own workbook; the two columns must be date and flow.
    path="data/station_a.xlsx",
    station="Station A",
    label="Station A, annual minimum 7-day mean",
    outdir="out/01_marginal_annual",
    season=None,          # annual blocks
    D=7,
    min_coverage=0.90,    # a block needs 90% of its days to be used
    start_month=1,        # calendar year; appropriate where minima fall in summer,
                          # far from the block boundary, so no drought is split
    zero_policy="keep",   # zeros are real low-flow observations, not errors
    T=(2, 5, 10, 20, 50, 100),
    B_gof=499,
    B_ci=999,
    scan_season=True,     # report which window would capture the annual minimum
)

if __name__ == "__main__":
    res = run_univariate(CFG)
    print(f"\nDone. {len(res['written'])} files in {res['outdir']}")
