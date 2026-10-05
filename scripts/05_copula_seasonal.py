"""Bivariate joint low-flow analysis, two gauges, SEASONAL block minima.

Port of ``copulaseas.R``. Identical to script 04 except for the block
definition -- which is the point: in the R pair, the seasonal copula script's
plotting block recomputed the observation cloud from the full-year data,
silently replacing the seasonal minima used for the fit. The curves were
seasonal and the scatter was annual, so the points sat systematically low
against the curves. Here the block definition lives in one place and is used
for the fit and the figure alike.

The season is chosen from the record, as in script 02: the window with the
highest capture rate, taken from the gauge named in ``SEASON_FROM``. Both
gauges must share one season definition -- a copula fitted to blocks defined
differently from its own marginals cannot be reconciled with anything.

Run:  python scripts/05_copula_seasonal.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lowflow import best_season, read_flow, scan_seasons
from lowflow.pipeline import BivariateConfig, run_bivariate

# Point these at your own workbooks; the two columns must be date and flow.
NAME_X, PATH_X = "Station A", "data/station_a.xlsx"
NAME_Y, PATH_Y = "Station B", "data/station_b.xlsx"

D = 7
MIN_COV = 0.90
SEASON = None          # None = choose from the record; or pin e.g. (7, 8, 9)
SEASON_FROM = "y"      # which gauge the scan is run on: "x" or "y"


def main():
    season = SEASON
    if season is None:
        path, name = (PATH_X, NAME_X) if SEASON_FROM == "x" else (PATH_Y, NAME_Y)
        rec = read_flow(path, name=name, D=D, zero_policy="keep")
        scan = scan_seasons(rec, min_coverage=MIN_COV)
        print(f"Season scan on {name}, best windows by capture rate:")
        print(scan.head(6).to_string(index=False))
        season = best_season(rec, min_coverage=MIN_COV)
    print(f"\nSeason in use: {season}   (keep this identical to script 02)")

    cfg = BivariateConfig(
        path_x=PATH_X,
        path_y=PATH_Y,
        name_x=NAME_X,
        name_y=NAME_Y,
        label=f"{NAME_X} and {NAME_Y}, seasonal minimum {D}-day mean, months {season}",
        outdir="out/05_copula_seasonal",
        season=season,
        D=D,
        min_coverage=MIN_COV,
        start_month=1,
        zero_policy="keep",
        T=(2, 5, 10, 20, 50, 100),
        B_gof=499,
        # Only copulas with upper-tail dependence on the droughtiness scale are
        # admissible here: Gaussian, Frank and plain Clayton all have
        # lambda_drought = 0. Add them back as a one-off sensitivity check.
        families=("Gumbel", "Joe", "Survival Clayton"),
    )
    res = run_bivariate(cfg)
    print(f"\nDone. {len(res['written'])} files in {res['outdir']}")


if __name__ == "__main__":
    main()
