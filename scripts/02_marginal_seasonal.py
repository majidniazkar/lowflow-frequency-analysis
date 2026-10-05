"""Univariate low-flow frequency analysis, SEASONAL block minima.

Port of ``Marginal_Seas.R``.

On the season choice
--------------------
The R script hard-coded a three-month window. It also built a table flagging
whether each year's annual minimum actually fell inside that window -- and then
never summarised it.

This script does not assume a season. It ranks every contiguous 3-to-5-month
window by **capture rate** -- the fraction of annual minima the window actually
contains -- and takes the top-ranked one via :func:`lowflow.best_season`. The
full ranking is printed and written to ``season_scan.csv``, so the choice is
visible and auditable. Set ``SEASON`` to a tuple of months only if you must
reproduce an earlier analysis; `best_season` raises if no window reaches a 0.8
capture rate, which means no contiguous season describes the drought on your
record and annual blocks (script 01) are the honest answer.

Run:  python scripts/02_marginal_seasonal.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lowflow import best_season, read_flow, scan_seasons
from lowflow.pipeline import UnivariateConfig, run_univariate

# Point these at your own workbook; the two columns must be date and flow.
PATH = "data/station_a.xlsx"
STATION = "Station A"

D = 7
MIN_COV = 0.90        # 90% of the days in the window, not of the year
SEASON = None         # None = choose from the record; or pin e.g. (7, 8, 9)


def main():
    rec = read_flow(PATH, name=STATION, D=D, zero_policy="keep")
    print(rec.summary())

    scan = scan_seasons(rec, min_coverage=MIN_COV)
    print("\nSeason scan, best windows by capture rate:")
    print(scan.head(6).to_string(index=False))

    season = SEASON or best_season(rec, min_coverage=MIN_COV)
    print(f"\nSeason in use: {season}")

    cfg = UnivariateConfig(
        path=PATH,
        station=STATION,
        label=f"{STATION}, seasonal minimum {D}-day mean, months {season}",
        outdir="out/02_marginal_seasonal",
        season=season,
        D=D,
        min_coverage=MIN_COV,
        start_month=1,
        zero_policy="keep",
        T=(2, 5, 10, 20, 50, 100),
        B_gof=499,
        B_ci=999,
        scan_season=True,
    )
    res = run_univariate(cfg)
    print(f"\nDone. {len(res['written'])} files in {res['outdir']}")


if __name__ == "__main__":
    main()
