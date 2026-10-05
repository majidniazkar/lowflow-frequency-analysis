"""Univariate low-flow frequency analysis of a CLIMATE PROJECTION series.

Port of ``Marginal_Proj.R``, with three changes that alter the conclusions:

1. **Its own output directory.** The R version wrote its figures into the
   observed-record folder under the observational filenames and titles,
   overwriting the observational results with projection results.
2. **Non-stationarity is handled, not merely tested.** A transient projection
   to 2100 is the case where a stationary GEV is least defensible. This run
   fits a trend in the GEV location parameter and compares it against the
   stationary fit by likelihood ratio, and additionally fits fixed epochs.
3. **The D-day switch is live.** In the R version the D-day mean was computed
   and then never used -- the minima came from the raw daily column -- so the
   projection and observed analyses were not comparable at any D > 1.

The season is taken from the **observed** record at the same gauge, not scanned
on the projection, so that the projected and observed analyses describe the same
event. Set ``OBS_PATH`` to that record.

Still missing, and not something code can supply: bias adjustment of the
simulated series against the observed record, targeted at the low tail rather
than the mean, and a statement of whether the file is one realisation or an
ensemble member. Quoting projected low flows as absolute values without the
former is not defensible; report relative change instead. See README.md.

Run:  python scripts/03_marginal_projection.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lowflow import best_season, read_flow
from lowflow.pipeline import UnivariateConfig, run_univariate

# Point these at your own workbooks; the two columns must be date and flow.
PATH = "data/projection.xlsx"        # the simulated series
OBS_PATH = "data/station_b.xlsx"     # the observed record at the same gauge
STATION = "Station B (projection)"

D = 7
MIN_COV = 0.90
SEASON = None          # None = take the season from OBS_PATH; or pin e.g. (7, 8, 9)
EPOCHS = ((2015, 2044), (2045, 2074), (2075, 2100))


def main():
    season = SEASON
    if season is None:
        obs = read_flow(OBS_PATH, name="observed", D=D, zero_policy="keep")
        season = best_season(obs, min_coverage=MIN_COV)
    print(f"Season in use: {season}   (from the observed record, so the two are comparable)")

    cfg = UnivariateConfig(
        path=PATH,
        station=STATION,
        label=f"{STATION}, seasonal minimum {D}-day mean, months {season}",
        outdir="out/03_marginal_projection",   # NOT the observed-record folder
        season=season,
        D=D,
        min_coverage=MIN_COV,
        start_month=1,
        zero_policy="keep",
        T=(2, 5, 10, 20, 50, 100),
        B_gof=499,
        B_ci=999,
        scan_season=True,
        nonstationary=True,                    # fit and test a trend
        epochs=EPOCHS,
    )
    res = run_univariate(cfg)
    print(f"\nDone. {len(res['written'])} files in {res['outdir']}")
    if res["nonstationary"]:
        t = res["nonstationary"]["test"]
        print(f"Trend in GEV location: deviance {t['deviance']:.2f} on {t['df']} df, "
              f"p = {t['p_value']:.4g} -> {t['verdict']}")


if __name__ == "__main__":
    main()
