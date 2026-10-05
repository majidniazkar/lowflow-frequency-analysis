"""Generate three synthetic daily-discharge workbooks in ``data/``.

These stand in for your own records so the pipeline runs end to end out of the
box: two gauges in a shared basin (a large river and a smaller tributary) plus
a transient projection series for script 03. The series carry the features the pipeline is built to survive: a seasonal
cycle with late-summer minima, AR(1) persistence, a climate signal shared
between the two rivers, a multi-month gauge outage, a year with only three
months of record, scattered missing days, and a handful of zero-flow days.

Replace ``data/*.xlsx`` with your own records -- same two columns, ``date`` and
``flow`` -- and point the scripts at them. Nothing downstream is specific to
these series.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SEED = 20240617


def seasonal_log(doy, lo_doy, hi, lo):
    ph = 2 * np.pi * (doy - lo_doy) / 365.25
    return np.log(lo) + (np.log(hi) - np.log(lo)) * 0.5 * (1 - np.cos(ph))


def ar1(n, rho, sd, rng):
    e = rng.normal(0, sd * np.sqrt(1 - rho**2), n)
    x = np.empty(n)
    x[0] = rng.normal(0, sd)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + e[i]
    return x


def make_river(dates, hi, lo, lo_doy, sd, rho, shared, share_w, rng):
    doy = dates.dayofyear.to_numpy()
    base = seasonal_log(doy, lo_doy, hi, lo)
    idio = ar1(len(dates), rho, sd, rng)
    yr = dates.year.to_numpy()
    ann = np.array([shared[y] for y in yr])
    return np.exp(base + share_w * ann + np.sqrt(max(0.0, 1 - share_w**2)) * idio)


def main(outdir="data"):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    dates = pd.date_range("1970-01-01", "2024-12-31", freq="D")
    years = sorted(set(dates.year))
    shared = pd.Series(ar1(len(years), 0.25, 0.38, rng), index=years)
    big = make_river(dates, 2600, 480, 213, 0.34, 0.93, shared, 0.75, rng)
    small = make_river(dates, 95, 4.2, 220, 0.52, 0.90, shared, 0.62, rng)

    a_df = pd.DataFrame({"date": dates, "flow": np.round(big, 1)})
    b_df = pd.DataFrame({"date": dates, "flow": np.round(small, 2)})
    b_df.loc[(dates >= "1983-04-10") & (dates <= "1983-11-02"), "flow"] = np.nan
    a_df.loc[(dates.year == 1991) & (dates.month > 3), "flow"] = np.nan
    b_df.loc[rng.choice(np.flatnonzero(b_df.flow.notna()), 140, replace=False), "flow"] = np.nan
    dry = b_df.flow.notna() & (b_df.flow < 1.0)
    b_df.loc[b_df.index[dry][:6], "flow"] = 0.0

    pdates = pd.date_range("2015-01-01", "2100-12-31", freq="D")
    pyears = sorted(set(pdates.year))
    trend = pd.Series(np.linspace(0, -0.55, len(pyears)) + ar1(len(pyears), 0.2, 0.33, rng),
                      index=pyears)
    proj = make_river(pdates, 90, 3.9, 224, 0.55, 0.90, trend, 0.70, rng)
    proj_df = pd.DataFrame({"date": pdates, "flow": np.round(proj, 2)})

    for df, name in [(a_df, "station_a"), (b_df, "station_b"),
                     (proj_df, "projection")]:
        df.to_excel(f"{outdir}/{name}.xlsx", index=False)
        print(f"wrote {outdir}/{name}.xlsx  ({len(df)} rows, "
              f"{int(df.flow.isna().sum())} missing)")


if __name__ == "__main__":
    main()
