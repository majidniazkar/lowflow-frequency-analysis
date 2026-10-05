"""Block-minimum extraction with explicit coverage screening and season checks.

Two findings from the R review are implemented here:

* **C2** -- a block is only used if it has enough valid days. ``min_coverage``
  is a fraction of the days *available in the block*, so the same threshold
  works for a 365-day annual block and a 92-day seasonal one.
* **C4** -- :func:`season_diagnostic` answers the question the R script set up
  but never evaluated: how often does the annual minimum actually fall inside
  the chosen season? :func:`scan_seasons` ranks candidate windows by that
  capture rate, so the season is chosen from the record instead of assumed.
"""

from __future__ import annotations

import calendar
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .dataio import FlowRecord

__all__ = [
    "BlockMinima",
    "block_minima",
    "season_diagnostic",
    "scan_seasons",
    "best_season",
    "MONTH_ABBR",
]

MONTH_ABBR = {m: calendar.month_abbr[m] for m in range(1, 13)}


def season_label(months) -> str:
    """``[6, 7, 8]`` -> ``'JJA'``; non-contiguous sets fall back to a join."""
    months = list(months)
    if not months:
        return "year"
    contiguous = all((months[i] % 12) + 1 == months[i + 1] for i in range(len(months) - 1))
    if contiguous:
        return "".join(MONTH_ABBR[m][0] for m in months)
    return "-".join(MONTH_ABBR[m] for m in months)


@dataclass
class BlockMinima:
    """Block minima plus the screening record that produced them."""

    table: pd.DataFrame  # one row per block, including rejected ones
    name: str
    season: tuple[int, ...] | None
    D: int
    block_kind: str  # "calendar year" or "water year starting <Mon>"

    @property
    def used(self) -> pd.DataFrame:
        """Only the blocks that passed the coverage screen."""
        return self.table.loc[self.table["used"]].reset_index(drop=True)

    @property
    def values(self) -> np.ndarray:
        """The minima actually available for fitting, in record order."""
        return self.used["min_flow"].to_numpy(float)

    @property
    def years(self) -> np.ndarray:
        return self.used["year"].to_numpy(int)

    @property
    def label(self) -> str:
        s = season_label(self.season) if self.season else "annual"
        return f"{self.name} {s} minimum of {self.D}-day mean"

    def summary(self) -> str:
        t = self.table
        rej = t.loc[~t["used"]]
        lines = [
            f"Block definition  : {self.block_kind}, season = "
            f"{season_label(self.season) if self.season else 'full year'}, D = {self.D}",
            f"Blocks in span    : {len(t)}",
            f"Blocks used       : {int(t['used'].sum())}",
            f"Blocks rejected   : {len(rej)}"
            + (f" ({', '.join(str(y) for y in rej['year'])})" if len(rej) else ""),
        ]
        if int(t["used"].sum()):
            u = self.used
            lines += [
                f"Coverage of used  : {u['coverage'].min():.2f} to {u['coverage'].max():.2f}"
                f" (median {u['coverage'].median():.2f})",
                f"Minimum of minima : {u['min_flow'].min():.4g} on {u.loc[u['min_flow'].idxmin(), 'date']}",
                f"Month of minimum  : "
                + ", ".join(
                    f"{MONTH_ABBR[m]}={c}"
                    for m, c in sorted(pd.to_datetime(u["date"]).dt.month.value_counts().items())
                ),
            ]
        return "\n".join(lines)


def _block_year(dates: pd.DatetimeIndex, start_month: int) -> np.ndarray:
    """Water-year label: the year in which the block *starts*."""
    if start_month == 1:
        return dates.year.to_numpy()
    return np.where(dates.month >= start_month, dates.year, dates.year - 1)


def block_minima(
    record: FlowRecord,
    *,
    season: "list[int] | tuple[int, ...] | None" = None,
    min_coverage: float = 0.9,
    start_month: int = 1,
    column: str = "flow_D",
) -> BlockMinima:
    """Extract one minimum per block, screening blocks by data coverage.

    Parameters
    ----------
    season
        Months defining the low-flow window, e.g. ``[7, 8, 9]``. ``None`` uses
        the whole block. Use :func:`scan_seasons` to pick this from the data
        rather than by assumption.
    min_coverage
        Fraction of the days present in the block that must carry a valid value
        for the block's minimum to be used. 0.9 of a 92-day season is ~83 days.
        Set to 0 to reproduce the unscreened R behaviour.
    start_month
        1 for calendar years. For a basin whose low-flow season straddles the
        new year, set this so no drought is split across two blocks. Where the
        minimum falls in late summer, 1 is appropriate and this argument
        exists only to make that choice explicit.
    column
        ``"flow_D"`` (the D-day mean, default) or ``"flow"`` for raw daily.

    Notes
    -----
    Coverage is computed on ``column``, so with ``D > 1`` the windows voided by
    gaps in :func:`~lowflow.dataio.read_flow` correctly count as missing.
    """
    daily = record.daily
    if column not in daily:
        raise KeyError(f"{column!r} not in record; have {list(daily.columns)}")
    d = daily.reset_index()
    idx = pd.DatetimeIndex(pd.to_datetime(d["date"], errors="coerce", utc=False))
    if idx.tz is not None:          # defensive: read_flow normally strips this
        idx = idx.tz_localize(None)
        d["date"] = idx
    d["block_year"] = _block_year(idx, start_month)
    if season is not None:
        season = tuple(int(m) for m in season)
        d = d.loc[pd.DatetimeIndex(d["date"]).month.isin(season)]
        if d.empty:
            raise ValueError(f"season {season} selects no days in {record.name}")

    rows = []
    for year, g in d.groupby("block_year", sort=True):
        vals = g[column]
        n_poss, n_valid = len(g), int(vals.notna().sum())
        cov = n_valid / n_poss if n_poss else 0.0
        if n_valid:
            i = vals.idxmin()
            mn, when = float(vals.loc[i]), g.loc[i, "date"]
        else:
            mn, when = np.nan, pd.NaT
        rows.append(
            dict(
                year=int(year),
                min_flow=mn,
                date=when,
                n_possible=n_poss,
                n_valid=n_valid,
                coverage=cov,
                used=bool(n_valid) and cov >= min_coverage,
            )
        )
    table = pd.DataFrame(rows)
    kind = "calendar year" if start_month == 1 else f"water year starting {MONTH_ABBR[start_month]}"
    return BlockMinima(table=table, name=record.name, season=season, D=record.D, block_kind=kind)


def season_diagnostic(
    record: FlowRecord,
    season,
    *,
    min_coverage: float = 0.9,
    start_month: int = 1,
    column: str = "flow_D",
) -> pd.DataFrame:
    """Compare each block's seasonal minimum against its unrestricted minimum.

    This is the table ``Marginal_Seas.R`` built and then never summarised. The
    returned frame carries ``in_season`` per block; its mean is the **capture
    rate** -- the fraction of annual minima that the season actually contains.
    A low capture rate means the seasonal analysis is not describing the annual
    drought, and no amount of downstream care will fix that.
    """
    ann = block_minima(
        record, season=None, min_coverage=min_coverage, start_month=start_month, column=column
    ).table
    sea = block_minima(
        record, season=season, min_coverage=min_coverage, start_month=start_month, column=column
    ).table
    out = ann.merge(sea, on="year", suffixes=("_annual", "_seasonal"))
    out["month_annual"] = pd.to_datetime(out["date_annual"]).dt.month
    out["in_season"] = out["month_annual"].isin(tuple(int(m) for m in season))
    out["ratio_seasonal_to_annual"] = out["min_flow_seasonal"] / out["min_flow_annual"]
    return out


def scan_seasons(
    record: FlowRecord,
    *,
    widths=(3, 4, 5),
    min_coverage: float = 0.9,
    start_month: int = 1,
    column: str = "flow_D",
) -> pd.DataFrame:
    """Rank every contiguous month window by how often it captures the annual minimum.

    Returns one row per candidate window, sorted by ``capture_rate`` descending.
    Use this before committing to a season. The low-flow window is an
    assumption; this function tests it against the record. See
    :func:`best_season` to take the top-ranked window directly.
    """
    ann = block_minima(
        record, season=None, min_coverage=min_coverage, start_month=start_month, column=column
    ).used
    month_of_min = pd.to_datetime(ann["date"]).dt.month.to_numpy()
    rows = []
    for w in widths:
        for start in range(1, 13):
            months = tuple(((start - 1 + k) % 12) + 1 for k in range(w))
            capture = float(np.isin(month_of_min, months).mean()) if len(month_of_min) else np.nan
            bm = block_minima(
                record,
                season=months,
                min_coverage=min_coverage,
                start_month=start_month,
                column=column,
            )
            rows.append(
                dict(
                    window=season_label(months),
                    months=months,
                    width=w,
                    capture_rate=capture,
                    n_blocks_used=len(bm.used),
                    median_min=float(np.median(bm.values)) if len(bm.values) else np.nan,
                )
            )
    return (
        pd.DataFrame(rows)
        .sort_values(["capture_rate", "width"], ascending=[False, True])
        .reset_index(drop=True)
    )


def best_season(
    record: FlowRecord,
    *,
    widths=(3, 4, 5),
    min_coverage: float = 0.9,
    start_month: int = 1,
    column: str = "flow_D",
    min_capture: float = 0.8,
    tol: float = 0.02,
) -> tuple:
    """Return the month window the record itself supports, as a tuple of months.

    Chooses from :func:`scan_seasons`: among the candidates whose capture rate is
    within ``tol`` of the best one, take the NARROWEST window. The tolerance
    matters because capture rate is monotone in window width -- a wider window
    is a superset and can only capture more minima -- so ranking on capture rate
    alone would almost always return the widest candidate offered, which defeats
    the point of a seasonal analysis. A 4-month window capturing 99% of minima is
    a better description of the drought season than a 5-month one capturing 100%.
    Set ``tol=0`` for the strict highest-capture-rate window.

    Raises ``ValueError`` if no candidate reaches ``min_capture``, because in
    that case no contiguous window describes the annual drought and the
    analysis should be run on annual blocks instead.
    """
    scan = scan_seasons(
        record,
        widths=widths,
        min_coverage=min_coverage,
        start_month=start_month,
        column=column,
    )
    if not len(scan) or not np.isfinite(scan["capture_rate"].to_numpy()).any():
        raise ValueError("season scan produced no usable candidate window")
    best = float(scan["capture_rate"].max())
    near = scan.loc[scan["capture_rate"] >= best - tol]
    top = near.sort_values(["width", "capture_rate"], ascending=[True, False]).iloc[0]
    rate = float(top["capture_rate"])
    if rate < min_capture:
        raise ValueError(
            f"the best window ({top['window']}, capture rate {rate:.2f}) still misses "
            f"the annual minimum in {(1 - rate) * 100:.0f}% of blocks, below the "
            f"min_capture={min_capture:.2f} threshold. No contiguous season describes "
            "the drought on this record -- use annual blocks (season=None), or lower "
            "min_capture deliberately."
        )
    return tuple(top["months"])
