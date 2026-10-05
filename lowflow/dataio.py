"""Reading, regularising and screening daily discharge records.

This module exists because three of the findings in the R review were data-layer
problems, not statistics problems:

* the D-day moving average was computed *after* invalid rows had been deleted, so
  the window silently averaged across calendar gaps (R finding C3);
* no per-year coverage requirement, so a year with 30 valid days contributed a
  "minimum" on equal footing with a year with 365 (C2);
* zero flows were dropped by ``flow > 0``, which removes the most informative
  observations in a drought study and biases the low tail upward (C1).

The order of operations here is deliberate and is the whole point of the module:
**regularise to a gap-free daily index first, then roll, then screen.**
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = ["FlowRecord", "read_flow", "ZeroPolicy"]


@dataclass(frozen=True)
class ZeroPolicy:
    """What to do with zero / below-detection discharge.

    Parameters
    ----------
    mode
        ``"keep"``   -- keep zeros as genuine observations (the default, and the
                        right choice for low-flow work: a zero is the most
                        extreme low flow the gauge can record).
        ``"censor"`` -- replace zeros by ``floor`` and record them as
                        left-censored. Use when the fitted distribution has
                        support on (0, inf) and a hard zero breaks the
                        likelihood, e.g. a two-parameter Weibull.
        ``"drop"``   -- delete them. This is what the R scripts did. Kept only
                        so the old behaviour can be reproduced on demand; it
                        biases every low-flow statistic upward.
    floor
        Replacement value used by ``"censor"``.
    """

    mode: str = "keep"
    floor: float = 0.001

    def __post_init__(self) -> None:
        if self.mode not in {"keep", "censor", "drop"}:
            raise ValueError(f"unknown zero policy {self.mode!r}")


@dataclass
class FlowRecord:
    """A regularised daily discharge record plus its screening diagnostics.

    Attributes
    ----------
    daily
        Gap-free daily frame indexed by date, columns ``flow`` and ``flow_D``
        (the D-day right-aligned mean; identical to ``flow`` when ``D == 1``).
    name, D
        Station label and averaging window actually used.
    diagnostics
        Everything the screening steps discovered, so it can be printed in a
        methods section rather than left implicit.
    """

    daily: pd.DataFrame
    name: str
    D: int
    diagnostics: dict = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover - convenience only
        d = self.diagnostics
        return (
            f"<FlowRecord {self.name}: {d.get('n_days_span')} days "
            f"{d.get('first_date')}..{d.get('last_date')}, "
            f"{d.get('n_valid')} valid, D={self.D}>"
        )

    def summary(self) -> str:
        d = self.diagnostics
        lines = [
            f"Record            : {self.name}",
            f"Span              : {d['first_date']} to {d['last_date']} "
            f"({d['n_days_span']} calendar days)",
            f"Valid observations: {d['n_valid']} "
            f"({100 * d['n_valid'] / d['n_days_span']:.1f}% of span)",
            f"Missing in span   : {d['n_missing_in_span']}",
            f"Duplicate dates   : {d['n_duplicate_dates']}",
            f"Negative flows    : {d['n_negative']} (set to NaN)",
            f"Zero flows        : {d['n_zero']} (policy: {d['zero_policy']})",
            f"Averaging window  : D = {self.D} day(s)",
        ]
        if d.get("timezone"):
            lines.append(f"Timezone          : {d['timezone']}")
        if d.get("n_unparseable_dates"):
            lines.append(f"Unparseable dates : {d['n_unparseable_dates']} (dropped)")
        if self.D > 1:
            lines.append(
                f"D-day means lost to gaps: {d['n_rolling_lost_to_gaps']} "
                "(windows spanning a gap are NaN, not averaged across it)"
            )
        return "\n".join(lines)


def read_flow(
    path,
    *,
    name: str | None = None,
    date_col: str = "date",
    flow_col: str = "flow",
    D: int = 1,
    zero_policy: ZeroPolicy | str = "keep",
    sheet_name: int | str = 0,
) -> FlowRecord:
    """Load a daily discharge series and prepare it for block-minimum extraction.

    Accepts ``.xlsx``/``.xls`` (via openpyxl) or ``.csv``. Unlike the R scripts
    this does **not** silently drop rows: invalid values become NaN so that the
    rolling mean can refuse to average across them, and every removal is counted
    in :attr:`FlowRecord.diagnostics`.

    Parameters
    ----------
    D
        Averaging window in days, right-aligned. ``D = 7`` gives the 7-day mean
        underlying the conventional MAM7 / 7Q10 low-flow index; ``D = 1`` is the
        single-day minimum. Windows that would span a missing day yield NaN
        rather than a mean over fewer days.
    zero_policy
        See :class:`ZeroPolicy`. Either an instance or one of the mode strings.
    """
    if isinstance(zero_policy, str):
        zero_policy = ZeroPolicy(zero_policy)
    path = str(path)
    if path.lower().endswith((".xlsx", ".xls")):
        raw = pd.read_excel(path, sheet_name=sheet_name)
    else:
        raw = pd.read_csv(path)

    missing = {date_col, flow_col} - set(raw.columns)
    if missing:
        raise KeyError(
            f"{path}: missing column(s) {sorted(missing)}; found {list(raw.columns)}. "
            "Pass date_col=/flow_col= if your workbook uses different headers."
        )
    df = raw[[date_col, flow_col]].rename(columns={date_col: "date", flow_col: "flow"})
    diag: dict = {}

    # --- timezone normalisation -------------------------------------------
    # A date column exported as ISO strings with a UTC offset
    # ("2003-01-01T00:00:00+01:00") parses to a tz-AWARE index, while a plain
    # Excel date column parses to a tz-naive one. Two records that disagree
    # cannot be aligned: their DatetimeIndexes never compare equal, an index
    # union silently degrades to dtype=object, and every join between them
    # returns zero rows. Normalising here, at the single point of entry, is
    # what keeps that from reaching the rest of the package.
    #
    # Timezones are DROPPED, not converted: a gauge record is timestamped in
    # the station's own local time and the calendar date as recorded is the
    # block label we want. Converting to UTC would shift dates across midnight.
    s = pd.to_datetime(df["date"], errors="coerce")
    tz_note = None
    if isinstance(s.dtype, pd.DatetimeTZDtype):
        tz_note = f"dropped fixed offset {s.dt.tz} (kept local wall-clock date)"
        s = s.dt.tz_localize(None)
    elif s.dtype == object:
        # mixed or inconsistent offsets: the only consistent reading is UTC
        s = pd.to_datetime(df["date"], errors="coerce", utc=True).dt.tz_localize(None)
        tz_note = "mixed UTC offsets in the date column; converted to UTC then dropped"
    diag["timezone"] = tz_note
    n_unparsed = int(s.isna().sum())
    diag["n_unparseable_dates"] = n_unparsed
    df["date"] = s.dt.normalize()
    df = df.loc[df["date"].notna()]
    if df.empty:
        raise ValueError(f"{path}: no parseable dates in column {date_col!r}")
    df["flow"] = pd.to_numeric(df["flow"], errors="coerce")
    n_dup = int(df.duplicated("date").sum())
    if n_dup:
        df = df.groupby("date", as_index=False)["flow"].mean()
    diag["n_duplicate_dates"] = n_dup

    neg = df["flow"] < 0
    diag["n_negative"] = int(neg.sum())
    df.loc[neg, "flow"] = np.nan

    zeros = df["flow"] == 0
    diag["n_zero"] = int(zeros.sum())
    diag["zero_policy"] = zero_policy.mode
    if zero_policy.mode == "drop":
        df.loc[zeros, "flow"] = np.nan
    elif zero_policy.mode == "censor":
        df.loc[zeros, "flow"] = zero_policy.floor

    # --- the fix for R finding C3: regularise BEFORE rolling -----------------
    full = pd.date_range(df["date"].min(), df["date"].max(), freq="D")
    daily = df.set_index("date").reindex(full).rename_axis("date")
    diag["first_date"] = str(full[0].date())
    diag["last_date"] = str(full[-1].date())
    diag["n_days_span"] = len(full)
    diag["n_valid"] = int(daily["flow"].notna().sum())
    diag["n_missing_in_span"] = int(daily["flow"].isna().sum())

    if D < 1 or int(D) != D:
        raise ValueError(f"D must be a positive integer, got {D!r}")
    D = int(D)
    if D == 1:
        daily["flow_D"] = daily["flow"]
        diag["n_rolling_lost_to_gaps"] = 0
    else:
        # min_periods=D is what makes a window containing a gap return NaN
        daily["flow_D"] = daily["flow"].rolling(window=D, min_periods=D).mean()
        enough_room = np.arange(len(daily)) >= D - 1
        lost = daily["flow_D"].isna() & daily["flow"].notna() & enough_room
        diag["n_rolling_lost_to_gaps"] = int(lost.sum())

    return FlowRecord(daily=daily, name=name or path, D=D, diagnostics=diag)
