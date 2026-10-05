# Input data

**No discharge records are committed to this repository.** Gauged series
generally come from a monitoring agency under terms that do not allow
redistribution, so `.gitignore` excludes every workbook in this directory:
your own records cannot be committed by accident.

## Getting a runnable dataset

    python scripts/00_make_synthetic_data.py

writes `station_a.xlsx`, `station_b.xlsx` and `projection.xlsx` here: two
gauges in a shared basin (a large river and a smaller tributary) plus a
transient projection series. These are **synthetic**, generated from a fixed
seed, and built to exercise every path in the pipeline: a seasonal cycle with
late-summer minima, AR(1) persistence, a climate signal shared between the two
gauges, a multi-month outage, a year with only three months of record,
scattered missing days and a handful of zero-flow days. Every number a script
prints from them is a property of the generator, not of any real river.

Replace them with your own workbooks and set the station names in the scripts
to match. Keep the `date` / `flow` column names, or pass `date_col=` and
`flow_col=`.

## Expected format

One row per day, in an `.xlsx` or `.csv` file with two columns:

| column | content |
| --- | --- |
| `date` | a date, any format `pandas.to_datetime` accepts |
| `flow` | mean daily discharge, m3/s; blank or `NaN` for a missing day |

Other column names are fine — pass `date_col=` / `flow_col=` to `read_flow`.
Gaps may be absent rows or present-but-empty rows; `read_flow` reindexes onto a
complete daily calendar either way, so a missing day is never silently treated
as a dry day.

One trap is worth knowing before you export. If your `date` column carries a
UTC offset (`2003-01-01T00:00:00+01:00`) it parses to a timezone-**aware**
index, while a plain Excel date column does not; two records that disagree can
never be aligned, and the joint analysis sees zero common days. `read_flow`
drops the timezone on entry, keeping the local wall-clock date rather than
converting to UTC, and reports what it did in `.summary()`.
