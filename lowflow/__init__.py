"""Low-flow (drought) frequency analysis: marginals, copulas, figures.

Univariate and bivariate frequency analysis of annual and seasonal minimum
river discharge. See README.md for the workflow and docs/CONVENTIONS.md for
the conventions and the failure modes they guard against.

Typical use::

    from lowflow import read_flow, block_minima, best_season, fit_all, compare

    rec = read_flow("my_gauge.xlsx", name="My gauge", D=7)
    bm  = block_minima(rec, season=best_season(rec), min_coverage=0.9)
    fits = fit_all(bm.values)
    print(compare(fits))
"""

from .blocks import (
    BlockMinima,
    best_season,
    block_minima,
    scan_seasons,
    season_diagnostic,
)
from .dataio import FlowRecord, ZeroPolicy, read_flow
from .marginals import (
    GammaDist,
    GEVMinima,
    GEVMinimaNS,
    PearsonIII,
    WeibullDist,
    ad_test,
    compare,
    fit_all,
    ks_test_naive,
    lmom_ratios,
    lr_test,
    trend_tests,
)

__all__ = [
    "read_flow",
    "FlowRecord",
    "ZeroPolicy",
    "block_minima",
    "BlockMinima",
    "season_diagnostic",
    "scan_seasons",
    "best_season",
    "GEVMinima",
    "GEVMinimaNS",
    "PearsonIII",
    "GammaDist",
    "WeibullDist",
    "fit_all",
    "compare",
    "ad_test",
    "ks_test_naive",
    "lr_test",
    "lmom_ratios",
    "trend_tests",
]

__version__ = "1.1.0"


def check_install(search_root=None, verbose=True):
    """Report WHICH copy of lowflow is imported and whether it carries the fixes.

    Several copies of this package commonly end up on one machine (a zip
    extracted twice, a notebooks/ bundle, an older download). ``import lowflow``
    silently takes whichever one is found first on ``sys.path``, and a stale copy
    reproduces bugs that were fixed in the current one -- with no error to say so.

    Parameters
    ----------
    search_root
        Directory to scan for OTHER copies. Defaults to the user's home.

    Returns a dict; prints a verdict when ``verbose``.
    """
    import inspect
    import pathlib
    import sys

    from . import dataio, pipeline

    here = pathlib.Path(__file__).resolve().parent
    checks = {
        "timezone normalised on read": "timezone" in inspect.getsource(dataio.read_flow),
        "empty-block guard raises": "No block is usable for BOTH rivers"
                                    in inspect.getsource(pipeline.run_bivariate),
        "joint records masked not reindexed": "set_axis"
                                              in inspect.getsource(pipeline._joint_records),
    }
    others = []
    root = pathlib.Path(search_root) if search_root else pathlib.Path.home()
    try:
        for cand in root.rglob("lowflow/__init__.py"):
            d = cand.resolve().parent
            if d != here:
                others.append(str(d))
    except (PermissionError, OSError):
        pass

    info = dict(version=__version__, imported_from=str(here),
                checks=checks, other_copies=sorted(others)[:20])
    if verbose:
        print(f"lowflow {__version__}")
        print(f"imported from : {here}")
        for k, v in checks.items():
            print(f"  [{'ok ' if v else 'NO '}] {k}")
        if not all(checks.values()):
            print("\n  >> This copy is STALE. Replace this directory with the current")
            print("     lowflow/ and restart the kernel (Kernel > Restart).")
        if others:
            print(f"\n  {len(others)} other copy/copies of lowflow on this machine:")
            for o in sorted(others)[:10]:
                print("    ", o)
            print("  `import lowflow` takes whichever is first on sys.path.")
        elif all(checks.values()):
            print("\n  Current, and the only copy found.")
    return info
