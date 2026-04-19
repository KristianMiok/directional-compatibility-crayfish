"""Task 1.5 — invasion contamination audit (real-data schema).

The master dataset already carries a Status column — Native / Alien /
Introduced / Type locality. Records flagged here are anything not in the
configured native set. This is intentionally trivial: with high-quality
status labels in the source, we trust them. The output goes to manual
review by Lucian / Mihaela / Dave for any cases that look mislabeled.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def flag_records(
    woc: pd.DataFrame,
    *,
    native_values: list[str] | None = None,
    status_column: str = "Status",
) -> pd.DataFrame:
    """Return rows whose Status is NOT in the native set, with a `reason` column."""
    nv = {v.lower() for v in (native_values or ["Native", "Type locality"])}
    is_native = woc[status_column].astype(str).str.lower().isin(nv)
    flagged = woc[~is_native].copy()
    flagged["reason"] = flagged[status_column].fillna("missing_status")
    return flagged
