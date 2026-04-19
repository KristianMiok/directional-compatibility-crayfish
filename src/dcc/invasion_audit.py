"""Task 1.5 — invasion contamination audit.

Flags records that may come from invaded ranges (rather than native ones)
so Lucian / Mihaela / Dave can manually validate them. The output is the
input to a human review loop, not an automated decision.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_native_ranges(path: str | Path) -> dict[str, dict[str, set[str]]]:
    """Load a native-ranges CSV.

    Expected columns: species_name, native_continents (semicolon-separated),
    native_basins (semicolon-separated, optional).
    """
    df = pd.read_csv(path)
    out: dict[str, dict[str, set[str]]] = {}
    for _, r in df.iterrows():
        out[r["species_name"]] = {
            "continents": set(str(r.get("native_continents", "")).split(";")) - {""},
            "basins": set(str(r.get("native_basins", "")).split(";")) - {""},
        }
    return out


def flag_records(
    woc: pd.DataFrame,
    native_ranges: dict[str, dict[str, set[str]]],
    *,
    flag_outside_continent: bool = True,
    post_year: int | None = 1900,
) -> pd.DataFrame:
    """Return rows flagged for manual review, with a `reason` column."""
    flagged = []
    for _, r in woc.iterrows():
        sp = r["species_name"]
        nat = native_ranges.get(sp)
        if nat is None:
            # No native-range info for this species — flag everything as 'unknown_native_range'
            flagged.append({**r.to_dict(), "reason": "unknown_native_range"})
            continue

        reasons = []
        if flag_outside_continent and nat["continents"] and r["continent"] not in nat["continents"]:
            reasons.append("outside_native_continent")
        if (
            nat["basins"]
            and r["basin_id"] not in nat["basins"]
            and post_year is not None
            and pd.notna(r.get("year"))
            and float(r["year"]) >= post_year
        ):
            reasons.append("outside_native_basin_post_year")
        if r.get("native_status", "").lower() in {"invasive", "non-native"}:
            reasons.append("source_marked_invasive")

        if reasons:
            flagged.append({**r.to_dict(), "reason": "|".join(reasons)})

    return pd.DataFrame(flagged)
