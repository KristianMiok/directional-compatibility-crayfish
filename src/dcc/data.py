"""Loaders for the World of Crayfish (WoC) export and the GeoFRESH river network.

Schema assumptions are intentionally minimal here. Once we see the real
exports, tighten these by adding column dtypes and a pandera/pydantic schema.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd


# --- WoC occurrences --------------------------------------------------------

# Expected (minimum) columns in the WoC export. Update once we have the real file.
WOC_REQUIRED_COLUMNS = {
    "Crayfish_scientific_name",  # species name
    "long_or",                   # original longitude
    "lat_or",                    # original latitude
    "Year_of_record",
    "Accuracy",                  # "High" / "Low"
    "distance_m",                # snap distance to GeoFRESH
    "subc_id",                   # GeoFRESH sub-catchment / segment ID
    "basin_id",
    "Status",                    # Native / Type locality / Alien / Introduced
}


def load_woc(path: str | Path) -> pd.DataFrame:
    """Load the WoC occurrence table and validate required columns."""
    df = pd.read_csv(path, low_memory=False)
    missing = WOC_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"WoC export at {path} is missing required columns: {sorted(missing)}. "
            "Update WOC_REQUIRED_COLUMNS in src/dcc/data.py if the schema has changed."
        )
    return df


def filter_records(
    df: pd.DataFrame,
    *,
    high_accuracy_only: bool = True,
    max_snap_distance_m: float | None = 200,
) -> pd.DataFrame:
    """Apply the standard quality filters used throughout Stage 1.

    Real-data column names: `Accuracy` ("High"/"Low") and `distance_m`.
    """
    out = df.copy()
    if high_accuracy_only:
        out = out[out["Accuracy"].astype(str).str.lower().eq("high")]
    if max_snap_distance_m is not None:
        out = out[out["distance_m"].le(max_snap_distance_m)]
    return out.reset_index(drop=True)


def deduplicate_by_segment(df: pd.DataFrame) -> pd.DataFrame:
    """One record per (species, subc_id). Override of the earlier signature
    so the rest of the codebase keeps working with real column names."""
    return df.drop_duplicates(subset=["Crayfish_scientific_name", "subc_id"]).reset_index(drop=True)





# --- GeoFRESH network -------------------------------------------------------

def load_geofresh_segments(path: str | Path) -> gpd.GeoDataFrame:
    """Load the GeoFRESH segment layer (gpkg/geojson/shp)."""
    return gpd.read_file(path)
