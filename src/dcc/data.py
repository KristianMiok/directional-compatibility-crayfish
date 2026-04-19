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
    "species_name",
    "family",
    "genus",
    "longitude",
    "latitude",
    "year",
    "accuracy_flag",      # high / medium / low (or boolean)
    "snap_distance_m",    # distance to nearest GeoFRESH segment after snapping
    "segment_id",         # GeoFRESH segment the record was snapped to
    "basin_id",
    "country",
    "continent",
    "source",             # citation / dataset of origin
    "native_status",      # native / invasive / unknown
}


def load_woc(path: str | Path) -> pd.DataFrame:
    """Load the WoC occurrence table and validate required columns."""
    df = pd.read_csv(path)
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
    """Apply the standard quality filters used throughout Stage 1."""
    out = df.copy()
    if high_accuracy_only:
        # Accept either boolean or string accuracy flags.
        if out["accuracy_flag"].dtype == bool:
            out = out[out["accuracy_flag"]]
        else:
            out = out[out["accuracy_flag"].astype(str).str.lower().eq("high")]
    if max_snap_distance_m is not None:
        out = out[out["snap_distance_m"].le(max_snap_distance_m)]
    return out.reset_index(drop=True)


def deduplicate_by_segment(df: pd.DataFrame) -> pd.DataFrame:
    """One record per (species, segment_id). Benchmark count for thresholds."""
    return df.drop_duplicates(subset=["species_name", "segment_id"]).reset_index(drop=True)


# --- GeoFRESH network -------------------------------------------------------

def load_geofresh_segments(path: str | Path) -> gpd.GeoDataFrame:
    """Load the GeoFRESH segment layer (gpkg/geojson/shp)."""
    return gpd.read_file(path)
