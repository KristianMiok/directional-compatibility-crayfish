"""Task 1.1 — species inventory.

Produces one row per species with counts at each filter stage plus
geographic and temporal summaries.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint

from .data import deduplicate_by_segment, filter_records


INVENTORY_COLUMNS = [
    "species_name",
    "family",
    "genus",
    "total_records",
    "records_high_accuracy",
    "records_snapped_le_200m",
    "records_deduplicated_segment",
    "n_basins",
    "n_countries",
    "n_continents",
    "range_bbox_area_km2",
    "native_only_flag",
    "temporal_span_years",
    "temporal_density",
    "ecological_strategy",  # filled later if GeoTraits is linkable
]


def _convex_hull_area_km2(lon: pd.Series, lat: pd.Series) -> float:
    """Convex hull area of a point cloud, returned in km²."""
    pts = list(zip(lon.dropna(), lat.dropna(), strict=True))
    if len(pts) < 3:
        return 0.0
    hull = MultiPoint(pts).convex_hull
    # Reproject to an equal-area CRS (World Cylindrical Equal Area, EPSG:6933)
    gdf = gpd.GeoSeries([hull], crs="EPSG:4326").to_crs("EPSG:6933")
    return float(gdf.area.iloc[0]) / 1e6  # m² -> km²


def _native_only_flag(group: pd.DataFrame) -> str:
    """Return 'True' / 'False' / 'uncertain' based on the native_status column."""
    statuses = set(group["native_status"].astype(str).str.lower().unique())
    if statuses <= {"native"}:
        return "True"
    if "invasive" in statuses or "non-native" in statuses:
        return "False"
    return "uncertain"


def build_inventory(woc: pd.DataFrame, *, max_snap_m: float = 200) -> pd.DataFrame:
    """Build the species inventory table from a raw WoC dataframe."""
    rows: list[dict] = []
    for species, group in woc.groupby("species_name", sort=True):
        total = len(group)

        high_acc = filter_records(group, high_accuracy_only=True, max_snap_distance_m=None)
        snapped = filter_records(group, high_accuracy_only=True, max_snap_distance_m=max_snap_m)
        dedup = deduplicate_by_segment(snapped)

        years = pd.to_numeric(group["year"], errors="coerce").dropna()
        span = int(years.max() - years.min()) if not years.empty else 0
        density = float(len(group) / span) if span > 0 else float("nan")

        rows.append({
            "species_name": species,
            "family": group["family"].mode().iat[0] if not group["family"].mode().empty else None,
            "genus": group["genus"].mode().iat[0] if not group["genus"].mode().empty else None,
            "total_records": total,
            "records_high_accuracy": len(high_acc),
            "records_snapped_le_200m": len(snapped),
            "records_deduplicated_segment": len(dedup),
            "n_basins": int(group["basin_id"].nunique()),
            "n_countries": int(group["country"].nunique()),
            "n_continents": int(group["continent"].nunique()),
            "range_bbox_area_km2": _convex_hull_area_km2(group["longitude"], group["latitude"]),
            "native_only_flag": _native_only_flag(group),
            "temporal_span_years": span,
            "temporal_density": density,
            "ecological_strategy": np.nan,  # populated later from GeoTraits if available
        })

    return pd.DataFrame(rows, columns=INVENTORY_COLUMNS)
