"""Task 1.1 — species inventory (real-data schema).

Outputs Lucian's spec columns (species_name, family, genus, etc.). Columns
that cannot be derived from the master dataset (family, genus, country,
continent, ecological_strategy) are written as NA so the structure of the
deliverable is preserved.

Adds `records_after_thinning` per Kristian's request: one record per
(species, subc_id) — i.e., the same as deduplicate_by_segment. With the
real WoC + GeoFRESH dataset, subc_id is the GeoFRESH stream segment and
this matches the common spatial-thinning unit. A more sophisticated grid-
cell thinning (matching the processing_report from prior work) can be
added later if Lucian prefers.
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
    "records_after_thinning",
    "n_basins",
    "n_countries",
    "n_continents",
    "range_bbox_area_km2",
    "native_only_flag",
    "temporal_span_years",
    "temporal_density",
    "ecological_strategy",
]


def _convex_hull_area_km2(lon: pd.Series, lat: pd.Series) -> float:
    """Convex hull area of a point cloud in km², using EPSG:6933 (equal-area)."""
    pts = list(zip(lon.dropna(), lat.dropna(), strict=True))
    if len(pts) < 3:
        return 0.0
    hull = MultiPoint(pts).convex_hull
    gdf = gpd.GeoSeries([hull], crs="EPSG:4326").to_crs("EPSG:6933")
    return float(gdf.area.iloc[0]) / 1e6


def _native_only_flag(group: pd.DataFrame, native_values: set[str]) -> str:
    """Return 'True' / 'False' / 'uncertain' based on Status column."""
    statuses = set(group["Status"].dropna().astype(str).unique())
    if not statuses:
        return "uncertain"
    if statuses <= native_values:
        return "True"
    if statuses & native_values and not statuses <= native_values:
        return "False"  # mix of native and non-native
    return "False"


def build_inventory(
    woc: pd.DataFrame,
    *,
    max_snap_m: float = 200,
    native_values: set[str] | None = None,
) -> pd.DataFrame:
    """Build the species inventory table from the real WoC+GeoFRESH master dataset."""
    if native_values is None:
        native_values = {"Native", "Type locality"}

    rows: list[dict] = []
    sp_col = "Crayfish_scientific_name"
    for species, group in woc.groupby(sp_col, sort=True):
        total = len(group)
        high_acc = filter_records(group, high_accuracy_only=True, max_snap_distance_m=None)
        snapped = filter_records(group, high_accuracy_only=True, max_snap_distance_m=max_snap_m)
        dedup = deduplicate_by_segment(snapped)
        thinned = dedup  # one record per (species, subc_id) = same as dedup with current rule

        years = pd.to_numeric(group["Year_of_record"], errors="coerce").dropna()
        span = int(years.max() - years.min()) if not years.empty else 0
        density = float(len(group) / span) if span > 0 else float("nan")

        rows.append({
            "species_name": species,
            "family": np.nan,
            "genus": np.nan,
            "total_records": total,
            "records_high_accuracy": len(high_acc),
            "records_snapped_le_200m": len(snapped),
            "records_deduplicated_segment": len(dedup),
            "records_after_thinning": len(thinned),
            "n_basins": int(group["basin_id"].nunique()),
            "n_countries": np.nan,
            "n_continents": np.nan,
            "range_bbox_area_km2": _convex_hull_area_km2(group["long_or"], group["lat_or"]),
            "native_only_flag": _native_only_flag(group, native_values),
            "temporal_span_years": span,
            "temporal_density": density,
            "ecological_strategy": np.nan,
        })

    return pd.DataFrame(rows, columns=INVENTORY_COLUMNS)
