"""Generate a small synthetic World of Crayfish (WoC) export for smoke-testing Stage 1.

This produces a CSV that matches `WOC_REQUIRED_COLUMNS` from src/dcc/data.py.
It is NOT real ecological data — coordinates are scattered uniformly inside
crude continental boxes and "basins" are fictitious. The only purpose is to
exercise every step of the pipeline end-to-end so we can find wiring bugs
before the real WoC export arrives.

Run from the project root:

    uv run python fixtures/generate_synthetic_woc.py

Writes to:
    data/raw/woc_occurrences.csv
    data/external/native_ranges.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

RNG = np.random.default_rng(20260419)  # deterministic — change seed if you want variation

# 8 fictitious species, designed to exercise different code paths:
#   - 3 species with ample data (>500 records) — survive the strict threshold
#   - 3 species in the 200–500 band — survive the default threshold
#   - 1 species with ~100 records — only survives the relaxed threshold
#   - 1 species with ~30 records — never survives any threshold (tests "drop" path)
# A few species are cross-continental (one native, one invaded) to exercise
# the invasion audit.

SPECIES_SPEC = [
    # (name, family, genus, native_continent, n_records, also_in_continent_or_None)
    ("Astacus astacus",        "Astacidae",   "Astacus",      "Europe",        650, None),
    ("Pacifastacus leniusculus","Astacidae",  "Pacifastacus", "North America", 800, "Europe"),   # invader in EU
    ("Procambarus clarkii",    "Cambaridae",  "Procambarus",  "North America", 720, "Europe"),   # invader in EU
    ("Cherax quadricarinatus", "Parastacidae","Cherax",       "Oceania",       300, None),
    ("Astacopsis gouldi",      "Parastacidae","Astacopsis",   "Oceania",       250, None),
    ("Cambarus bartonii",      "Cambaridae",  "Cambarus",     "North America", 220, None),
    ("Austropotamobius pallipes","Astacidae", "Austropotamobius","Europe",     110, None),  # only at threshold 80
    ("Engaeus orientalis",     "Parastacidae","Engaeus",      "Oceania",        30, None),  # never retained
]

# Crude bounding boxes (lon_min, lon_max, lat_min, lat_max) for "native" continent samples.
CONTINENT_BBOX = {
    "Europe":        (-10.0,  30.0, 36.0, 60.0),
    "North America": (-125.0,-70.0, 25.0, 50.0),
    "Oceania":       (113.0, 154.0,-43.0,-12.0),
}

# Roughly 6 fictitious basins per continent, just IDs — no geometry needed for Stage 1.
BASINS_PER_CONTINENT = {cont: [f"{cont[:2].upper()}_BAS_{i:02d}" for i in range(1, 7)]
                        for cont in CONTINENT_BBOX}

COUNTRIES_PER_CONTINENT = {
    "Europe":        ["DE", "FR", "PL", "RO", "ES", "IT"],
    "North America": ["US", "CA", "MX"],
    "Oceania":       ["AU", "NZ"],
}


def _sample_records(n: int, continent: str, base_segment_id: int) -> pd.DataFrame:
    """Sample n records uniformly inside the continent bbox."""
    lon_min, lon_max, lat_min, lat_max = CONTINENT_BBOX[continent]
    lon = RNG.uniform(lon_min, lon_max, n)
    lat = RNG.uniform(lat_min, lat_max, n)

    # 85% high accuracy, 15% medium
    accuracy = RNG.choice(["high", "medium"], size=n, p=[0.85, 0.15])
    # Snap distance: most are close, a tail extends past 200 m
    snap = np.clip(RNG.exponential(scale=80, size=n), 1, 5000)
    # One segment per (lon, lat) pseudo-snap — share segments occasionally to test dedup
    segment_ids = base_segment_id + RNG.integers(0, max(n // 2, 1), size=n)
    basins = RNG.choice(BASINS_PER_CONTINENT[continent], size=n)
    countries = RNG.choice(COUNTRIES_PER_CONTINENT[continent], size=n)
    years = RNG.integers(1950, 2025, size=n)

    # Synthetic environmental values per record. Not realistic ranges for the
    # species — just spread enough to give Task 1.3 non-trivial completeness.
    mean_ann_temp_c = RNG.normal(loc=12 + (lat_max - lat) / (lat_max - lat_min) * 8, scale=2)
    ann_precip_mm = RNG.uniform(400, 2200, n)
    q_mean_m3s = np.exp(RNG.normal(loc=2.0, scale=1.2, size=n))
    elev_m = np.clip(RNG.normal(loc=400, scale=300, size=n), 0, 3500)
    slope_pct = np.clip(RNG.exponential(scale=4, size=n), 0, 60)
    lc_dominant = RNG.choice(
        ["forest", "grassland", "cropland", "wetland", "urban"],
        size=n, p=[0.35, 0.2, 0.25, 0.15, 0.05],
    )

    return pd.DataFrame({
        "longitude": lon,
        "latitude": lat,
        "year": years,
        "accuracy_flag": accuracy,
        "snap_distance_m": snap,
        "segment_id": segment_ids,
        "basin_id": basins,
        "country": countries,
        "continent": continent,
        "mean_ann_temp_c": mean_ann_temp_c,
        "ann_precip_mm": ann_precip_mm,
        "q_mean_m3s": q_mean_m3s,
        "elev_m": elev_m,
        "slope_pct": slope_pct,
        "lc_dominant": lc_dominant,
    })


def build_dataset() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (occurrences_df, native_ranges_df)."""
    all_records: list[pd.DataFrame] = []
    native_ranges: list[dict] = []
    segment_offset = 1_000_000  # arbitrary

    for name, family, genus, native, n, invaded in SPECIES_SPEC:
        # Native records — these are the bulk
        n_native = int(n * (0.85 if invaded else 1.0))
        df_native = _sample_records(n_native, native, segment_offset)
        df_native["species_name"] = name
        df_native["family"] = family
        df_native["genus"] = genus
        df_native["native_status"] = "native"
        df_native["source"] = f"WoC_synthetic_{name.split()[0].lower()}"
        all_records.append(df_native)
        segment_offset += n_native + 100

        # Invaded-range records — should be flagged in Task 1.5
        if invaded:
            n_inv = n - n_native
            df_inv = _sample_records(n_inv, invaded, segment_offset)
            df_inv["species_name"] = name
            df_inv["family"] = family
            df_inv["genus"] = genus
            df_inv["native_status"] = RNG.choice(["invasive", "unknown"], size=n_inv, p=[0.7, 0.3])
            df_inv["source"] = f"WoC_synthetic_{name.split()[0].lower()}_invaded"
            df_inv["year"] = RNG.integers(1950, 2025, size=n_inv)  # mix pre/post 1900 cutoff
            all_records.append(df_inv)
            segment_offset += n_inv + 100

        native_ranges.append({
            "species_name": name,
            "native_continents": native,  # semicolon-separated if multiple — only one here
            "native_basins": ";".join(BASINS_PER_CONTINENT[native]),
        })

    woc = pd.concat(all_records, ignore_index=True)
    # Re-order to match WOC_REQUIRED_COLUMNS expectations (loose order — pandas is name-based)
    cols = ["species_name", "family", "genus", "longitude", "latitude", "year",
            "accuracy_flag", "snap_distance_m", "segment_id", "basin_id",
            "country", "continent", "source", "native_status",
            "mean_ann_temp_c", "ann_precip_mm", "q_mean_m3s",
            "elev_m", "slope_pct", "lc_dominant"]
    woc = woc[cols]

    nat = pd.DataFrame(native_ranges)
    return woc, nat


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    raw_dir = project_root / "data" / "raw"
    ext_dir = project_root / "data" / "external"
    raw_dir.mkdir(parents=True, exist_ok=True)
    ext_dir.mkdir(parents=True, exist_ok=True)

    woc, nat = build_dataset()
    woc_path = raw_dir / "woc_occurrences.csv"
    nat_path = ext_dir / "native_ranges.csv"
    woc.to_csv(woc_path, index=False)
    nat.to_csv(nat_path, index=False)

    print(f"Wrote {woc_path}")
    print(f"  {len(woc):,} records across {woc['species_name'].nunique()} species")
    print(f"  Records per species:")
    for sp, n in woc["species_name"].value_counts().items():
        print(f"    {sp:35s} {n:5d}")
    print(f"Wrote {nat_path} ({len(nat)} species)")


if __name__ == "__main__":
    main()