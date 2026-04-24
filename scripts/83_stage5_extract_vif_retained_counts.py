#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "data" / "processed" / "stage3" / "sdm_artifacts"
PRIMARY_MATRIX_FP = ROOT / "data" / "processed" / "stage3" / "matrices_full_geofresh_gbm.h5"
OUTDIR = ROOT / "data" / "processed" / "stage5" / "supplementary_analyses"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = OUTDIR / "vif_retained_per_species.csv"
OUT_JSON = OUTDIR / "vif_retained_per_species_summary.json"


def load_metadata(fp: Path) -> dict:
    return json.loads(fp.read_text())


def load_primary_species(h5_fp: Path) -> list[str]:
    with h5py.File(h5_fp, "r") as h5:
        return [
            s.decode("utf-8") if isinstance(s, (bytes, bytes)) else str(s)
            for s in h5["species"][:]
        ]


def main() -> None:
    primary_species = load_primary_species(PRIMARY_MATRIX_FP)
    primary_species_set = set(primary_species)

    rows = []
    species_dirs = sorted([p for p in ARTIFACTS.iterdir() if p.is_dir()])

    for sp_dir in species_dirs:
        climate_fp = sp_dir / "climate_local__gbm_metadata.json"
        full_fp = sp_dir / "full_geofresh__gbm_metadata.json"

        if not climate_fp.exists() or not full_fp.exists():
            continue

        climate = load_metadata(climate_fp)
        full = load_metadata(full_fp)

        species_name = str(climate.get("species") or full.get("species") or sp_dir.name)

        if species_name not in primary_species_set:
            continue

        rows.append(
            {
                "species": species_name,
                "n_vars_climate_local": int(climate["n_predictors"]),
                "n_vars_full_geofresh": int(full["n_predictors"]),
            }
        )

    if not rows:
        raise RuntimeError("No matching species metadata found for the primary 155-species cohort.")

    df = pd.DataFrame(rows)

    # preserve the exact primary matrix species order
    df["species"] = pd.Categorical(df["species"], categories=primary_species, ordered=True)
    df = df.sort_values("species").reset_index(drop=True)
    df["species"] = df["species"].astype(str)

    df.to_csv(OUT_CSV, index=False)

    summary = {
        "n_species_primary_matrix": int(len(primary_species)),
        "n_species_with_both_metadata": int(len(df)),
        "n_vars_climate_local_unique": sorted([int(x) for x in df["n_vars_climate_local"].dropna().unique()]),
        "n_vars_full_geofresh_unique": sorted([int(x) for x in df["n_vars_full_geofresh"].dropna().unique()]),
        "n_vars_climate_local_mean": float(df["n_vars_climate_local"].mean()),
        "n_vars_full_geofresh_mean": float(df["n_vars_full_geofresh"].mean()),
        "n_vars_climate_local_sd": float(df["n_vars_climate_local"].std(ddof=1)),
        "n_vars_full_geofresh_sd": float(df["n_vars_full_geofresh"].std(ddof=1)),
        "all_primary_species_matched": bool(len(df) == len(primary_species)),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_CSV}")
    print(f"  {OUT_JSON}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nFirst 10 rows:")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
