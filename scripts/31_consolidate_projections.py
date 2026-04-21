"""Stage 3b consolidation — merge per-source-species parquets into HDF5.

Produces one HDF5 file per (predictor_set, core) containing the full
source-to-target projection tensor. Structure:

  /source_species              (N,) fixed-length strings
  /target_species              (M,) fixed-length strings
  /target_cell_species_idx     (K,) int — which target species each cell belongs to
  /target_cell_subc_id         (K,) int64 — the GeoFRESH subc_id
  /predictions                 (N, K) float32 — predicted suitability
                                A[i, j] = species i's SDM applied to cell j

N = number of source species that have a model for this (predictor_set, core)
M = number of unique target species (same across cores but varies by predictor_set
    because local_full is a superset handling C. carolinus)
K = total target cells (sum of each species' presence-env-vector count)

Each (predictor_set, core) gets one file:
  data/processed/stage3/projections_full_geofresh_gbm.h5
  data/processed/stage3/projections_climate_local_gbm.h5
  data/processed/stage3/projections_climate_local_glm.h5
  (and optionally a local_full_gbm.h5 if we include C. carolinus there too)

Stage 3c reads these and builds N x N compatibility matrices.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE3 = PROJECT_ROOT / "data" / "processed" / "stage3"
PROJECTIONS_DIR = STAGE3 / "projections"
OUTPUT_DIR = STAGE3


def consolidate(predictor_set: str, core: str) -> Path | None:
    core_lower = core.lower()
    pattern = f"*__{predictor_set}__{core_lower}.parquet"
    files = sorted(PROJECTIONS_DIR.glob(pattern))
    if not files:
        print(f"  {predictor_set} / {core}: no parquet files, skipping")
        return None

    print(f"\n{predictor_set} / {core}: {len(files)} source species")

    # First pass: collect unique target species, target cells (in stable order
    # taken from the first file — all projections should share the same target
    # env CSV so ordering is consistent)
    first = pd.read_parquet(files[0])
    target_species = first[["target_species", "target_cell_idx", "subc_id"]].drop_duplicates()
    target_species = target_species.sort_values(["target_species", "target_cell_idx"]).reset_index(drop=True)

    K = len(target_species)
    M = target_species["target_species"].nunique()
    print(f"  target cells: {K} across {M} target species")

    # Build the canonical (target_species, target_cell_idx) -> global_cell_idx map
    key_to_idx = {(r.target_species, r.target_cell_idx): i
                  for i, r in target_species.iterrows()}

    # Build predictions matrix
    source_list = []
    predictions = np.full((len(files), K), np.nan, dtype=np.float32)

    for i, f in enumerate(files):
        df = pd.read_parquet(f)
        source_sp = df["source_species"].iloc[0]
        source_list.append(source_sp)

        # Map each row to its global cell index
        keys = list(zip(df["target_species"], df["target_cell_idx"]))
        indices = np.array([key_to_idx.get(k, -1) for k in keys], dtype=np.int64)
        valid = indices >= 0
        predictions[i, indices[valid]] = df.loc[valid, "predicted_suitability"].astype(np.float32).values

        if i % 25 == 0 or i == len(files) - 1:
            print(f"  [{i+1}/{len(files)}] {source_sp:40s}  "
                  f"non-nan={np.sum(~np.isnan(predictions[i])):6d}/{K}")

    out_path = OUTPUT_DIR / f"projections_{predictor_set}_{core_lower}.h5"
    with h5py.File(out_path, "w") as h5:
        h5.attrs["predictor_set"] = predictor_set
        h5.attrs["core"] = core
        h5.attrs["n_source_species"] = len(source_list)
        h5.attrs["n_target_species"] = M
        h5.attrs["n_target_cells"] = K

        h5.create_dataset("source_species", data=np.array(source_list, dtype="S80"))
        h5.create_dataset("target_species",
                          data=np.array(target_species["target_species"].astype(str).values, dtype="S80"))
        h5.create_dataset("target_cell_idx",
                          data=target_species["target_cell_idx"].astype(np.int32).values)
        h5.create_dataset("target_cell_subc_id",
                          data=target_species["subc_id"].astype(np.int64).values)
        h5.create_dataset("predictions", data=predictions,
                          compression="gzip", compression_opts=4)

    print(f"  Wrote {out_path}  ({out_path.stat().st_size / 1e6:.1f} MB)")
    return out_path


def main() -> int:
    combos = [
        ("full_geofresh", "GBM"),
        ("climate_local", "GBM"),
        ("climate_local", "GLM"),
        ("local_full",    "GBM"),  # C. carolinus and any spillover
        ("local_full",    "GLM"),
    ]
    produced = []
    for ps, core in combos:
        p = consolidate(ps, core)
        if p is not None:
            produced.append(p)

    print("\n" + "=" * 70)
    print(f"Wrote {len(produced)} HDF5 files:")
    for p in produced:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())