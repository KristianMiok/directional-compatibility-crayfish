"""Stage 3c — Build directional compatibility matrices from projection HDF5 files.

Reads the per-(predictor_set, core) projection tensors from Stage 3b and
collapses them into N x N compatibility matrices. Four metrics per variant:

  mean_suitability[A, B] — mean predicted suitability of species A's SDM
                           over the target cells of species B
  fraction_above[A, B]   — fraction of B's cells where A's prediction is
                           >= A's per-species maxSSS threshold
  schoener_D[A, B]       — directional Schoener's D between A's predictions
                           over B's cells and B's predictions over B's cells
  warren_I[A, B]         — directional Warren's I, same structure

Variants produced:
  1. full_geofresh_gbm   — primary matrix, 155 x 155
  2. climate_local_gbm   — Piece 2 comparison partner, 156 x 156
  3. climate_local_ensemble — mean(GBM, GLM) where both pass QC per species,
                               GBM alone otherwise. 156 x 156.

Plus the common 155-subset of variants (1) and (2) for the piece 2 comparison
(Lucian's constraint: identical rosters when comparing climate_local vs
full_geofresh restructuring).

Per-species maxSSS thresholds come from the Stage 3a metadata JSONs.

Output: one HDF5 per variant at data/processed/stage3/matrices_*.h5
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE3 = PROJECT_ROOT / "data" / "processed" / "stage3"
ARTIFACTS = STAGE3 / "sdm_artifacts"
OUTPUT_DIR = STAGE3


def slugify(name: str) -> str:
    return re.sub(r"^_+|_+$", "", re.sub(r"[^a-z0-9]+", "_", name.lower()))


def load_maxsss_lookup(predictor_set: str, core: str) -> dict[str, float]:
    """Species name -> mean maxSSS threshold for a given (predictor_set, core).

    For Cambarus carolinus on full_geofresh (failed QC), falls back to the
    local_full metadata since that's where its SDM actually lives.
    """
    lookup = {}
    for meta_file in ARTIFACTS.glob(f"*/{predictor_set}__{core.lower()}_metadata.json"):
        meta = json.loads(meta_file.read_text())
        if meta.get("qc_pass") and meta.get("maxsss_mean") is not None:
            lookup[meta["species"]] = float(meta["maxsss_mean"])

    # Special case: Cambarus carolinus on full_geofresh uses local_full fallback
    if predictor_set == "full_geofresh":
        fallback = ARTIFACTS / "cambarus_carolinus" / "local_full__gbm_metadata.json"
        if fallback.exists():
            meta = json.loads(fallback.read_text())
            if meta.get("qc_pass") and meta.get("maxsss_mean") is not None:
                lookup["Cambarus carolinus"] = float(meta["maxsss_mean"])

    return lookup


def build_target_cell_index(h5):
    """Return dict: target_species -> np.ndarray of column indices belonging to it."""
    target_species_arr = h5["target_species"][:].astype(str)
    idx = {}
    for i, sp in enumerate(target_species_arr):
        idx.setdefault(sp, []).append(i)
    return {sp: np.array(cols, dtype=np.int64) for sp, cols in idx.items()}


def directional_schoener_warren(p_source: np.ndarray,
                                 p_target_self: np.ndarray) -> tuple[float, float]:
    """Directional niche overlap between source A and target B on B's cells.

    p_source        = A's predicted suitability over B's cells, shape (K_B,)
    p_target_self   = B's own predicted suitability over B's cells, shape (K_B,)

    Schoener's D = 1 - 0.5 * sum(|P_A - P_B|)  after normalizing each to sum to 1
    Warren's I   = 1 - 0.5 * sum((sqrt(P_A) - sqrt(P_B))^2) / 2  (Hellinger-based)

    We apply the normalization so the values live in [0, 1] and are directly
    comparable to the literature. NaN inputs are skipped; if fewer than 2 valid
    cells remain, returns (NaN, NaN).
    """
    mask = np.isfinite(p_source) & np.isfinite(p_target_self)
    if mask.sum() < 2:
        return float("nan"), float("nan")

    a = p_source[mask].astype(np.float64)
    b = p_target_self[mask].astype(np.float64)

    sum_a = a.sum()
    sum_b = b.sum()
    if sum_a <= 0 or sum_b <= 0:
        return float("nan"), float("nan")

    pa = a / sum_a
    pb = b / sum_b

    D = 1.0 - 0.5 * np.abs(pa - pb).sum()
    H_sq = 0.5 * ((np.sqrt(pa) - np.sqrt(pb)) ** 2).sum()
    I = 1.0 - H_sq

    return float(D), float(I)


def build_matrices(h5_path: Path, variant_name: str,
                   predictor_set: str, core: str,
                   maxsss_lookup: dict[str, float]) -> Path:
    """Build the 4 matrices for one (predictor_set, core) variant."""
    print(f"\n=== Building matrices for {variant_name} ===")
    print(f"  source: {h5_path.name}")

    with h5py.File(h5_path, "r") as h5:
        source_species = [s.decode() for s in h5["source_species"][:]]
        all_target_species = [s.decode() for s in h5["target_species"][:]]
        target_cell_idx_map = build_target_cell_index(h5)
        predictions = h5["predictions"][:]  # (N_source, K_total_cells)

    N = len(source_species)
    # Target species roster = source roster (matrix is square; source and target
    # are the same set of species, each being one row and one column)
    target_list = source_species[:]
    assert set(target_list) <= set(target_cell_idx_map.keys()), \
        "source species missing from target cells"

    # For directional D/I we need p_target_self — each species' self-projection
    # (row A_i's prediction on cells belonging to target i). Cache once:
    self_predictions = {}
    for i, sp in enumerate(target_list):
        cols = target_cell_idx_map[sp]
        self_predictions[sp] = predictions[i, cols]

    mean_suit = np.full((N, N), np.nan, dtype=np.float32)
    frac_above = np.full((N, N), np.nan, dtype=np.float32)
    schoener = np.full((N, N), np.nan, dtype=np.float32)
    warren = np.full((N, N), np.nan, dtype=np.float32)

    thresholds = np.array([maxsss_lookup.get(sp, np.nan) for sp in source_species],
                          dtype=np.float32)
    missing_thresh = [sp for sp in source_species if sp not in maxsss_lookup]
    if missing_thresh:
        print(f"  WARN: {len(missing_thresh)} species have no maxSSS threshold "
              f"(fraction_above will be NaN for those rows):")
        for sp in missing_thresh[:5]:
            print(f"    {sp}")

    for i, source_sp in enumerate(source_species):
        row = predictions[i]  # (K_total_cells,)
        thresh = thresholds[i]
        for j, target_sp in enumerate(target_list):
            cols = target_cell_idx_map[target_sp]
            p_A_on_B = row[cols]  # A's predictions on B's cells
            p_B_on_B = self_predictions[target_sp]

            valid = np.isfinite(p_A_on_B)
            n_valid = int(valid.sum())
            if n_valid == 0:
                continue

            mean_suit[i, j] = float(np.nanmean(p_A_on_B))
            if np.isfinite(thresh):
                frac_above[i, j] = float((p_A_on_B[valid] >= thresh).mean())
            D, I = directional_schoener_warren(p_A_on_B, p_B_on_B)
            schoener[i, j] = D
            warren[i, j] = I

        if (i + 1) % 25 == 0 or i == N - 1:
            print(f"  [{i+1}/{N}] {source_sp}")

    out_path = OUTPUT_DIR / f"matrices_{variant_name}.h5"
    with h5py.File(out_path, "w") as h5:
        h5.attrs["variant"] = variant_name
        h5.attrs["predictor_set"] = predictor_set
        h5.attrs["core"] = core
        h5.attrs["n_species"] = N
        h5.attrs["description"] = (
            "Directional compatibility matrices. "
            "Row = source species A (whose SDM was trained). "
            "Col = target species B (whose native cells are the projection domain). "
            "M[A, B] quantifies how compatible A's niche is with B's realized environment."
        )
        h5.create_dataset("species", data=np.array(source_species, dtype="S80"))
        h5.create_dataset("mean_suitability", data=mean_suit,
                          compression="gzip", compression_opts=4)
        h5.create_dataset("fraction_above", data=frac_above,
                          compression="gzip", compression_opts=4)
        h5.create_dataset("schoener_D", data=schoener,
                          compression="gzip", compression_opts=4)
        h5.create_dataset("warren_I", data=warren,
                          compression="gzip", compression_opts=4)
        h5.create_dataset("maxsss_threshold", data=thresholds)

    size_mb = out_path.stat().st_size / 1e6
    print(f"  Wrote {out_path.name} ({size_mb:.2f} MB)")
    return out_path


def build_ensemble_matrices(gbm_h5_path: Path, glm_h5_path: Path) -> Path:
    """climate_local ensemble = per-cell mean(GBM, GLM) where both species pass
    QC for both cores, GBM alone for species where GLM failed QC.

    Lucian's rule (P2+checkpoint): AUC/TSS drives ensemble membership. GLM is
    dropped for species that failed AUC or TSS on GLM/climate_local; Boyce
    reliability does NOT drive ensemble membership.
    """
    print("\n=== Building climate_local ensemble matrices ===")

    # Which species have qc_pass on climate_local/GLM?
    glm_qc_pass = set()
    for meta_file in ARTIFACTS.glob("*/climate_local__glm_metadata.json"):
        meta = json.loads(meta_file.read_text())
        if meta.get("qc_pass"):
            glm_qc_pass.add(meta["species"])

    # Load GBM projections + GLM projections, average per-cell where both exist
    with h5py.File(gbm_h5_path, "r") as g:
        source_gbm = [s.decode() for s in g["source_species"][:]]
        target_species_gbm = [s.decode() for s in g["target_species"][:]]
        preds_gbm = g["predictions"][:]

    with h5py.File(glm_h5_path, "r") as g:
        source_glm = [s.decode() for s in g["source_species"][:]]
        preds_glm = g["predictions"][:]

    # Align source species order (GLM may have different roster)
    assert source_gbm == source_glm, \
        "GBM and GLM have different source species rosters — need reindexing"

    use_glm = np.array([sp in glm_qc_pass for sp in source_gbm], dtype=bool)
    print(f"  Species using ensemble: {use_glm.sum()}/{len(source_gbm)}")
    print(f"  Species using GBM alone: {(~use_glm).sum()}/{len(source_gbm)}")

    ensemble = preds_gbm.copy()
    # Where GLM also valid, take mean
    ensemble[use_glm] = (preds_gbm[use_glm] + preds_glm[use_glm]) / 2.0

    # Write as temp projections file so build_matrices can consume it
    temp_path = OUTPUT_DIR / "projections_climate_local_ensemble.h5"
    with h5py.File(gbm_h5_path, "r") as src:
        with h5py.File(temp_path, "w") as dst:
            for k, v in src.attrs.items():
                dst.attrs[k] = v
            dst.attrs["core"] = "ENSEMBLE"
            for key in ["source_species", "target_species", "target_cell_idx",
                        "target_cell_subc_id"]:
                dst.create_dataset(key, data=src[key][:])
            dst.create_dataset("predictions", data=ensemble.astype(np.float32),
                               compression="gzip", compression_opts=4)
    print(f"  Wrote intermediate ensemble projections: {temp_path.name}")

    # Use GBM's maxSSS as the ensemble threshold (Lucian approved equal-weight
    # ensemble; using GBM threshold for fraction_above is a defensible choice
    # since GBM drives for 1 of 156 species anyway)
    maxsss_gbm = load_maxsss_lookup("climate_local", "GBM")
    return build_matrices(temp_path, "climate_local_ensemble",
                          "climate_local", "ENSEMBLE", maxsss_gbm)


def build_common_subset(full_gf_h5: Path, climate_h5: Path) -> Path:
    """Piece 2 comparison partner: climate_local_gbm restricted to the same
    155-species roster as full_geofresh_gbm. Stage 4 Mantel test and
    restructuring analyses consume this + the full_geofresh matrix.
    """
    print("\n=== Building climate_local_gbm_common155 (piece 2 comparison) ===")

    with h5py.File(full_gf_h5, "r") as g:
        species_155 = [s.decode() for s in g["species"][:]]

    with h5py.File(climate_h5, "r") as g:
        species_156 = [s.decode() for s in g["species"][:]]
        index_map = np.array([species_156.index(sp) for sp in species_155],
                             dtype=np.int64)

        out_path = OUTPUT_DIR / "matrices_climate_local_gbm_common155.h5"
        with h5py.File(out_path, "w") as dst:
            for k, v in g.attrs.items():
                dst.attrs[k] = v
            dst.attrs["variant"] = "climate_local_gbm_common155"
            dst.attrs["n_species"] = 155
            dst.attrs["description"] = (
                "Piece 2 comparison partner: climate_local_gbm matrices reduced to "
                "the common 155-species subset shared with full_geofresh_gbm. Use "
                "these matrices (not the 156-species version) for Mantel tests and "
                "climate-vs-full restructuring analyses."
            )
            dst.create_dataset("species", data=np.array(species_155, dtype="S80"))
            for metric in ["mean_suitability", "fraction_above", "schoener_D", "warren_I"]:
                M = g[metric][:]
                M_sub = M[np.ix_(index_map, index_map)]
                dst.create_dataset(metric, data=M_sub,
                                   compression="gzip", compression_opts=4)
            thresh_sub = g["maxsss_threshold"][:][index_map]
            dst.create_dataset("maxsss_threshold", data=thresh_sub)

    size_mb = out_path.stat().st_size / 1e6
    print(f"  Wrote {out_path.name} ({size_mb:.2f} MB) — 155 species")
    return out_path


def main() -> int:
    maxsss_fg = load_maxsss_lookup("full_geofresh", "GBM")
    maxsss_cl_gbm = load_maxsss_lookup("climate_local", "GBM")
    print(f"maxSSS thresholds loaded: full_geofresh/GBM={len(maxsss_fg)}, "
          f"climate_local/GBM={len(maxsss_cl_gbm)}")

    full_gf_matrix = build_matrices(
        STAGE3 / "projections_full_geofresh_gbm.h5",
        "full_geofresh_gbm", "full_geofresh", "GBM", maxsss_fg)

    climate_matrix = build_matrices(
        STAGE3 / "projections_climate_local_gbm.h5",
        "climate_local_gbm", "climate_local", "GBM", maxsss_cl_gbm)

    ensemble_matrix = build_ensemble_matrices(
        STAGE3 / "projections_climate_local_gbm.h5",
        STAGE3 / "projections_climate_local_glm.h5")

    common155_matrix = build_common_subset(full_gf_matrix, climate_matrix)

    print("\n" + "=" * 70)
    print("Stage 3c outputs:")
    for p in [full_gf_matrix, climate_matrix, ensemble_matrix, common155_matrix]:
        print(f"  {p.name}  ({p.stat().st_size / 1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())