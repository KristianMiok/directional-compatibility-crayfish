#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
BOOTROOT = STAGE3 / "stage4d"
PRIMARY_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"
OUTDIR = BOOTROOT / "summary"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def offdiag_row_means(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    vals = np.empty(n, dtype=np.float64)
    for i in range(n):
        vals[i] = np.nanmean(np.delete(mat[i, :], i))
    return vals


def offdiag_col_means(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    vals = np.empty(n, dtype=np.float64)
    for i in range(n):
        vals[i] = np.nanmean(np.delete(mat[:, i], i))
    return vals


def top_species(species: list[str], scores: np.ndarray, n: int = 20) -> list[str]:
    order = np.argsort(-scores)
    return [species[i] for i in order[:n]]


def matrix_spearman(a: np.ndarray, b: np.ndarray) -> float:
    mask = ~np.eye(a.shape[0], dtype=bool)
    return float(spearmanr(a[mask], b[mask]).statistic)


def main() -> None:
    primary_species, primary_mat = load_matrix(PRIMARY_FP)
    primary_donors = top_species(primary_species, offdiag_row_means(primary_mat), n=20)
    primary_acceptors = top_species(primary_species, offdiag_col_means(primary_mat), n=20)

    iter_dirs = sorted([p for p in BOOTROOT.glob("iter_*") if p.is_dir()])
    if not iter_dirs:
        raise RuntimeError("No bootstrap iteration directories found.")

    rows = []
    donor_presence = {sp: 0 for sp in primary_species}
    acceptor_presence = {sp: 0 for sp in primary_species}

    n_ok = 0

    for it_dir in iter_dirs:
        mat_fp = it_dir / "matrices_full_geofresh_gbm.h5"
        if not mat_fp.exists():
            print(f"Skipping {it_dir.name}: matrix file missing")
            continue

        species, mat = load_matrix(mat_fp)
        if species != primary_species:
            print(f"Skipping {it_dir.name}: species roster mismatch")
            continue

        donors = top_species(species, offdiag_row_means(mat), n=20)
        acceptors = top_species(species, offdiag_col_means(mat), n=20)

        for sp in donors:
            donor_presence[sp] += 1
        for sp in acceptors:
            acceptor_presence[sp] += 1

        donor_overlap = len(set(primary_donors) & set(donors))
        acceptor_overlap = len(set(primary_acceptors) & set(acceptors))
        rho = matrix_spearman(primary_mat, mat)

        rows.append({
            "iteration": it_dir.name,
            "matrix_spearman_vs_primary": rho,
            "top20_donor_overlap_count": donor_overlap,
            "top20_acceptor_overlap_count": acceptor_overlap,
        })
        n_ok += 1

    if n_ok == 0:
        raise RuntimeError("No bootstrap matrices were found/usable for summary.")

    iter_df = pd.DataFrame(rows).sort_values("iteration").reset_index(drop=True)
    donor_stability = pd.DataFrame({
        "species": primary_species,
        "donor_top20_frequency": [donor_presence[s] / n_ok for s in primary_species],
        "acceptor_top20_frequency": [acceptor_presence[s] / n_ok for s in primary_species],
    })

    unstable = donor_stability[
        (donor_stability["donor_top20_frequency"] > 0) & (donor_stability["donor_top20_frequency"] < 0.5)
        | (donor_stability["acceptor_top20_frequency"] > 0) & (donor_stability["acceptor_top20_frequency"] < 0.5)
    ].copy()

    iter_df.to_csv(OUTDIR / "table_stage4d_iteration_stability.csv", index=False)
    donor_stability.to_csv(OUTDIR / "table_stage4d_species_top20_frequencies.csv", index=False)
    unstable.to_csv(OUTDIR / "table_stage4d_unstable_species_flags.csv", index=False)

    summary = {
        "n_iterations_used": int(n_ok),
        "matrix_spearman_mean": float(iter_df["matrix_spearman_vs_primary"].mean()),
        "matrix_spearman_sd": float(iter_df["matrix_spearman_vs_primary"].std(ddof=1) if len(iter_df) > 1 else 0.0),
        "donor_overlap_mean": float(iter_df["top20_donor_overlap_count"].mean()),
        "acceptor_overlap_mean": float(iter_df["top20_acceptor_overlap_count"].mean()),
    }
    (OUTDIR / "stage4d_summary.json").write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUTDIR / 'table_stage4d_iteration_stability.csv'}")
    print(f"  {OUTDIR / 'table_stage4d_species_top20_frequencies.csv'}")
    print(f"  {OUTDIR / 'table_stage4d_unstable_species_flags.csv'}")
    print(f"  {OUTDIR / 'stage4d_summary.json'}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nIteration stability:")
    print(iter_df.to_string(index=False))


if __name__ == "__main__":
    main()