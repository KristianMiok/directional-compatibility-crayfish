#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
OUTDIR = STAGE3 / "tables"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_matrix(path: Path):
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mean_suitability = h5["mean_suitability"][:].astype(np.float64)
    return species, mean_suitability


def offdiag_row_col_means(mat: np.ndarray):
    n = mat.shape[0]
    row_means = np.empty(n, dtype=np.float64)
    col_means = np.empty(n, dtype=np.float64)

    for i in range(n):
        row_vals = np.delete(mat[i, :], i)
        col_vals = np.delete(mat[:, i], i)
        row_means[i] = np.nanmean(row_vals)
        col_means[i] = np.nanmean(col_vals)

    return row_means, col_means


def top_donors_acceptors(species, mat, top_n=20):
    row_means, col_means = offdiag_row_col_means(mat)

    donor_order = np.argsort(-row_means)
    acceptor_order = np.argsort(-col_means)

    donors = pd.DataFrame({
        "rank": np.arange(1, top_n + 1),
        "species": np.array(species)[donor_order][:top_n],
        "score": row_means[donor_order][:top_n],
        "role": "donor",
    })

    acceptors = pd.DataFrame({
        "rank": np.arange(1, top_n + 1),
        "species": np.array(species)[acceptor_order][:top_n],
        "score": col_means[acceptor_order][:top_n],
        "role": "acceptor",
    })

    return donors, acceptors, row_means, col_means


def top_asymmetric_pairs(species, mat, top_n=20):
    n = len(species)
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            a_to_b = float(mat[i, j])
            b_to_a = float(mat[j, i])
            rows.append({
                "species_a": species[i],
                "species_b": species[j],
                "a_to_b": a_to_b,
                "b_to_a": b_to_a,
                "abs_diff": abs(a_to_b - b_to_a),
                "signed_diff_a_minus_b": a_to_b - b_to_a,
            })

    df = (
        pd.DataFrame(rows)
        .sort_values("abs_diff", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def top_repositioning(species_fg, mat_fg, species_cl, mat_cl, top_n=20):
    if species_fg != species_cl:
        raise ValueError("Species rosters do not match between full_geofresh and climate_local_common155.")

    fg_out, _ = offdiag_row_col_means(mat_fg)
    cl_out, _ = offdiag_row_col_means(mat_cl)

    df = pd.DataFrame({
        "species": species_fg,
        "full_geofresh_outgoing_mean": fg_out,
        "climate_local_outgoing_mean": cl_out,
        "shift_full_minus_climate": fg_out - cl_out,
        "abs_shift": np.abs(fg_out - cl_out),
    }).sort_values("abs_shift", ascending=False).head(top_n).reset_index(drop=True)

    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def main():
    fg_species, fg_mat = load_matrix(STAGE3 / "matrices_full_geofresh_gbm.h5")
    cl155_species, cl155_mat = load_matrix(STAGE3 / "matrices_climate_local_gbm_common155.h5")

    donors, acceptors, _, _ = top_donors_acceptors(fg_species, fg_mat, top_n=20)
    asym_pairs = top_asymmetric_pairs(fg_species, fg_mat, top_n=20)
    reposition = top_repositioning(fg_species, fg_mat, cl155_species, cl155_mat, top_n=20)

    donors_acceptors = pd.concat([donors, acceptors], ignore_index=True)

    donors_acceptors.to_csv(OUTDIR / "table_stage3_top20_donors_acceptors.csv", index=False)
    asym_pairs.to_csv(OUTDIR / "table_stage3_top20_asymmetric_pairs.csv", index=False)
    reposition.to_csv(OUTDIR / "table_stage3_top20_repositioning.csv", index=False)

    print("Wrote:")
    print(f"  {OUTDIR / 'table_stage3_top20_donors_acceptors.csv'}")
    print(f"  {OUTDIR / 'table_stage3_top20_asymmetric_pairs.csv'}")
    print(f"  {OUTDIR / 'table_stage3_top20_repositioning.csv'}")

    print("\nTop 10 donors:")
    print(donors.head(10).to_string(index=False))
    print("\nTop 10 acceptors:")
    print(acceptors.head(10).to_string(index=False))
    print("\nTop 10 asymmetric pairs:")
    print(asym_pairs.head(10).to_string(index=False))
    print("\nTop 10 repositioning species:")
    print(reposition.head(10).to_string(index=False))


if __name__ == "__main__":
    main()