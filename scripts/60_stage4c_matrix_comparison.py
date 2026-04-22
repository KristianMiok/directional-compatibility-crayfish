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
OUTDIR = STAGE3 / "stage4c"
OUTDIR.mkdir(parents=True, exist_ok=True)

FULL_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"
CLIM155_FP = STAGE3 / "matrices_climate_local_gbm_common155.h5"


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def summarize_perm_test(observed: float, null_vals: np.ndarray, n_perm: int) -> dict:
    p_ge = (np.sum(null_vals >= observed) + 1.0) / (n_perm + 1.0)
    p_le = (np.sum(null_vals <= observed) + 1.0) / (n_perm + 1.0)
    p_two = min(1.0, 2.0 * min(p_ge, p_le))
    return {
        "observed_spearman": float(observed),
        "null_mean": float(np.mean(null_vals)),
        "null_sd": float(np.std(null_vals, ddof=1)),
        "null_q025": float(np.quantile(null_vals, 0.025)),
        "null_q50": float(np.quantile(null_vals, 0.5)),
        "null_q975": float(np.quantile(null_vals, 0.975)),
        "p_ge": float(p_ge),
        "p_le": float(p_le),
        "p_two_sided": float(p_two),
        "n_perm": int(n_perm),
    }


def directed_offdiag_correlation_test(
    a: np.ndarray,
    b: np.ndarray,
    n_perm: int = 999,
    seed: int = 42,
) -> dict:
    if a.shape != b.shape:
        raise ValueError("Matrices must have same shape.")

    mask = ~np.eye(a.shape[0], dtype=bool)
    av = a[mask]
    bv = b[mask]
    obs = spearmanr(av, bv).statistic

    rng = np.random.default_rng(seed)
    null_vals = np.empty(n_perm, dtype=np.float64)

    for i in range(n_perm):
        perm = rng.permutation(a.shape[0])
        bp = b[perm][:, perm]
        null_vals[i] = spearmanr(av, bp[mask]).statistic
        if (i + 1) % 100 == 0 or (i + 1) == n_perm:
            print(f"  Directed permutations: {i+1}/{n_perm}")

    out = summarize_perm_test(obs, null_vals, n_perm)
    out["comparison_type"] = "directed_offdiag"
    return out


def symmetric_uppertriangle_correlation_test(
    a: np.ndarray,
    b: np.ndarray,
    n_perm: int = 999,
    seed: int = 42,
) -> dict:
    if a.shape != b.shape:
        raise ValueError("Matrices must have same shape.")

    a_sym = 0.5 * (a + a.T)
    b_sym = 0.5 * (b + b.T)

    tri = np.triu_indices_from(a_sym, k=1)
    av = a_sym[tri]
    bv = b_sym[tri]
    obs = spearmanr(av, bv).statistic

    rng = np.random.default_rng(seed)
    null_vals = np.empty(n_perm, dtype=np.float64)

    for i in range(n_perm):
        perm = rng.permutation(a.shape[0])
        bp = b_sym[perm][:, perm]
        null_vals[i] = spearmanr(av, bp[tri]).statistic
        if (i + 1) % 100 == 0 or (i + 1) == n_perm:
            print(f"  Symmetric permutations: {i+1}/{n_perm}")

    out = summarize_perm_test(obs, null_vals, n_perm)
    out["comparison_type"] = "symmetric_uppertriangle"
    return out


def pairwise_residuals(
    species: list[str],
    full_mat: np.ndarray,
    clim_mat: np.ndarray,
    top_n: int = 50,
) -> pd.DataFrame:
    rows = []
    n = len(species)

    for i in range(n):
        for j in range(i + 1, n):
            full_ij = float(full_mat[i, j])
            clim_ij = float(clim_mat[i, j])
            diff_ij = full_ij - clim_ij

            full_ji = float(full_mat[j, i])
            clim_ji = float(clim_mat[j, i])
            diff_ji = full_ji - clim_ji

            rows.append({
                "species_a": species[i],
                "species_b": species[j],
                "direction": f"{species[i]} -> {species[j]}",
                "full_geofresh": full_ij,
                "climate_local": clim_ij,
                "diff_full_minus_climate": diff_ij,
                "abs_diff": abs(diff_ij),
            })
            rows.append({
                "species_a": species[j],
                "species_b": species[i],
                "direction": f"{species[j]} -> {species[i]}",
                "full_geofresh": full_ji,
                "climate_local": clim_ji,
                "diff_full_minus_climate": diff_ji,
                "abs_diff": abs(diff_ji),
            })

    df = pd.DataFrame(rows).sort_values("abs_diff", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df.head(top_n).copy()


def undirected_pair_residuals(
    species: list[str],
    full_mat: np.ndarray,
    clim_mat: np.ndarray,
    top_n: int = 50,
) -> pd.DataFrame:
    rows = []
    n = len(species)

    for i in range(n):
        for j in range(i + 1, n):
            full_mean = float(np.mean([full_mat[i, j], full_mat[j, i]]))
            clim_mean = float(np.mean([clim_mat[i, j], clim_mat[j, i]]))
            diff = full_mean - clim_mean
            rows.append({
                "species_a": species[i],
                "species_b": species[j],
                "full_geofresh_pair_mean": full_mean,
                "climate_local_pair_mean": clim_mean,
                "diff_full_minus_climate": diff,
                "abs_diff": abs(diff),
            })

    df = pd.DataFrame(rows).sort_values("abs_diff", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df.head(top_n).copy()


def offdiag_row_means(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        out[i] = np.nanmean(np.delete(mat[i, :], i))
    return out


def repositioners(
    species: list[str],
    full_mat: np.ndarray,
    clim_mat: np.ndarray,
) -> pd.DataFrame:
    full_out = offdiag_row_means(full_mat)
    clim_out = offdiag_row_means(clim_mat)

    df = pd.DataFrame({
        "species": species,
        "full_geofresh_outgoing_mean": full_out,
        "climate_local_outgoing_mean": clim_out,
        "shift_full_minus_climate": full_out - clim_out,
        "abs_shift": np.abs(full_out - clim_out),
    }).sort_values("abs_shift", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def species_residual_involvement(
    species: list[str],
    full_mat: np.ndarray,
    clim_mat: np.ndarray,
) -> pd.DataFrame:
    abs_resid = np.abs(full_mat - clim_mat)
    n = len(species)

    scores = []
    for i in range(n):
        vals = np.delete(abs_resid[i, :], i)
        scores.append(np.nanmean(vals))

    df = pd.DataFrame({
        "species": species,
        "mean_abs_pairwise_residual": scores,
    }).sort_values("mean_abs_pairwise_residual", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def compare_reposition_vs_residual(
    reposition_df: pd.DataFrame,
    residual_df: pd.DataFrame,
    top_n: int = 10,
) -> dict:
    merged = reposition_df.merge(
        residual_df,
        on="species",
        how="inner",
        suffixes=("_reposition", "_residual"),
    )

    rho = spearmanr(
        merged["abs_shift"].to_numpy(),
        merged["mean_abs_pairwise_residual"].to_numpy(),
    ).statistic

    top_rep = set(reposition_df.head(top_n)["species"])
    top_res = set(residual_df.head(top_n)["species"])
    overlap = sorted(top_rep & top_res)

    return {
        "spearman_abs_shift_vs_mean_abs_pairwise_residual": float(rho),
        "top_n": int(top_n),
        "overlap_count": int(len(overlap)),
        "overlap_species": overlap,
    }


def main() -> None:
    full_species, full_mat = load_matrix(FULL_FP)
    clim_species, clim_mat = load_matrix(CLIM155_FP)

    if full_species != clim_species:
        raise RuntimeError("Species rosters do not match between full_geofresh and climate_local_common155.")

    print("Running directed off-diagonal permutation correlation test...")
    directed_test = directed_offdiag_correlation_test(full_mat, clim_mat, n_perm=999, seed=42)

    print("\nRunning symmetric upper-triangle permutation correlation test...")
    symmetric_test = symmetric_uppertriangle_correlation_test(full_mat, clim_mat, n_perm=999, seed=4242)

    print("\nBuilding residual tables...")
    residual_pairs = pairwise_residuals(full_species, full_mat, clim_mat, top_n=50)
    undirected_pairs = undirected_pair_residuals(full_species, full_mat, clim_mat, top_n=50)
    reposition_df = repositioners(full_species, full_mat, clim_mat)
    residual_species_df = species_residual_involvement(full_species, full_mat, clim_mat)
    relation = compare_reposition_vs_residual(reposition_df, residual_species_df, top_n=10)

    residual_pairs.to_csv(OUTDIR / "table_stage4c_top50_shifted_pairs_directed.csv", index=False)
    undirected_pairs.to_csv(OUTDIR / "table_stage4c_top50_shifted_pairs_undirected.csv", index=False)
    reposition_df.head(50).to_csv(OUTDIR / "table_stage4c_top50_repositioners.csv", index=False)
    residual_species_df.head(50).to_csv(OUTDIR / "table_stage4c_top50_species_pairwise_residuals.csv", index=False)

    summary = {
        "directed_offdiag_test": directed_test,
        "symmetric_uppertriangle_test": symmetric_test,
        "reposition_vs_pairwise_residual_relation": relation,
    }
    (OUTDIR / "stage4c_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {OUTDIR / 'table_stage4c_top50_shifted_pairs_directed.csv'}")
    print(f"  {OUTDIR / 'table_stage4c_top50_shifted_pairs_undirected.csv'}")
    print(f"  {OUTDIR / 'table_stage4c_top50_repositioners.csv'}")
    print(f"  {OUTDIR / 'table_stage4c_top50_species_pairwise_residuals.csv'}")
    print(f"  {OUTDIR / 'stage4c_summary.json'}")

    print("\nDirected off-diagonal summary:")
    for k, v in directed_test.items():
        print(f"  {k}: {v}")

    print("\nSymmetric upper-triangle summary:")
    for k, v in symmetric_test.items():
        print(f"  {k}: {v}")

    print("\nTop 10 shifted directed pairs:")
    print(residual_pairs.head(10).to_string(index=False))

    print("\nTop 10 shifted undirected pairs:")
    print(undirected_pairs.head(10).to_string(index=False))

    print("\nTop 10 species by pairwise residual involvement:")
    print(residual_species_df.head(10).to_string(index=False))

    print("\nReposition vs residual relation:")
    for k, v in relation.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()