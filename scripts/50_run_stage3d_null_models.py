#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
STAGE3_DIR = ROOT / "data" / "processed" / "stage3"
OUT_DIR = STAGE3_DIR / "null_models"
OUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class MatrixBundle:
    path: Path
    species: list[str]
    mean_suitability: np.ndarray
    fraction_above: np.ndarray
    schoener_D: np.ndarray
    warren_I: np.ndarray


def load_matrix_bundle(path: Path) -> MatrixBundle:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        return MatrixBundle(
            path=path,
            species=species,
            mean_suitability=h5["mean_suitability"][:].astype(np.float32),
            fraction_above=h5["fraction_above"][:].astype(np.float32),
            schoener_D=h5["schoener_D"][:].astype(np.float32),
            warren_I=h5["warren_I"][:].astype(np.float32),
        )


def matrix_stats(mat: np.ndarray) -> dict[str, float]:
    n = mat.shape[0]
    diag_mask = np.eye(n, dtype=bool)
    off_mask = ~diag_mask

    diag = np.diag(mat)
    off = mat[off_mask]
    asym = mat - mat.T

    row_gt_col = np.sum(mat > mat.T, axis=1)
    donor_share = float(np.mean(row_gt_col > 0.8 * n))
    acceptor_share = float(np.mean(row_gt_col < 0.2 * n))

    return {
        "diag_mean": float(np.nanmean(diag)),
        "offdiag_mean": float(np.nanmean(off)),
        "diag_minus_offdiag": float(np.nanmean(diag) - np.nanmean(off)),
        "asym_abs_mean": float(np.nanmean(np.abs(asym))),
        "asym_abs_median": float(np.nanmedian(np.abs(asym))),
        "donor_share_gt80": donor_share,
        "acceptor_share_lt20": acceptor_share,
    }


def permute_rows_cols_independently(mat: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    row_perm = rng.permutation(mat.shape[0])
    col_perm = rng.permutation(mat.shape[1])
    return mat[row_perm][:, col_perm]


def summarize_null(observed: float, null_values: np.ndarray) -> dict[str, float]:
    null_values = np.asarray(null_values, dtype=float)
    p_ge = (np.sum(null_values >= observed) + 1.0) / (len(null_values) + 1.0)
    p_le = (np.sum(null_values <= observed) + 1.0) / (len(null_values) + 1.0)
    p_two = min(1.0, 2.0 * min(p_ge, p_le))
    z = (observed - np.mean(null_values)) / np.std(null_values, ddof=1) if np.std(null_values, ddof=1) > 0 else np.nan
    return {
        "observed": float(observed),
        "null_mean": float(np.mean(null_values)),
        "null_sd": float(np.std(null_values, ddof=1)),
        "null_q025": float(np.quantile(null_values, 0.025)),
        "null_q50": float(np.quantile(null_values, 0.5)),
        "null_q975": float(np.quantile(null_values, 0.975)),
        "z_score": float(z),
        "p_ge": float(p_ge),
        "p_le": float(p_le),
        "p_two_sided": float(p_two),
    }


def run_single_matrix_null(
    name: str,
    mat: np.ndarray,
    n_perm: int,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    obs = matrix_stats(mat)

    null_store = {k: np.empty(n_perm, dtype=np.float32) for k in obs.keys()}

    for i in range(n_perm):
        pm = permute_rows_cols_independently(mat, rng)
        st = matrix_stats(pm)
        for k, v in st.items():
            null_store[k][i] = v

        if (i + 1) % 100 == 0 or (i + 1) == n_perm:
            print(f"  [{i+1}/{n_perm}] {name}")

    return {
        "matrix_name": name,
        "n_perm": n_perm,
        "results": {k: summarize_null(obs[k], v) for k, v in null_store.items()},
    }


def aligned_matrix_comparison_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    valid = np.isfinite(a) & np.isfinite(b)
    av = a[valid].ravel()
    bv = b[valid].ravel()

    rho = spearmanr(av, bv).statistic
    return {
        "spearman_rho": float(rho),
        "mean_abs_diff": float(np.mean(np.abs(av - bv))),
        "rmse": float(np.sqrt(np.mean((av - bv) ** 2))),
    }


def run_pairwise_alignment_null(
    name: str,
    a: np.ndarray,
    b: np.ndarray,
    n_perm: int,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    obs = aligned_matrix_comparison_stats(a, b)

    null_store = {k: np.empty(n_perm, dtype=np.float32) for k in obs.keys()}

    for i in range(n_perm):
        perm = rng.permutation(b.shape[0])
        bp = b[perm][:, perm]
        st = aligned_matrix_comparison_stats(a, bp)
        for k, v in st.items():
            null_store[k][i] = v

        if (i + 1) % 100 == 0 or (i + 1) == n_perm:
            print(f"  [{i+1}/{n_perm}] {name}")

    return {
        "comparison_name": name,
        "n_perm": n_perm,
        "results": {k: summarize_null(obs[k], v) for k, v in null_store.items()},
    }


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    fg = load_matrix_bundle(STAGE3_DIR / "matrices_full_geofresh_gbm.h5")
    cl = load_matrix_bundle(STAGE3_DIR / "matrices_climate_local_gbm.h5")
    cl155 = load_matrix_bundle(STAGE3_DIR / "matrices_climate_local_gbm_common155.h5")
    ens = load_matrix_bundle(STAGE3_DIR / "matrices_climate_local_ensemble.h5")

    print("=== Single-matrix null models ===")
    single_jobs = [
        ("full_geofresh_gbm__mean_suitability", fg.mean_suitability),
        ("full_geofresh_gbm__fraction_above", fg.fraction_above),
        ("climate_local_gbm__mean_suitability", cl.mean_suitability),
        ("climate_local_gbm__fraction_above", cl.fraction_above),
        ("climate_local_ensemble__mean_suitability", ens.mean_suitability),
        ("climate_local_ensemble__fraction_above", ens.fraction_above),
    ]

    single_results = []
    for i, (name, mat) in enumerate(single_jobs):
        print(f"Running {name}")
        res = run_single_matrix_null(name, mat, n_perm=args.n_perm, seed=args.seed + i)
        write_json(OUT_DIR / f"{name}__null_summary.json", res)
        single_results.append(res)

    print("=== Pairwise alignment null models ===")
    if fg.species != cl155.species:
        raise RuntimeError("Species rosters do not match between full_geofresh and climate_local common155 matrices.")

    pair_jobs = [
        ("piece2_full_geofresh_vs_climate_local_common155__mean_suitability", fg.mean_suitability, cl155.mean_suitability),
        ("piece2_full_geofresh_vs_climate_local_common155__fraction_above", fg.fraction_above, cl155.fraction_above),
        ("climate_local_gbm_vs_ensemble__mean_suitability", cl.mean_suitability, ens.mean_suitability),
        ("climate_local_gbm_vs_ensemble__fraction_above", cl.fraction_above, ens.fraction_above),
    ]

    pair_results = []
    for i, (name, a, b) in enumerate(pair_jobs, start=100):
        print(f"Running {name}")
        res = run_pairwise_alignment_null(name, a, b, n_perm=args.n_perm, seed=args.seed + i)
        write_json(OUT_DIR / f"{name}__null_summary.json", res)
        pair_results.append(res)

    combined = {
        "n_perm": args.n_perm,
        "seed": args.seed,
        "single_matrix_jobs": [x["matrix_name"] for x in single_results],
        "pairwise_jobs": [x["comparison_name"] for x in pair_results],
    }
    write_json(OUT_DIR / "stage3d_manifest.json", combined)

    print("======================================================================")
    print(f"Wrote null-model summaries to: {OUT_DIR}")


if __name__ == "__main__":
    main()
