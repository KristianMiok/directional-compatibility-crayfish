#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
PROJDIR = STAGE3 / "projections"
OUTDIR = ROOT / "data" / "processed" / "stage5" / "supplementary_analyses"
OUTDIR.mkdir(parents=True, exist_ok=True)

MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 2:
        return np.nan

    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()

    sx = rx.std(ddof=1)
    sy = ry.std(ddof=1)
    if sx == 0 or sy == 0:
        return np.nan

    return float(np.corrcoef(rx, ry)[0, 1])


def top20_overlap(v1: pd.Series, v2: pd.Series, top_n: int = 20) -> int:
    s1 = set(v1.sort_values(ascending=False).head(top_n).index)
    s2 = set(v2.sort_values(ascending=False).head(top_n).index)
    return int(len(s1 & s2))


def available_datasets(h5: h5py.File) -> list[str]:
    keys = []

    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset):
            keys.append(name)

    h5.visititems(visitor)
    return sorted(keys)


def resolve_metric_dataset(h5: h5py.File, candidates: list[str]) -> str:
    all_keys = available_datasets(h5)
    lower_map = {k.lower(): k for k in all_keys}

    for cand in candidates:
        if cand in h5:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]

    raise KeyError(
        f"Could not find any of these dataset names: {candidates}\n"
        f"Available datasets in H5:\n" + "\n".join(all_keys)
    )


def load_matrix_metrics(h5_fp: Path) -> tuple[list[str], dict[str, np.ndarray], dict[str, str]]:
    metric_aliases = {
        "mean_suitability": [
            "mean_suitability",
        ],
        "fraction_above": [
            "fraction_above",
            "fraction_above_threshold",
            "fraction_above_thr",
            "fraction_above_0.5",
        ],
        "directional_schoeners_d": [
            "directional_schoeners_d",
            "directional_schoener_d",
            "schoeners_d_directional",
            "schoener_d_directional",
            "schoeners_d",
            "schoener_d",
            "schoener_D",
        ],
        "directional_warrens_i": [
            "directional_warrens_i",
            "directional_warren_i",
            "warrens_i_directional",
            "warren_i_directional",
            "warrens_i",
            "warren_i",
            "warren_I",
        ],
    }

    with h5py.File(h5_fp, "r") as h5:
        species = [
            s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)
            for s in h5["species"][:]
        ]

        metrics = {}
        resolved = {}
        for metric_name, candidates in metric_aliases.items():
            ds_name = resolve_metric_dataset(h5, candidates)
            arr = np.array(h5[ds_name][:], dtype=float)
            if arr.shape != (len(species), len(species)):
                raise ValueError(
                    f"Dataset '{ds_name}' has shape {arr.shape}, expected {(len(species), len(species))}"
                )
            metrics[metric_name] = arr
            resolved[metric_name] = ds_name

    return species, metrics, resolved


def load_q90_metric(species: list[str], proj_dir: Path) -> pd.Series:
    q90_vals = {}

    for source_species in species:
        slug = (
            source_species.lower()
            .replace("(", "")
            .replace(")", "")
            .replace(",", "")
            .replace(".", "")
            .replace("-", "_")
            .replace("/", "_")
            .replace(" ", "_")
        )
        fp = proj_dir / f"{slug}__full_geofresh__gbm.parquet"
        if not fp.exists():
            raise FileNotFoundError(f"Missing projection parquet: {fp}")

        df = pd.read_parquet(fp)
        if "target_species" not in df.columns or "predicted_suitability" not in df.columns:
            raise ValueError(f"Parquet missing required columns: {fp}")

        grp = df.groupby("target_species")["predicted_suitability"].quantile(0.90)
        row = pd.Series([float(grp.get(sp, np.nan)) for sp in species], index=species)
        q90_vals[source_species] = row

    q90_mat = pd.DataFrame(q90_vals).T.loc[species, species]
    q90_mat.index = species
    q90_mat.columns = species

    return q90_mat


def build_species_metric_tables(
    species: list[str],
    matrix_metrics: dict[str, np.ndarray],
    q90_mat: pd.DataFrame,
) -> tuple[dict[str, pd.Series], dict[str, pd.Series]]:
    outgoing = {}
    incoming = {}

    for metric_name, mat in matrix_metrics.items():
        outgoing[metric_name] = pd.Series(np.nanmean(mat, axis=1), index=species, name=metric_name)
        incoming[metric_name] = pd.Series(np.nanmean(mat, axis=0), index=species, name=metric_name)

    outgoing["quantile_90_suitability"] = pd.Series(
        np.nanmean(q90_mat.to_numpy(dtype=float), axis=1),
        index=species,
        name="quantile_90_suitability",
    )
    incoming["quantile_90_suitability"] = pd.Series(
        np.nanmean(q90_mat.to_numpy(dtype=float), axis=0),
        index=species,
        name="quantile_90_suitability",
    )

    return outgoing, incoming


def matrix_from_series_dict(series_dict: dict[str, pd.Series], mode: str, top_n: int = 20) -> pd.DataFrame:
    metric_names = list(series_dict.keys())
    out = pd.DataFrame(index=metric_names, columns=metric_names, dtype=float)

    for m1 in metric_names:
        for m2 in metric_names:
            if mode == "overlap":
                out.loc[m1, m2] = top20_overlap(series_dict[m1], series_dict[m2], top_n=top_n)
            elif mode == "spearman":
                out.loc[m1, m2] = spearman_corr(
                    series_dict[m1].loc[series_dict[m2].index].to_numpy(),
                    series_dict[m2].to_numpy(),
                )
            else:
                raise ValueError(f"Unknown mode: {mode}")

    return out


def offdiag_mean(df: pd.DataFrame) -> float:
    arr = df.to_numpy(dtype=float)
    mask = ~np.eye(arr.shape[0], dtype=bool)
    vals = arr[mask]
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if len(vals) else np.nan


def main() -> None:
    species, matrix_metrics, resolved_names = load_matrix_metrics(MATRIX_FP)

    print("Loading q90 from projection parquets...")
    q90_mat = load_q90_metric(species, PROJDIR)

    outgoing, incoming = build_species_metric_tables(species, matrix_metrics, q90_mat)

    donor_top20_overlap = matrix_from_series_dict(outgoing, mode="overlap", top_n=20).astype(int)
    acceptor_top20_overlap = matrix_from_series_dict(incoming, mode="overlap", top_n=20).astype(int)
    spearman_outgoing = matrix_from_series_dict(outgoing, mode="spearman")
    spearman_incoming = matrix_from_series_dict(incoming, mode="spearman")

    donor_fp = OUTDIR / "metric_sensitivity5_donor_top20_overlap.csv"
    acceptor_fp = OUTDIR / "metric_sensitivity5_acceptor_top20_overlap.csv"
    out_spear_fp = OUTDIR / "metric_sensitivity5_spearman_outgoing.csv"
    in_spear_fp = OUTDIR / "metric_sensitivity5_spearman_incoming.csv"
    summary_fp = OUTDIR / "metric_sensitivity5_summary.json"

    donor_top20_overlap.to_csv(donor_fp)
    acceptor_top20_overlap.to_csv(acceptor_fp)
    spearman_outgoing.to_csv(out_spear_fp, float_format="%.6f")
    spearman_incoming.to_csv(in_spear_fp, float_format="%.6f")

    summary = {
        "input_h5": str(MATRIX_FP),
        "projection_dir": str(PROJDIR),
        "n_species": int(len(species)),
        "resolved_metric_datasets": resolved_names,
        "added_metric": "quantile_90_suitability",
        "top_n": 20,
        "donor_top20_overlap_offdiag_mean": float(offdiag_mean(donor_top20_overlap)),
        "acceptor_top20_overlap_offdiag_mean": float(offdiag_mean(acceptor_top20_overlap)),
        "spearman_outgoing_offdiag_mean": float(offdiag_mean(spearman_outgoing)),
        "spearman_incoming_offdiag_mean": float(offdiag_mean(spearman_incoming)),
    }
    summary_fp.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {donor_fp}")
    print(f"  {acceptor_fp}")
    print(f"  {out_spear_fp}")
    print(f"  {in_spear_fp}")
    print(f"  {summary_fp}")

    print("\nQuick summary:")
    print(f"  donor top20 overlap offdiag mean: {summary['donor_top20_overlap_offdiag_mean']:.3f}")
    print(f"  acceptor top20 overlap offdiag mean: {summary['acceptor_top20_overlap_offdiag_mean']:.3f}")
    print(f"  outgoing Spearman offdiag mean: {summary['spearman_outgoing_offdiag_mean']:.3f}")
    print(f"  incoming Spearman offdiag mean: {summary['spearman_incoming_offdiag_mean']:.3f}")

    print("\nquantile_90_suitability vs other metrics (outgoing Spearman):")
    print(spearman_outgoing.loc["quantile_90_suitability"].to_string())

    print("\nquantile_90_suitability vs other metrics (incoming Spearman):")
    print(spearman_incoming.loc["quantile_90_suitability"].to_string())


if __name__ == "__main__":
    main()
