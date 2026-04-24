#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
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
    keys: list[str] = []

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


def load_metrics_from_h5(h5_fp: Path) -> tuple[list[str], dict[str, np.ndarray], dict[str, str]]:
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
        ],
        "directional_warrens_i": [
            "directional_warrens_i",
            "directional_warren_i",
            "warrens_i_directional",
            "warren_i_directional",
            "warrens_i",
            "warren_i",
        ],
    }

    with h5py.File(h5_fp, "r") as h5:
        if "species" not in h5:
            raise KeyError("Missing 'species' dataset in H5 file")

        species = [
            s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)
            for s in h5["species"][:]
        ]

        resolved_names: dict[str, str] = {}
        metrics: dict[str, np.ndarray] = {}

        for metric_name, candidates in metric_aliases.items():
            ds_name = resolve_metric_dataset(h5, candidates)
            arr = np.array(h5[ds_name][:], dtype=float)
            if arr.shape != (len(species), len(species)):
                raise ValueError(
                    f"Dataset '{ds_name}' for metric '{metric_name}' has shape {arr.shape}, "
                    f"expected {(len(species), len(species))}"
                )
            resolved_names[metric_name] = ds_name
            metrics[metric_name] = arr

    return species, metrics, resolved_names


def build_species_metric_tables(
    species: list[str],
    metrics: dict[str, np.ndarray],
) -> tuple[dict[str, pd.Series], dict[str, pd.Series]]:
    outgoing: dict[str, pd.Series] = {}
    incoming: dict[str, pd.Series] = {}

    for metric_name, mat in metrics.items():
        outgoing[metric_name] = pd.Series(
            np.nanmean(mat, axis=1),
            index=species,
            name=metric_name,
        )
        incoming[metric_name] = pd.Series(
            np.nanmean(mat, axis=0),
            index=species,
            name=metric_name,
        )

    return outgoing, incoming


def matrix_from_series_dict(
    series_dict: dict[str, pd.Series],
    mode: str,
    top_n: int = 20,
) -> pd.DataFrame:
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
    species, metrics, resolved_names = load_metrics_from_h5(MATRIX_FP)

    outgoing, incoming = build_species_metric_tables(species, metrics)

    donor_top20_overlap = matrix_from_series_dict(outgoing, mode="overlap", top_n=20).astype(int)
    acceptor_top20_overlap = matrix_from_series_dict(incoming, mode="overlap", top_n=20).astype(int)
    spearman_outgoing = matrix_from_series_dict(outgoing, mode="spearman")
    spearman_incoming = matrix_from_series_dict(incoming, mode="spearman")

    donor_top20_overlap_fp = OUTDIR / "metric_sensitivity_donor_top20_overlap.csv"
    acceptor_top20_overlap_fp = OUTDIR / "metric_sensitivity_acceptor_top20_overlap.csv"
    spearman_outgoing_fp = OUTDIR / "metric_sensitivity_spearman_outgoing.csv"
    spearman_incoming_fp = OUTDIR / "metric_sensitivity_spearman_incoming.csv"
    summary_fp = OUTDIR / "metric_sensitivity_summary.json"

    donor_top20_overlap.to_csv(donor_top20_overlap_fp)
    acceptor_top20_overlap.to_csv(acceptor_top20_overlap_fp)
    spearman_outgoing.to_csv(spearman_outgoing_fp, float_format="%.6f")
    spearman_incoming.to_csv(spearman_incoming_fp, float_format="%.6f")

    summary = {
        "input_h5": str(MATRIX_FP),
        "n_species": int(len(species)),
        "resolved_metric_datasets": resolved_names,
        "top_n": int(20),
        "donor_top20_overlap_offdiag_mean": float(offdiag_mean(donor_top20_overlap)),
        "acceptor_top20_overlap_offdiag_mean": float(offdiag_mean(acceptor_top20_overlap)),
        "spearman_outgoing_offdiag_mean": float(offdiag_mean(spearman_outgoing)),
        "spearman_incoming_offdiag_mean": float(offdiag_mean(spearman_incoming)),
        "donor_top20_overlap_diagonal": int(donor_top20_overlap.iloc[0, 0]),
        "acceptor_top20_overlap_diagonal": int(acceptor_top20_overlap.iloc[0, 0]),
    }
    summary_fp.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {donor_top20_overlap_fp}")
    print(f"  {acceptor_top20_overlap_fp}")
    print(f"  {spearman_outgoing_fp}")
    print(f"  {spearman_incoming_fp}")
    print(f"  {summary_fp}")

    print("\nResolved datasets:")
    for k, v in resolved_names.items():
        print(f"  {k}: {v}")

    print("\nQuick summary:")
    print(f"  n_species: {len(species)}")
    print(f"  donor top20 overlap offdiag mean: {summary['donor_top20_overlap_offdiag_mean']:.3f}")
    print(f"  acceptor top20 overlap offdiag mean: {summary['acceptor_top20_overlap_offdiag_mean']:.3f}")
    print(f"  outgoing Spearman offdiag mean: {summary['spearman_outgoing_offdiag_mean']:.3f}")
    print(f"  incoming Spearman offdiag mean: {summary['spearman_incoming_offdiag_mean']:.3f}")


if __name__ == "__main__":
    main()
