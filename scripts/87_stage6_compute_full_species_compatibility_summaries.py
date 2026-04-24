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
OUTDIR = ROOT / "data" / "processed" / "stage6_invasion_followup"
OUTDIR.mkdir(parents=True, exist_ok=True)

MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"

OUT_CSV = OUTDIR / "species_full_compatibility_summaries.csv"
OUT_JSON = OUTDIR / "species_full_compatibility_summaries_summary.json"


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
        f"Available datasets:\n" + "\n".join(all_keys)
    )


def load_matrix_metrics(h5_fp: Path) -> tuple[list[str], dict[str, np.ndarray], dict[str, str]]:
    metric_aliases = {
        "mean_suitability": ["mean_suitability"],
        "fraction_above": ["fraction_above", "fraction_above_threshold", "fraction_above_thr", "fraction_above_0.5"],
        "schoener_D": ["schoener_D", "directional_schoeners_d", "directional_schoener_d", "schoeners_d"],
        "warren_I": ["warren_I", "directional_warrens_i", "directional_warren_i", "warrens_i"],
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
                    f"Dataset '{ds_name}' shape {arr.shape} != expected {(len(species), len(species))}"
                )
            metrics[metric_name] = arr
            resolved[metric_name] = ds_name

    return species, metrics, resolved


def species_slug(species_name: str) -> str:
    return (
        species_name.lower()
        .replace("(", "")
        .replace(")", "")
        .replace(",", "")
        .replace(".", "")
        .replace("-", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )


def load_q90_matrix(species: list[str], proj_dir: Path) -> pd.DataFrame:
    rows = []
    for sp in species:
        fp = proj_dir / f"{species_slug(sp)}__full_geofresh__gbm.parquet"
        if not fp.exists():
            raise FileNotFoundError(f"Missing projection parquet: {fp}")

        df = pd.read_parquet(fp)
        grp = df.groupby("target_species")["predicted_suitability"].quantile(0.90)
        rows.append([float(grp.get(target, np.nan)) for target in species])

    out = pd.DataFrame(rows, index=species, columns=species)
    return out


def rank_desc(series: pd.Series) -> pd.Series:
    return series.rank(method="average", ascending=False)


def main() -> None:
    species, matrix_metrics, resolved = load_matrix_metrics(MATRIX_FP)

    print("Loading q90 matrix from projection parquets...")
    q90_mat = load_q90_matrix(species, PROJDIR)

    out = pd.DataFrame({"species": species})

    for metric_name, mat in matrix_metrics.items():
        outgoing = pd.Series(np.nanmean(mat, axis=1), index=species)
        incoming = pd.Series(np.nanmean(mat, axis=0), index=species)

        out[f"mean_outgoing_{metric_name}"] = out["species"].map(outgoing)
        out[f"mean_incoming_{metric_name}"] = out["species"].map(incoming)
        out[f"rank_outgoing_{metric_name}"] = out["species"].map(rank_desc(outgoing))
        out[f"rank_incoming_{metric_name}"] = out["species"].map(rank_desc(incoming))

    q90_outgoing = pd.Series(np.nanmean(q90_mat.to_numpy(dtype=float), axis=1), index=species)
    q90_incoming = pd.Series(np.nanmean(q90_mat.to_numpy(dtype=float), axis=0), index=species)

    out["mean_outgoing_q90_suitability"] = out["species"].map(q90_outgoing)
    out["mean_incoming_q90_suitability"] = out["species"].map(q90_incoming)
    out["rank_outgoing_q90_suitability"] = out["species"].map(rank_desc(q90_outgoing))
    out["rank_incoming_q90_suitability"] = out["species"].map(rank_desc(q90_incoming))

    out.to_csv(OUT_CSV, index=False)

    summary = {
        "n_species": int(len(out)),
        "resolved_metric_datasets": resolved,
        "columns": list(out.columns),
        "output_csv": str(OUT_CSV),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_CSV}")
    print(f"  {OUT_JSON}")

    print("\nFirst 10 rows:")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
