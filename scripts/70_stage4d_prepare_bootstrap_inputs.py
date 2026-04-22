#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import h5py


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
OUTROOT = STAGE3 / "stage4d"
INPUTS_SRC = STAGE3 / "biomod_inputs"
MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"

OUTROOT.mkdir(parents=True, exist_ok=True)


def load_primary_species() -> list[str]:
    with h5py.File(MATRIX_FP, "r") as h5:
        return [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]


def slugify_species(name: str) -> str:
    return (
        name.lower()
        .replace("(", "")
        .replace(")", "")
        .replace(",", "")
        .replace(".", "")
        .replace("-", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )


def species_from_filename(path: Path) -> str:
    # filenames look like:
    # cambarus_carolinus__local_full_biomod_pilot.csv
    stem = path.stem
    parts = stem.split("__")
    if len(parts) < 2:
        return stem
    return parts[0]


def bootstrap_single_input(
    df: pd.DataFrame,
    rng: np.random.Generator,
    presence_col: str = "resp",
    keep_frac: float = 0.9,
) -> tuple[pd.DataFrame, dict]:
    if presence_col not in df.columns:
        raise ValueError(f"Missing presence column: {presence_col}")

    pres = df[df[presence_col] == 1].copy()
    bg = df[df[presence_col] != 1].copy()

    n_pres = len(pres)
    if n_pres == 0:
        raise ValueError("No presences found.")

    n_keep = max(1, math.ceil(keep_frac * n_pres))
    keep_idx = rng.choice(pres.index.to_numpy(), size=n_keep, replace=False)
    pres_boot = pres.loc[sorted(keep_idx)].copy()

    out = pd.concat([pres_boot, bg], axis=0).reset_index(drop=True)

    meta = {
        "n_pres_original": int(n_pres),
        "n_pres_bootstrap": int(n_keep),
        "n_background": int(len(bg)),
        "dropped_presence_fraction": float(1.0 - n_keep / n_pres),
    }
    return out, meta


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--n-iter", type=int, default=3, help="Number of bootstrap iterations to prepare.")
    ap.add_argument("--keep-frac", type=float, default=0.9, help="Fraction of presence records to keep.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    primary_species = load_primary_species()
    primary_slugs = {slugify_species(s) for s in primary_species}

    csv_files = sorted(INPUTS_SRC.glob("*_biomod_pilot.csv"))
    if not csv_files:
        raise RuntimeError(f"No biomod input CSVs found in {INPUTS_SRC}")

    # Only keep files for the 155-species primary cohort plus any fallback local_full files that match those species
    selected = []
    for fp in csv_files:
        sp_slug = species_from_filename(fp)
        if sp_slug in primary_slugs:
            selected.append(fp)

    if not selected:
        raise RuntimeError("No input CSVs matched the primary matrix species set.")

    manifest_rows = []

    for it in range(1, args.n_iter + 1):
        iter_dir = OUTROOT / f"iter_{it:03d}" / "biomod_inputs"
        iter_dir.mkdir(parents=True, exist_ok=True)

        rng = np.random.default_rng(args.seed + it)

        print(f"Preparing bootstrap iteration {it}/{args.n_iter} ...")

        for fp in selected:
            df = pd.read_csv(fp)
            boot_df, meta = bootstrap_single_input(
                df=df,
                rng=rng,
                presence_col="resp",
                keep_frac=args.keep_frac,
            )

            out_fp = iter_dir / fp.name
            boot_df.to_csv(out_fp, index=False)

            manifest_rows.append({
                "iteration": it,
                "source_file": str(fp),
                "bootstrap_file": str(out_fp),
                **meta,
            })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUTROOT / "bootstrap_input_manifest.csv", index=False)

    summary = {
        "n_iter": int(args.n_iter),
        "keep_frac": float(args.keep_frac),
        "seed": int(args.seed),
        "n_selected_input_files": int(len(selected)),
        "manifest_file": str(OUTROOT / "bootstrap_input_manifest.csv"),
    }
    (OUTROOT / "bootstrap_input_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {OUTROOT / 'bootstrap_input_manifest.csv'}")
    print(f"  {OUTROOT / 'bootstrap_input_summary.json'}")

    print("\nFirst 10 manifest rows:")
    print(manifest.head(10).to_string(index=False))


if __name__ == "__main__":
    main()