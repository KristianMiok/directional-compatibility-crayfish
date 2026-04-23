#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import json

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
BOOTROOT = STAGE3 / "stage4d"

R_FIT = ROOT / "R" / "02_run_stage3_production.R"
R_PROJECT = ROOT / "R" / "03_project_sdm.R"
PRIMARY_MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"


def load_primary_species() -> list[str]:
    with h5py.File(PRIMARY_MATRIX_FP, "r") as h5:
        return [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]


def run_cmd(cmd: list[str], log_path: Path | None = None) -> None:
    if log_path is None:
        print("RUN:", " ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)


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


def species_name_from_input_csv(fp: Path) -> str:
    df = pd.read_csv(fp, nrows=3)
    if "Crayfish_scientific_name" in df.columns:
        vals = df["Crayfish_scientific_name"].dropna().unique().tolist()
        if len(vals) >= 1:
            return str(vals[0])
    return fp.stem.split("__")[0].replace("_", " ")


def build_target_env_csv(iter_dir: Path, predictor_set: str) -> Path:
    in_dir = iter_dir / "biomod_inputs"
    out_dir = iter_dir / "target_env_vectors"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fp = out_dir / f"all_targets_env__{predictor_set}.csv"

    rows = []
    csvs = sorted(in_dir.glob(f"*__{predictor_set}_biomod_pilot.csv"))

    target_cell_idx = 0
    for fp in csvs:
        df = pd.read_csv(fp)
        if "resp" not in df.columns:
            raise ValueError(f"Missing resp column in {fp}")

        pres = df[df["resp"] == 1].copy()
        species_name = species_name_from_input_csv(fp)

        pres["target_species"] = species_name
        pres["target_cell_idx"] = np.arange(target_cell_idx, target_cell_idx + len(pres))
        target_cell_idx += len(pres)

        if "subc_id" in pres.columns:
            pres["target_cell_subc_id"] = pres["subc_id"]
        else:
            pres["target_cell_subc_id"] = np.arange(len(pres))

        rows.append(pres)

    if not rows:
        raise RuntimeError(f"No bootstrap biomod inputs found for predictor_set={predictor_set}")

    target_df = pd.concat(rows, axis=0).reset_index(drop=True)
    target_df.to_csv(out_fp, index=False)
    return out_fp


def fit_models(
    iter_dir: Path,
    predictor_set: str,
    core: str,
    species_subset: list[str] | None = None,
) -> list[str]:
    input_dir = iter_dir / "biomod_inputs"
    workspace = iter_dir / "biomod_workspace"
    artifacts = iter_dir / "sdm_artifacts"
    log_dir = iter_dir / "logs" / "fit"

    workspace.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    csvs = sorted(input_dir.glob(f"*__{predictor_set}_biomod_pilot.csv"))

    species_all = []
    for fp in csvs:
        species_name = species_name_from_input_csv(fp)
        if species_subset is not None and species_name not in species_subset:
            continue
        species_all.append(species_name)

    species_run = []
    total = len(species_all)

    for i, species_name in enumerate(species_all, start=1):
        species_run.append(species_name)
        print(f"[FIT {i}/{total}] {species_name}")

        cmd = [
            "Rscript",
            str(R_FIT),
            "--species", species_name,
            "--predictor-set", predictor_set,
            "--core", core,
            "--input-dir", str(input_dir),
            "--workspace", str(workspace),
            "--artifacts-dir", str(artifacts),
        ]
        run_cmd(cmd, log_dir / f"{slugify_species(species_name)}.log")

    return species_run


def project_models(
    iter_dir: Path,
    predictor_set: str,
    core: str,
    species_list: list[str],
    target_env_csv: Path,
) -> None:
    workspace = iter_dir / "biomod_workspace"
    proj_dir = iter_dir / "projections"
    log_dir = iter_dir / "logs" / "project"

    proj_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    total = len(species_list)

    for i, species_name in enumerate(species_list, start=1):
        source_slug = slugify_species(species_name)
        out_fp = proj_dir / f"{source_slug}__{predictor_set}__{core.lower()}.parquet"

        print(f"[PROJECT {i}/{total}] {species_name}")

        cmd = [
            "Rscript",
            str(R_PROJECT),
            "--source-species", species_name,
            "--source-slug", source_slug,
            "--predictor-set", predictor_set,
            "--target-env-csv", str(target_env_csv),
            "--output-parquet", str(out_fp),
            "--workspace", str(workspace),
        ]
        run_cmd(cmd, log_dir / f"{source_slug}.log")


def aggregate_mean_matrix(
    iter_dir: Path,
    predictor_set: str,
    core: str,
    species_order: list[str],
) -> Path:
    proj_dir = iter_dir / "projections"
    out_fp = iter_dir / f"matrices_{predictor_set}_{core.lower()}.h5"

    rows = []
    for species_name in species_order:
        slug = slugify_species(species_name)
        pq_fp = proj_dir / f"{slug}__{predictor_set}__{core.lower()}.parquet"
        if not pq_fp.exists():
            raise FileNotFoundError(f"Missing projection parquet: {pq_fp}")

        df = pd.read_parquet(pq_fp)

        if "target_species" not in df.columns:
            raise ValueError(f"target_species column missing in {pq_fp}")
        if "predicted_suitability" not in df.columns:
            raise ValueError(f"predicted_suitability column missing in {pq_fp}")

        grp = df.groupby("target_species")["predicted_suitability"].mean()
        row = [float(grp.get(sp, np.nan)) for sp in species_order]
        rows.append(row)

    mat = np.array(rows, dtype=np.float32)

    with h5py.File(out_fp, "w") as h5:
        h5.create_dataset("species", data=np.array(species_order, dtype="S"))
        h5.create_dataset("mean_suitability", data=mat)

    return out_fp


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iteration", required=True, help="e.g. iter_001")
    ap.add_argument("--predictor-set", default="full_geofresh")
    ap.add_argument("--core", default="GBM")
    ap.add_argument("--species", nargs="*", default=None, help="Optional subset for smoke test")
    args = ap.parse_args()

    iter_dir = BOOTROOT / args.iteration
    if not iter_dir.exists():
        raise FileNotFoundError(f"Iteration directory not found: {iter_dir}")

    primary_species = load_primary_species()
    species_subset = args.species if args.species else primary_species

    print("Building bootstrap target env CSV ...")
    target_env_csv = build_target_env_csv(iter_dir, args.predictor_set)
    print(f"Target env CSV: {target_env_csv}")

    print("\nFitting bootstrap SDMs ...")
    species_run = fit_models(
        iter_dir=iter_dir,
        predictor_set=args.predictor_set,
        core=args.core,
        species_subset=species_subset,
    )

    print("\nProjecting bootstrap SDMs ...")
    project_models(
        iter_dir=iter_dir,
        predictor_set=args.predictor_set,
        core=args.core,
        species_list=species_run,
        target_env_csv=target_env_csv,
    )

    if set(species_run) == set(primary_species):
        print("\nAggregating bootstrap mean_suitability matrix ...")
        out_fp = aggregate_mean_matrix(
            iter_dir=iter_dir,
            predictor_set=args.predictor_set,
            core=args.core,
            species_order=primary_species,
        )
        print(f"Matrix written: {out_fp}")
    else:
        out_fp = None
        print("\nSpecies subset run only; skipping final matrix aggregation.")

    summary = {
        "iteration": args.iteration,
        "predictor_set": args.predictor_set,
        "core": args.core,
        "n_species_run": len(species_run),
        "target_env_csv": str(target_env_csv),
        "matrix_output": str(out_fp) if out_fp is not None else None,
    }
    (iter_dir / "bootstrap_run_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {iter_dir / 'bootstrap_run_summary.json'}")


if __name__ == "__main__":
    main()