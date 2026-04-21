"""Stage 3a production SDM batch runner.

For each species in the locked >=30 cohort (minus P. virginalis):
  - Generate biomod-ready CSVs for the relevant predictor sets
  - Train ESM-GBM (always) and ESM-GLM (where Lucian's P2 rule allows)
  - Persist trained models to data/processed/stage3/biomod_workspace/
  - Write per-(species, predictor_set, core) metadata JSON

After all species complete, aggregate metadata into a per-species quality CSV
for Lucian's soft checkpoint review.

Resume-on-restart: skips any (species, predictor_set, core) for which the
metadata JSON already exists with fit_status == "ok".

Usage:
    python scripts/20_run_stage3_production_sdms.py
    python scripts/20_run_stage3_production_sdms.py --species "Astacus astacus" "Procambarus clarkii"
    python scripts/20_run_stage3_production_sdms.py --dry-run
    python scripts/20_run_stage3_production_sdms.py --aggregate-only
"""

from __future__ import annotations

import re

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
COHORT_CSV = PROJECT_ROOT / "data" / "processed" / "cohort_threshold30.csv"
STAGE3_DIR = PROJECT_ROOT / "data" / "processed" / "stage3"
BIOMOD_INPUTS_DIR = STAGE3_DIR / "biomod_inputs"
BIOMOD_WORKSPACE = STAGE3_DIR / "biomod_workspace"
ARTIFACTS_DIR = STAGE3_DIR / "sdm_artifacts"
QC_DIR = STAGE3_DIR / "qc"
PROGRESS_LOG = STAGE3_DIR / "batch_progress" / "stage3a_sdm_progress.jsonl"

EXCLUDED_SPECIES = {"Procambarus virginalis"}

CHERAX_EXCEPTION = {
    "Cherax quadricarinatus": ["climate_local", "local_full"],
}

DEFAULT_PREDICTOR_SETS = ["full_geofresh", "climate_local"]
GLM_ALLOWED_SETS = {"climate_local", "local_full"}


def slugify(name: str) -> str:
    """Normalize species name to filesystem-safe slug.

    Same regex as R/02_run_stage3_production.R so Python and R produce
    identical slugs (avoids file-not-found bugs across the two languages).
    Collapses any run of non-alphanumerics to a single underscore and trims
    leading/trailing underscores.
    Examples:
      "Cambarellus (Pandicambarus) puer" -> "cambarellus_puer"
      "Astacus astacus"                  -> "astacus_astacus"
    """
    return re.sub(r"^_+|_+$", "", re.sub(r"[^a-z0-9]+", "_", name.lower()))


def predictor_sets_for_species(species: str) -> list[str]:
    return CHERAX_EXCEPTION.get(species, DEFAULT_PREDICTOR_SETS)


def cores_for_predictor_set(predictor_set: str) -> list[str]:
    if predictor_set in GLM_ALLOWED_SETS:
        return ["GBM", "GLM"]
    return ["GBM"]


def log_event(event: dict) -> None:
    PROGRESS_LOG.parent.mkdir(parents=True, exist_ok=True)
    event["ts"] = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    with PROGRESS_LOG.open("a") as f:
        f.write(json.dumps(event) + "\n")


def run_cmd(cmd: list[str], timeout: int = 1800) -> tuple[bool, str]:
    start = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        ok = result.returncode == 0
        output = (result.stdout + result.stderr)[-2000:]
        return ok, output
    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT after {time.time()-start:.0f}s"
    except Exception as e:
        return False, f"EXCEPTION: {e}"


def export_path(slug: str, predictor_set: str) -> Path:
    return BIOMOD_INPUTS_DIR / f"{slug}__{predictor_set}_biomod_pilot.csv"


def metadata_path(slug: str, predictor_set: str, core: str) -> Path:
    return ARTIFACTS_DIR / slug / f"{predictor_set}__{core.lower()}_metadata.json"


def metadata_is_ok(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        meta = json.loads(path.read_text())
        return meta.get("fit_status") == "ok"
    except Exception:
        return False


def aggregate_quality_report() -> Path:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for json_file in sorted(ARTIFACTS_DIR.rglob("*_metadata.json")):
        try:
            meta = json.loads(json_file.read_text())
        except Exception:
            continue
        rows.append({
            "species": meta.get("species"),
            "species_slug": meta.get("species_slug"),
            "predictor_set": meta.get("predictor_set"),
            "core": meta.get("core"),
            "fit_status": meta.get("fit_status"),
            "n_records": meta.get("n_records"),
            "n_presences": meta.get("n_presences"),
            "n_cv_runs": meta.get("n_cv_runs"),
            "auc_mean": meta.get("auc_mean"),
            "auc_std": meta.get("auc_std"),
            "tss_mean": meta.get("tss_mean"),
            "tss_std": meta.get("tss_std"),
            "boyce_mean": meta.get("boyce_mean"),
            "boyce_std": meta.get("boyce_std"),
            "boyce_n_valid": meta.get("boyce_n_valid"),
            "maxsss_mean": meta.get("maxsss_mean"),
            "maxsss_std": meta.get("maxsss_std"),
            "qc_pass": meta.get("qc_pass"),
            "unreliable_boyce": meta.get("unreliable_boyce"),
        })
    df = pd.DataFrame(rows)
    out = QC_DIR / "per_species_quality.csv"
    df.to_csv(out, index=False)

    if len(df) and "qc_pass" in df.columns:
        excluded = df[df["qc_pass"] == False][["species", "predictor_set", "core",
                                                  "auc_mean", "tss_mean", "fit_status"]]
        excluded.to_csv(QC_DIR / "qc_exclusions.csv", index=False)
    if len(df) and "unreliable_boyce" in df.columns:
        unreliable = df[df["unreliable_boyce"] == True][["species", "predictor_set",
                                                            "core", "boyce_n_valid"]]
        unreliable.to_csv(QC_DIR / "unreliable_boyce_species.csv", index=False)

    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", default=None,
                        help="Subset of species to run (default = full cohort minus excluded)")
    parser.add_argument("--predictor-sets", nargs="+", default=None,
                        help="Restrict to subset (default = full_geofresh + climate_local per Lucian's P2)")
    parser.add_argument("--cores", nargs="+", default=None,
                        help="Restrict to subset (default determined per predictor_set)")
    parser.add_argument("--force-export", action="store_true")
    parser.add_argument("--force-run", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="Skip running anything, just aggregate existing metadata into QC CSV")
    args = parser.parse_args()

    if args.aggregate_only:
        out = aggregate_quality_report()
        print(f"Wrote {out}")
        return 0

    cohort = pd.read_csv(COHORT_CSV)
    cohort = cohort[~cohort["species_name"].isin(EXCLUDED_SPECIES)].copy()
    if args.species:
        cohort = cohort[cohort["species_name"].isin(args.species)]
        if cohort.empty:
            print(f"ERROR: none of {args.species} in cohort", file=sys.stderr)
            return 2

    species_list = cohort["species_name"].tolist()
    print(f"Cohort: {len(species_list)} species (excluded: {sorted(EXCLUDED_SPECIES)})")
    print(f"Cherax exception: {CHERAX_EXCEPTION}")
    print(f"GLM allowed on: {sorted(GLM_ALLOWED_SETS)}")
    print()

    BIOMOD_INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    BIOMOD_WORKSPACE.mkdir(parents=True, exist_ok=True)
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    n_export_skip = n_export_ok = n_export_fail = 0
    n_run_skip = n_run_ok = n_run_fail = 0
    failures: list[dict] = []

    for i, species in enumerate(species_list, 1):
        slug = slugify(species)
        ps_list = predictor_sets_for_species(species)
        if args.predictor_sets:
            ps_list = [p for p in ps_list if p in args.predictor_sets]

        print(f"\n[{i}/{len(species_list)}] {species}  ({len(ps_list)} predictor sets)")

        for ps in ps_list:
            ex_path = export_path(slug, ps)
            if ex_path.exists() and not args.force_export:
                n_export_skip += 1
                print(f"  export {ps:14s} SKIP")
            else:
                cmd = ["python", "scripts/10_export_stage2_pilot_species.py",
                       "--species", species, "--predictor-set", ps,
                       "--output-dir", str(BIOMOD_INPUTS_DIR),
                       "--output-slug", slug]
                if args.dry_run:
                    print(f"  export {ps:14s} DRY  {' '.join(cmd)}")
                    n_export_ok += 1
                else:
                    ok, tail = run_cmd(cmd)
                    log_event({"kind": "export", "species": species, "predictor_set": ps,
                               "ok": ok, "tail": tail if not ok else None})
                    if ok:
                        n_export_ok += 1
                        print(f"  export {ps:14s} OK")
                    else:
                        n_export_fail += 1
                        failures.append({"species": species, "predictor_set": ps, "stage": "export"})
                        print(f"  export {ps:14s} FAIL\n{tail[-400:]}")
                        continue

            cores = cores_for_predictor_set(ps)
            if args.cores:
                cores = [c for c in cores if c in args.cores]

            for core in cores:
                meta_path = metadata_path(slug, ps, core)
                if metadata_is_ok(meta_path) and not args.force_run:
                    n_run_skip += 1
                    print(f"    run  {ps:14s} {core:3s}  SKIP")
                    continue

                cmd = ["Rscript", "R/02_run_stage3_production.R",
                       "--species", species,
                       "--predictor-set", ps,
                       "--core", core,
                       "--input-dir", str(BIOMOD_INPUTS_DIR),
                       "--workspace", str(BIOMOD_WORKSPACE),
                       "--artifacts-dir", str(ARTIFACTS_DIR)]
                if args.dry_run:
                    print(f"    run  {ps:14s} {core:3s}  DRY  {' '.join(cmd)}")
                    n_run_ok += 1
                else:
                    ok, tail = run_cmd(cmd, timeout=3600)
                    log_event({"kind": "run", "species": species, "predictor_set": ps,
                               "core": core, "ok": ok, "tail": tail if not ok else None})
                    if ok:
                        n_run_ok += 1
                        print(f"    run  {ps:14s} {core:3s}  OK")
                    else:
                        n_run_fail += 1
                        failures.append({"species": species, "predictor_set": ps,
                                         "core": core, "stage": "run"})
                        print(f"    run  {ps:14s} {core:3s}  FAIL")
                        print(tail[-500:])

    print()
    print("=" * 70)
    print(f"Exports: {n_export_ok} ok, {n_export_skip} skipped, {n_export_fail} failed")
    print(f"Runs:    {n_run_ok} ok, {n_run_skip} skipped, {n_run_fail} failed")
    if failures:
        print(f"\n{len(failures)} failures:")
        for f in failures:
            print(f"  {f}")

    if not args.dry_run:
        out = aggregate_quality_report()
        print(f"\nWrote QC report: {out}")

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())