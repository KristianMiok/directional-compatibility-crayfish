"""Stage 3b — Directional projection batch runner.

For each source species A with a trained SDM passing QC, project A's model onto
every other target species B's native occurrence environmental vectors.
Writes per-(source_species, predictor_set, core) parquet files with long-form
predictions. Resume-on-restart via output parquet existence.

Per-species predictor_set selection:
  - Default: full_geofresh (primary) + climate_local (piece 2 comparison)
  - Cambarus carolinus: local_full (fallback) + climate_local
  - Cherax quadricarinatus: EXCLUDED from matrix (insufficient native N)
  - Procambarus virginalis: EXCLUDED (from Stage 3a)

Cores retained per P2 + QC:
  - GBM: primary, retained where qc_pass
  - GLM: only on climate_local (and local_full for C. carolinus), retained where qc_pass

Target env vectors = presence rows (resp==1) from each species' biomod_input CSV
                     for the PRIMARY predictor set (full_geofresh, or local_full
                     for C. carolinus). This is B's "realized niche env space"
                     and is fixed across source species and cores so comparisons
                     are apples-to-apples.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE3_DIR = PROJECT_ROOT / "data" / "processed" / "stage3"
WORKSPACE = STAGE3_DIR / "biomod_workspace"
BIOMOD_INPUTS = STAGE3_DIR / "biomod_inputs"
ARTIFACTS = STAGE3_DIR / "sdm_artifacts"
PROJECTIONS_DIR = STAGE3_DIR / "projections"
TARGET_ENV_DIR = STAGE3_DIR / "target_env_vectors"
PROGRESS_LOG = STAGE3_DIR / "batch_progress" / "stage3b_progress.jsonl"

EXCLUDED_SPECIES = {"Procambarus virginalis", "Cherax quadricarinatus"}
CAMBARUS_CAROLINUS_PRIMARY = "local_full"  # fallback from full_geofresh


def slugify(name: str) -> str:
    """Same regex as R/02_run_stage3_production.R for cross-language consistency."""
    return re.sub(r"^_+|_+$", "", re.sub(r"[^a-z0-9]+", "_", name.lower()))


def log_event(event: dict) -> None:
    PROGRESS_LOG.parent.mkdir(parents=True, exist_ok=True)
    event["ts"] = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    with PROGRESS_LOG.open("a") as f:
        f.write(json.dumps(event) + "\n")


def primary_predictor_set(species: str) -> str:
    """Which predictor set defines B's 'realized niche env space' for this species."""
    if species == "Cambarus carolinus":
        return "local_full"
    return "full_geofresh"


def load_metadata(slug: str, predictor_set: str, core: str) -> dict | None:
    path = ARTIFACTS / slug / f"{predictor_set}__{core.lower()}_metadata.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def species_passes_qc(species: str, predictor_set: str, core: str) -> bool:
    meta = load_metadata(slugify(species), predictor_set, core)
    if meta is None:
        return False
    return bool(meta.get("qc_pass", False))


def build_target_env_csv(cohort: pd.DataFrame) -> Path:
    """For every species B in the cohort, extract presence env vectors from its
    primary-predictor biomod_input CSV. Concatenate into one long-form CSV used
    by every projection.

    Using only the primary predictor's columns means a source species trained on
    climate_local can still project onto B's primary-set env vectors — because
    climate_local's 6 predictors are a subset of full_geofresh's 16 and local_full's 8.
    """
    out_path = TARGET_ENV_DIR / "all_targets_env.csv"
    TARGET_ENV_DIR.mkdir(parents=True, exist_ok=True)

    frames = []
    for _, row in cohort.iterrows():
        sp = row["species_name"]
        ps = primary_predictor_set(sp)
        slug = slugify(sp)
        inp = BIOMOD_INPUTS / f"{slug}__{ps}_biomod_pilot.csv"
        if not inp.exists():
            print(f"  WARN: missing biomod_input for {sp} ({ps}): {inp}")
            continue
        df = pd.read_csv(inp)
        presences = df[df["resp"] == 1].copy()
        if presences.empty:
            continue
        presences["target_species"] = sp
        presences["target_cell_idx"] = range(len(presences))
        frames.append(presences)

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_path, index=False)
    print(f"Wrote target env vectors: {out_path} ({len(combined)} rows, "
          f"{combined['target_species'].nunique()} species)")
    return out_path


def run_projection(species: str, predictor_set: str, core: str,
                   target_env_csv: Path) -> tuple[bool, str]:
    slug = slugify(species)
    out_parquet = PROJECTIONS_DIR / f"{slug}__{predictor_set}__{core.lower()}.parquet"

    cmd = [
        "Rscript", "R/03_project_sdm.R",
        "--source-species", species,
        "--source-slug", slug,
        "--predictor-set", predictor_set,
        "--core", core,
        "--workspace", str(WORKSPACE.absolute()),
        "--target-env-csv", str(target_env_csv.absolute()),
        "--output-parquet", str(out_parquet.absolute()),
    ]

    start = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                                cwd=PROJECT_ROOT)
        dur = time.time() - start
        ok = result.returncode == 0 and out_parquet.exists()
        tail = (result.stdout + result.stderr)[-1500:]
        return ok, tail
    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT after {time.time()-start:.0f}s"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", default=None,
                        help="Subset of source species to project. Default = all valid.")
    parser.add_argument("--rebuild-target-env", action="store_true",
                        help="Rebuild the combined target env CSV even if it exists.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    # Load the locked cohort and apply exclusions
    cohort = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "cohort_threshold30.csv")
    cohort = cohort[~cohort["species_name"].isin(EXCLUDED_SPECIES)].copy()
    print(f"Matrix cohort: {len(cohort)} species (excluded: {sorted(EXCLUDED_SPECIES)})")

    # Build the combined target env CSV (one-time, reused for every projection)
    target_env = TARGET_ENV_DIR / "all_targets_env.csv"
    if args.rebuild_target_env or not target_env.exists():
        print("\nBuilding combined target env vectors CSV...")
        target_env = build_target_env_csv(cohort)
    else:
        print(f"\nReusing existing target env CSV: {target_env}")

    # Work out every (species, predictor_set, core) combo that needs projecting
    source_species = args.species if args.species else cohort["species_name"].tolist()
    source_species = [s for s in source_species if s not in EXCLUDED_SPECIES]

    combos: list[dict] = []
    for sp in source_species:
        slug = slugify(sp)
        # Determine valid predictor sets per species rules
        if sp == "Cambarus carolinus":
            sets = ["local_full", "climate_local"]
        else:
            sets = ["full_geofresh", "climate_local"]

        for ps in sets:
            # GBM always
            for core in ["GBM", "GLM"]:
                if core == "GLM" and ps not in ("climate_local", "local_full"):
                    continue
                if species_passes_qc(sp, ps, core):
                    combos.append({"species": sp, "slug": slug,
                                   "predictor_set": ps, "core": core})

    print(f"\nTotal projection jobs: {len(combos)}")

    PROJECTIONS_DIR.mkdir(parents=True, exist_ok=True)

    n_ok = n_skip = n_fail = 0
    failures = []

    for i, c in enumerate(combos, 1):
        out = PROJECTIONS_DIR / f"{c['slug']}__{c['predictor_set']}__{c['core'].lower()}.parquet"
        tag = f"{c['species']:40s} {c['predictor_set']:15s} {c['core']:3s}"

        if out.exists():
            n_skip += 1
            print(f"[{i:3d}/{len(combos)}] {tag}  SKIP")
            continue

        if args.dry_run:
            print(f"[{i:3d}/{len(combos)}] {tag}  DRY")
            continue

        print(f"[{i:3d}/{len(combos)}] {tag}  ...", end="", flush=True)
        ok, tail = run_projection(c["species"], c["predictor_set"], c["core"], target_env)
        if ok:
            n_ok += 1
            print("  OK")
            log_event({"kind": "projection", **c, "ok": True})
        else:
            n_fail += 1
            failures.append(c)
            print("  FAIL")
            print(tail[-600:])
            log_event({"kind": "projection", **c, "ok": False, "tail": tail[-1500:]})

    print("\n" + "=" * 70)
    print(f"Projections: {n_ok} ok, {n_skip} skipped, {n_fail} failed")
    if failures:
        print(f"\n{len(failures)} failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())