"""Stage 2 pilot batch runner with resume-on-restart.

For each (species, predictor_set) pair, calls the Python exporter if the
biomod-ready CSV doesn't exist yet. For each (species, predictor_set, core),
calls the R ESM runner unless the output evaluations CSV already exists.

Progress is logged line-by-line to a JSONL file, so crashed runs don't lose
history and a restart resumes cleanly.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PILOT_DIR = PROJECT_ROOT / "data" / "processed" / "stage2_pilot"
PROGRESS_LOG = PILOT_DIR / "batch_progress.jsonl"
SHORTLIST = PROJECT_ROOT / "data" / "processed" / "stage2_pilot_shortlist_frozen.csv"

DEFAULT_PREDICTOR_SETS = ["climate_local", "local_full", "upstream_full", "full_geofresh"]
DEFAULT_CORES = ["GLM", "GBM"]


def slugify(name: str) -> str:
    return name.lower().replace(" ", "_").replace("/", "_")


def log_event(event: dict) -> None:
    PILOT_DIR.mkdir(parents=True, exist_ok=True)
    event["ts"] = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    with PROGRESS_LOG.open("a") as f:
        f.write(json.dumps(event) + "\n")


def run_cmd(cmd: list[str], label: str) -> tuple[bool, str]:
    """Run a command, capturing stdout+stderr. Returns (ok, tail_of_output)."""
    start = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        dur = time.time() - start
        ok = result.returncode == 0
        output = (result.stdout + result.stderr)[-2000:]
        return ok, output
    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT after {time.time()-start:.0f}s"
    except Exception as e:
        return False, f"EXCEPTION: {e}"


def export_exists(slug: str, predictor_set: str) -> Path | None:
    p = PILOT_DIR / f"{slug}__{predictor_set}_biomod_pilot.csv"
    return p if p.exists() else None


def eval_exists(slug: str, predictor_set: str, core: str) -> Path | None:
    p = PILOT_DIR / f"{slug}__{predictor_set}_esm_{core.lower()}_evaluations.csv"
    return p if p.exists() else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", default=None,
                        help="Subset of species to run. Default = all from frozen shortlist.")
    parser.add_argument("--predictor-sets", nargs="+", default=DEFAULT_PREDICTOR_SETS)
    parser.add_argument("--cores", nargs="+", default=DEFAULT_CORES)
    parser.add_argument("--force-export", action="store_true",
                        help="Re-export biomod CSVs even if they exist")
    parser.add_argument("--force-run", action="store_true",
                        help="Re-run ESM even if evaluations CSV exists")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would run, but don't execute")
    args = parser.parse_args()

    shortlist = pd.read_csv(SHORTLIST)
    if args.species:
        shortlist = shortlist[shortlist["species_name"].isin(args.species)]
        if shortlist.empty:
            print(f"ERROR: none of {args.species} found in frozen shortlist", file=sys.stderr)
            return 2

    species_list = shortlist["species_name"].tolist()
    total_exports = len(species_list) * len(args.predictor_sets)
    total_runs = total_exports * len(args.cores)

    print(f"Frozen shortlist: {len(species_list)} species")
    print(f"Predictor sets: {args.predictor_sets}")
    print(f"Cores: {args.cores}")
    print(f"Total: {total_exports} exports + {total_runs} ESM runs")
    print()

    n_export_skipped = n_export_ok = n_export_fail = 0
    n_run_skipped = n_run_ok = n_run_fail = 0
    failures: list[dict] = []

    for i, species in enumerate(species_list, 1):
        slug = slugify(species)
        print(f"\n[{i}/{len(species_list)}] {species}")

        for ps in args.predictor_sets:
            existing = export_exists(slug, ps)
            if existing and not args.force_export:
                n_export_skipped += 1
                print(f"  export {ps:16s} SKIP  (exists: {existing.name})")
            else:
                cmd = ["python", "scripts/10_export_stage2_pilot_species.py",
                       "--species", species, "--predictor-set", ps]
                if args.dry_run:
                    print(f"  export {ps:16s} DRY   {' '.join(cmd)}")
                    n_export_ok += 1
                else:
                    ok, tail = run_cmd(cmd, f"export {slug} {ps}")
                    log_event({"kind": "export", "species": species, "predictor_set": ps,
                               "ok": ok, "tail": tail if not ok else None})
                    if ok:
                        n_export_ok += 1
                        print(f"  export {ps:16s} OK")
                    else:
                        n_export_fail += 1
                        failures.append({"species": species, "predictor_set": ps, "stage": "export"})
                        print(f"  export {ps:16s} FAIL\n{tail}")
                        continue

            for core in args.cores:
                existing_eval = eval_exists(slug, ps, core)
                if existing_eval and not args.force_run:
                    n_run_skipped += 1
                    print(f"    run  {ps:16s} {core:3s}  SKIP  (exists)")
                    continue

                cmd = ["Rscript", "R/01_run_stage2_esm_single_core.R",
                       "--species", species, "--predictor-set", ps, "--core", core]
                if args.dry_run:
                    print(f"    run  {ps:16s} {core:3s}  DRY   {' '.join(cmd)}")
                    n_run_ok += 1
                else:
                    ok, tail = run_cmd(cmd, f"run {slug} {ps} {core}")
                    log_event({"kind": "run", "species": species, "predictor_set": ps,
                               "core": core, "ok": ok, "tail": tail if not ok else None})
                    if ok:
                        n_run_ok += 1
                        print(f"    run  {ps:16s} {core:3s}  OK")
                    else:
                        n_run_fail += 1
                        failures.append({"species": species, "predictor_set": ps,
                                         "core": core, "stage": "run"})
                        print(f"    run  {ps:16s} {core:3s}  FAIL")
                        print(tail[-400:])

    print()
    print("=" * 70)
    print(f"Exports: {n_export_ok} ok, {n_export_skipped} skipped, {n_export_fail} failed")
    print(f"Runs:    {n_run_ok} ok, {n_run_skipped} skipped, {n_run_fail} failed")
    if failures:
        print(f"\n{len(failures)} failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
