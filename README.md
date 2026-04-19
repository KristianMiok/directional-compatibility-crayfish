# directional-compatibility-crayfish

Code and data pipeline for **Manuscript #5 — Directional Niche Compatibility**
(target journal: *Ecography*).

A directional ecological compatibility framework for freshwater crayfish: for
each species A we train an SDM on native records only, then project A onto the
native environmental space of every other species B, producing an asymmetric
compatibility matrix. The matrix is analysed as a macroecological object via
asymmetry, clustering, and network metrics.

> Projections are restricted to the **crayfish-occupied freshwater network
> envelope**. Global extrapolation is out of scope for this paper.

## Project status

Currently in **Stage 1 — Data Audit and Species Selection** (≈2 weeks).
See `docs/PLAN.md` for the full work plan and decision points (P1–P3).

## Repo layout

```
configs/         YAML configs (paths, thresholds, variable lists)
data/            Local data (gitignored)
  raw/           Untouched source dumps (WoC export, GeoFRESH pulls)
  interim/       Cleaned but not yet analysis-ready
  processed/     Analysis-ready tables (per-species, per-segment, etc.)
  external/      Reference layers (basin polygons, climate rasters)
docs/            Plan, decision-point notes, methods drafts
notebooks/       Exploratory work — disposable, not pipeline
reports/         Stage reports + figures (committed)
  stage1/        Stage 1 deliverables (CSVs go to data/processed; reports here)
scripts/         CLI entry points that orchestrate src/dcc modules
src/dcc/         Importable package (the real code lives here)
```

## Stage 1 deliverables

| Task | Output                                  | Script                            |
|------|-----------------------------------------|-----------------------------------|
| 1.1  | `species_inventory.csv`                 | `scripts/01_species_inventory.py` |
| 1.2  | `threshold_sensitivity.{csv,png}`       | `scripts/02_threshold_sensitivity.py` |
| 1.3  | `env_coverage.csv` + per-species plots  | `scripts/03_env_coverage.py`      |
| 1.4  | `basin_overlap_matrix.csv`              | `scripts/04_basin_overlap.py`     |
| 1.5  | `invasion_contamination_flagged.csv`    | `scripts/05_invasion_audit.py`    |
| —    | `etapa1_raport.md`                      | written by hand in `reports/stage1/` |

## Quickstart

```bash
# Environment
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install

# Drop data exports into data/raw/ (see configs/data.yaml for expected names)

# Run Stage 1 end-to-end
make stage1
```

## Reproducibility notes

- All paths and parameters live in `configs/` — no hard-coded values in `src/`.
- Intermediate data: CSV for tables, GeoJSON for geometries, HDF5 for large arrays.
- Pipeline is **frozen after Stage 2** (Decision Point P2). Until then, expect churn.

## Authors

- Lucian Pârvulescu — concept, manuscript lead
- Kristian Miok — implementation, modelling
