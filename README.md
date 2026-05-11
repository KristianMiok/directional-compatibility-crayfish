# directional-compatibility-crayfish

Code, supplementary tables, and pipeline for the manuscript

**"Directional ecological asymmetry in a globally distributed, network-constrained freshwater taxon"**
 — submitted to *Ecology* (Manuscript ID: ECY26-0558).

A directional ecological compatibility framework for freshwater crayfish.
For each species *A* in a 157-species global cohort we train a species
distribution model on native-only records, then project *A*'s realised niche
onto the native environmental space of every other species *B*, producing an
**asymmetric N×N compatibility matrix**. The matrix is analysed as a
macroecological object via asymmetry tests, hierarchical clustering, network
community detection, and bootstrap stability analysis.

Projections are restricted to the crayfish-occupied freshwater network
envelope. Global extrapolation is explicitly out of scope.

## Status

Manuscript submitted to *Ecology* (May 2026, in editorial screening).
The pipeline is frozen at the version corresponding to the submitted
manuscript. Post-acceptance changes will be tagged.

## Supplementary tables (publicly archived)

All tables referenced in the manuscript appendix are available under
[`reports/supplementary/`](reports/supplementary/):

| Table | Content                                                  |
|-------|----------------------------------------------------------|
| S1    | Stage 3 compact workbook (XLSX, per-species SDM summary) |
| S2    | Cluster membership (k=6)                                 |
| S3    | Cluster concordance metrics                              |
| S4    | Per-species network metrics (155 × 12)                   |
| S5    | Top-50 shifted pairs (directed)                          |
| S6    | Top-50 shifted pairs (undirected)                        |
| S7    | Top-50 repositioners                                     |
| S8    | Bootstrap iteration stability                            |
| S9    | Species top-20 bootstrap frequencies                     |
| S10   | Community-1 within-community matrix                      |
| S11   | Community-1 reach into Cambaridae clusters               |
| S12   | Community-1 cross-continental asymmetric pairs           |
| S13   | Metric sensitivity — donor top-20 overlap                |
| S14   | Metric sensitivity — acceptor top-20 overlap             |
| S15   | Metric sensitivity — Spearman outgoing                   |
| S16   | Metric sensitivity — Spearman incoming                   |
| S17   | Bootstrap structural stability                           |
| S18   | VIF retained per species                                 |

## Repository layout

```
config/, configs/    Predictor set definitions, paths, thresholds
data/                Local data (gitignored except cohort/shortlist anchors)
  raw/               Untouched source dumps (WoC export, GeoFRESH pulls)
  processed/         Analysis-ready tables, per-species artifacts,
                     compatibility matrices (HDF5), bootstrap outputs
docs/                Methods drafts and decision-point notes
notebooks/           Exploratory work — disposable, not pipeline
R/                   biomod2 production runner and ESM exploratory scripts
reports/             Stage reports, figures, and supplementary tables
  stage1/            Cohort selection and threshold sensitivity
  stage2/            ESM pilot report (calibration and predictor-set comparison)
  stage3/            Stage 3 plan and SDM production outputs
  stage5/            Final figures and tables for the manuscript
  supplementary/     Publicly archived supplementary tables (see above)
  figures/           Main-text figures
scripts/             CLI entry points for each pipeline stage
src/dcc/             Importable package (core logic)
```

## Pipeline stages

| Stage | Content                                                       | Status   |
|-------|---------------------------------------------------------------|----------|
| 1     | Cohort selection and quality audit (threshold N ≥ 30)         | Complete |
| 2     | SDM calibration pilot, predictor-set comparison, P2 decision  | Complete |
| 3     | Production SDMs (157 species), directional projection,        | Complete |
|       | matrix construction (4 metrics), null-model asymmetry tests   |          |
| 4     | Clustering, network community detection, sensitivity analyses | Complete |
| 5     | Bootstrap stability, VIF, figures, supplementary tables       | Complete |
| 6     | Invasion follow-up, manuscript                                | Complete |

Per-stage reports under `reports/stageN/`.

## Methods summary

- **SDMs**: ESM (ensemble of small models) via `biomod2` 4.3.4.5, dual-core
  GBM + GLM, calibrated per Stage 2 pilot (predictor set: `full_geofresh`
  primary, `climate_local` for the climate-vs-network-aware comparison).
- **Cohort**: 157 freshwater crayfish species with ≥ 30 native records after
  deduplication by hydrographic segment. *Procambarus virginalis* excluded
  for parthenogenetic CV-degeneracy (documented in Stage 2 report).
- **Predictors**: 16 GeoFRESH variables (local + upstream BIOCLIM and
  topography). Per-species `maxSSS` thresholds used for `fraction_above`.
- **Compatibility metrics**: mean predicted suitability, fraction above
  per-species maxSSS, directional Schoener's *D*, directional Warren's *I*.
- **Asymmetry test**: 1000-permutation null preserving each species'
  marginal cell distribution.
- **Quality policy**: species retained if AUC ≥ 0.7 and TSS ≥ 0.4.
  Per-species Boyce flagged unreliable if fewer than 3 valid CV folds
  (does not trigger exclusion).

Full methods in the manuscript and `reports/stage3/stage3_plan.md`.

## Reproducibility

- All paths and parameters live in `config/` and `configs/` — no hard-coded
  values in `src/`.
- Intermediate data: CSV for tables, HDF5 for large compatibility arrays.
- The submitted manuscript corresponds to commit
  [`5288e52`](https://github.com/KristianMiok/directional-compatibility-crayfish/commit/5288e52).
- A Zenodo DOI will be minted at manuscript acceptance for permanent version
  reference.

## Requirements

- Python 3.12+ (managed via `uv`; lockfile at `uv.lock`)
- R 4.3+ with `biomod2` 4.3.4.5, `ecospat`, `gbm`
- See `pyproject.toml` and `R/00_install_packages.R` for full dependencies.

## Quickstart

```bash
# Environment
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"

# Stage 1 (cohort selection)
make stage1

# Stage 2 pilot (after data drop into data/raw/)
make stage2

# Stage 3 production (long-running)
python scripts/20_run_stage3_production_sdms.py
```

The full pipeline takes ~15–25 hours of wall time, dominated by the Stage 3d
null model. See `Makefile` and per-stage scripts for entry points.

## Citation

If you use this code or data, please cite the manuscript (citation will be
added on acceptance) and link to this repository.


## License

Code is released under the MIT License (see `LICENSE`). Data are subject to
the licences of their original sources (WoC, GeoFRESH, GeoTraits — see
`docs/data_sources.md`).
