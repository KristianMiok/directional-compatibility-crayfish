# Manuscript #5 — Stage 3 Plan

**Author:** Kristian Miok
**Status:** Draft for Lucian's review before Stage 3 execution
**Preceded by:** Stage 2 pilot report + P2 closure email (Lucian confirmed the locked pipeline)

---

## 1. Objective

Build the directional compatibility matrix that is the paper's core data product.

For each ordered pair (A → B) of species in the locked 157-species cohort, quantify the compatibility of species A's realized niche with species B's environmental space. A → B and B → A are computed separately, producing an asymmetric N×N matrix per quality metric × predictor set × core. The asymmetry — which species are "universal donors" (niches compatible with many others) versus "universal acceptors" (habitats welcoming to many other species) — is the biological novelty the paper rests on.

## 2. Deliverables

1. **Per-species SDM artifacts** for 157 species × predictor sets × cores as per P2 lock. Stored with per-species Boyce validity flags and quality metadata.
2. **Four directional compatibility matrices** per (predictor_set, core) combination:
   - `mean_suitability[A, B]` — mean predicted suitability of A's SDM across B's environmental cells
   - `fraction_above[A, B]` — fraction of B's cells where A's predicted suitability ≥ threshold (default 0.5, confirmable)
   - `schoener_D[A, B]` — directional Schoener's D between A's predictions on B's space and B's own predictions on B's space
   - `warren_I[A, B]` — directional Warren's I, same structure
3. **Null model output** — 1000-permutation randomization of the above matrices to test asymmetry. Stored as percentile ranks per cell.
4. **Quality report** — per-species AUC / TSS / Boyce aggregates, unreliable-Boyce flags, any species hitting AUC < 0.7 or TSS < 0.4 exclusion thresholds.
5. **Stage 3 closure report** — `reports/stage3/stage3_report.md` for Lucian's P3 review.

## 3. Locked inputs from P2 (carried forward, no new decisions)

| Item | Locked value |
|---|---|
| Cohort | 157 species (158 at ≥30 threshold minus *Procambarus virginalis*) |
| Primary predictor set | `full_geofresh` (6 local CLI + 2 local TOP + 6 upstream CLI + 2 upstream TOP) |
| Piece 2 comparison set | `climate_local` (6 local BIOCLIM) |
| Primary core | GBM, tuned as in Stage 2: `n.trees=500, interaction.depth=2, shrinkage=0.01, n.minobsinnode=3, bag.fraction=0.75` |
| Secondary core | GLM on `climate_local` and `local_full` only; excluded from ensemble on `upstream_full` and `full_geofresh` |
| Species exceptions | *C. quadricarinatus*: restricted to `climate_local` + `local_full` (no upstream sets); *S. spinifrons*: ≥20 sensitivity run only, not in main matrix |
| N floor | ≥30 confirmed |
| Boyce NaN policy | Per-species Boyce aggregate requires ≥3 valid CV runs; below → flagged unreliable. Species exclusion driven by AUC < 0.7 or TSS < 0.4. Unreliable-Boyce species stay in matrix, rely on TSS+AUC. Supplementary table documents unreliable-Boyce cohort + sensitivity check. |

## 4. Pipeline architecture

Stage 3 decomposes into four substages. Each has its own driver script and output directory under `data/processed/stage3/`. Each emits a resume-on-restart progress log (JSONL) as in Stage 2.

### 4a. SDM training (Stage 3a)

Scale-up of the Stage 2 pilot batch runner to the full cohort.

- For each of 157 species, generate biomod-ready inputs for `full_geofresh` and `climate_local` (with the *Cherax* exception substituting `local_full` for `full_geofresh`).
- Run ESM-GBM (always) and ESM-GLM where allowed by the P2 rule.
- Boyce sentinel filter applied at source (R script already patched in Stage 2).
- Save per-species SDM artifact:
  - `data/processed/stage3/sdm_artifacts/<slug>/<predictor_set>__<core>_biomod2.rds` — the serialized biomod2 object
  - `data/processed/stage3/sdm_artifacts/<slug>/<predictor_set>__<core>_metadata.json` — n_records, AUC/TSS/Boyce per CV run + aggregates, Boyce validity flag, fit_status

**Estimated runtime**: ~5–6 hours wall time on a laptop. Dominated by high-N Cambaridae × GBM on `full_geofresh`.

Script: `scripts/20_run_stage3_production_sdms.py` (extends `scripts/12_run_stage2_pilot_batch.py`).

### 4b. Directional projection (Stage 3b) — novel compute

For each trained SDM (species A, predictor_set, core) and each other species B: predict suitability on B's environmental cells.

**Definition of B's environmental space** (methodological choice, flagged to Lucian):
B's environmental cells = the environmental feature vectors at each of B's native-range occurrence records, after the same filter chain used at training (Accuracy = High, distance_m ≤ 200, Status ∈ {Native, Type locality}, deduplicated by subc_id). This is the simplest, most consistent choice: it asks "does A's SDM predict high suitability at the specific environments where B is actually found?" — which matches how Schoener's D and Warren's I are computed in the niche-overlap literature. Alternatives (full gridded range clipped to range polygon; basin-level aggregates) were considered and are discussed in §8.

Output: one parquet file per (A, predictor_set, core):
`data/processed/stage3/projections/<A_slug>__<predictor_set>__<core>_projections.parquet`

Columns: `target_species`, `record_idx`, `pred_suitability`. Rows cover every cell of every B in the cohort (minus A itself and B's that failed QC).

**Estimated**: 157 × 156 = 24,492 projection calls per (predictor_set, core). Each is cheap (single-model predict on 30–10,000 rows). Total runtime estimated 1–3 hours.

Script: `scripts/21_project_directional.py`.

### 4c. Matrix aggregation (Stage 3c)

Collapse the projection parquet files into the four N×N matrices per (predictor_set, core).

- `mean_suitability[A, B]` — mean of A's predictions on B's cells
- `fraction_above[A, B]` — mean of (A's predictions on B's cells ≥ τ), τ = 0.5 by default
- `schoener_D[A, B]` — 1 - 0.5 × Σ|p_A(B_cells) - p_B(B_cells)|, using normalized prediction distributions
- `warren_I[A, B]` — Hellinger-based similarity, same structure

Output: HDF5 file(s) per (predictor_set, core):
`data/processed/stage3/matrices/<predictor_set>__<core>_matrices.h5`

Internal layout:
```
/mean_suitability          (157, 157) float32
/fraction_above            (157, 157) float32
/schoener_D                (157, 157) float32
/warren_I                  (157, 157) float32
/species_order             (157,) string — row/column species names in order
/attrs:
    predictor_set = "full_geofresh"
    core = "GBM"
    threshold = 0.5
    qc_pass_mask = (157,) bool — which species passed AUC/TSS exclusion
    boyce_valid_mask = (157,) bool — which species have reliable Boyce
    n_records_per_species = (157,) int32
    build_timestamp = "..."
```

HDF5 was chosen because (a) the N×N matrix is the natural 2D structure, (b) compression matters at this size, (c) single-file distribution is methods-paper-friendly, (d) Python (`h5py`), R (`rhdf5`), and MATLAB/Julia all read it.

Script: `scripts/22_build_matrices.py`.

### 4d. Null model (Stage 3d)

Randomization-based test for matrix asymmetry. For each of 1000 permutations: shuffle each species' environmental cells among target species (permute the "target identity" label while preserving each species' own cell distribution), rebuild the four matrices, record their asymmetry statistic.

Null distribution for each matrix cell → percentile rank for observed value → p-value for asymmetry.

Output: percentile-rank HDF5 file per (predictor_set, core):
`data/processed/stage3/matrices/<predictor_set>__<core>_null_percentiles.h5`

**Estimated**: 1000 perms × 4 matrices × 2 (predictor_set, core) combinations × ~30 s per perm = ~30–40 hours if naively serial. Would parallelize across permutations (embarrassingly parallel) on ~8 cores for ~4–5 hours. **This is the compute bottleneck of Stage 3.**

Script: `scripts/23_null_model.py`.

## 5. Data schemas and storage

### Directory layout (new)

```
data/processed/stage3/
├── sdm_artifacts/
│   └── <species_slug>/
│       ├── full_geofresh__GBM_biomod2.rds
│       ├── full_geofresh__GBM_metadata.json
│       ├── climate_local__GBM_biomod2.rds
│       ├── climate_local__GBM_metadata.json
│       ├── climate_local__GLM_biomod2.rds (where applicable)
│       └── climate_local__GLM_metadata.json
├── projections/
│   └── <A_slug>__<predictor_set>__<core>_projections.parquet
├── matrices/
│   ├── full_geofresh__GBM_matrices.h5
│   ├── full_geofresh__GBM_null_percentiles.h5
│   ├── climate_local__GBM_matrices.h5
│   ├── climate_local__GBM_null_percentiles.h5
│   ├── ensemble__climate_local_matrices.h5         # GBM + GLM averaged, where allowed
│   └── ensemble__local_full_matrices.h5
├── qc/
│   ├── per_species_quality.csv
│   ├── excluded_species.csv
│   └── unreliable_boyce_species.csv
└── batch_progress/
    ├── stage3a_sdm_progress.jsonl
    ├── stage3b_projection_progress.jsonl
    ├── stage3c_matrix_progress.jsonl
    └── stage3d_null_progress.jsonl
```

### What gets committed to git

- Scripts and R code: yes
- `data/processed/stage3/` as a whole: gitignored (existing rule)
- Exceptions (travel with the repo):
  - `data/processed/stage3/qc/per_species_quality.csv`
  - `data/processed/stage3/qc/excluded_species.csv`
  - `data/processed/stage3/qc/unreliable_boyce_species.csv`
  - The `.h5` matrix files are too large (and binary-diff-unfriendly) for git; they go in a separate release artifact or cloud storage for distribution.

## 6. Compute plan

| Substage | Rough runtime | Parallelizable? |
|---|---|---|
| 3a SDM training | 5–6 h | Yes, by species (batch runner already parallel-safe in principle) |
| 3b Directional projection | 1–3 h | Yes, by A species |
| 3c Matrix aggregation | 15–30 min | Trivially, by (predictor_set, core) |
| 3d Null model | 4–5 h if 8-core, 30–40 h if serial | Yes, by permutation |
| **Total** | **~10–15 h wall time on the laptop if parallelized** | |

If 3d proves too slow, consider renting 16-core cloud compute for the null model run (~€5–15 for a few hours). Not needed unless laptop timing is prohibitive.

## 7. Quality control and exclusions

Per Lucian's P2 email:

- **Species-level exclusion** at Stage 3 only occurs if AUC < 0.7 **or** TSS < 0.4 on the primary predictor set (`full_geofresh` or `local_full` depending on species exception). A species hitting either threshold is flagged excluded and does not enter the matrix as either source (A) or target (B).
- **Boyce unreliability** (aggregate from < 3 valid CV runs) does *not* trigger exclusion. Affected species stay in the matrix with an `unreliable_boyce_flag = True` in the per-species quality report.
- **Pre-known exceptions** from P2: *P. virginalis* excluded; *Cherax quadricarinatus* uses `local_full` not `full_geofresh`; *Samastacus spinifrons* in sensitivity set only.

**Quality gate before advancing to matrix construction** (end of Stage 3a): produce `data/processed/stage3/qc/per_species_quality.csv` and review with Lucian before building the matrices. This is a soft checkpoint, not a hard gate — but if 20+ species fail QC it would restructure the analysis and we want to know early.

## 8. Open questions — decisions needed from Lucian before execution

These are the only substantive choices not already locked by P2. I've put a recommended default on each, but they are Lucian's call.

1. **B's environmental space definition (§4b above).** Recommended default: B's native-range occurrence env vectors. Alternatives: full gridded env clipped to a range polygon (requires range maps we don't have); basin-level aggregation (one vector per unique basin_id). The choice changes the ecological interpretation of the compatibility score and should be Lucian's decision.
2. **Suitability threshold τ for `fraction_above` matrix.** Recommended default: 0.5 (conventional biomod2 threshold). Alternatives: per-species TSS-optimized cutoff (requires a second pass on the SDMs); fixed 0.4 or 0.6 (simpler but arbitrary).
3. **Ensemble weighting for GLM+GBM averaging** on `climate_local` and `local_full`. Recommended default: equal-weighted mean of the two cores' predictions. Alternatives: AUC-weighted; only use GLM where Boyce is reliable for that species. This is a minor Stage 3c choice.
4. **Null model permutation scheme.** Recommended default: permute target-species labels on each species' cell vectors (preserves each species' marginal cell distribution, randomizes pair identity). Alternative: full random shuffle of A↔B labels per iteration. The default is more conservative (preserves more structure); the alternative is closer to standard randomization.
5. **Do we need to re-run the Stage 2 pilot comparison with the corrected Boyce filter and issue a revised Stage 2 report, or is P2's closure sufficient given the corrections were computed and shared via email?** Not blocking Stage 3; just confirming whether a formal revised report is expected before Stage 4.

## 9. P3 closure checklist

Before proposing P3 closure to Lucian, the deliverables in §2 must be produced, plus:

- [ ] Per-species quality report reviewed, any unexpected exclusions discussed
- [ ] Unreliable-Boyce supplementary table drafted
- [ ] Null-model asymmetry test run and p-values distributed sensibly (not all significant, not all non-significant — either would suggest the null is malformed)
- [ ] Three figures ready for Stage 4 inspection: raw N×N matrix heatmap (primary), asymmetry histogram (A→B vs B→A), climate-only vs full-GeoFRESH comparison for one or two illustrative species pairs
- [ ] Runtime log and memory profile (useful for the methods section later)

## 10. Risk register

| Risk | Likelihood | Mitigation |
|---|---|---|
| Low-N species produce unstable SDMs at full-cohort scale that weren't caught in the 26-species pilot | Moderate | QC gate after Stage 3a; review before building matrices |
| Boyce NaN rate higher than pilot (pilot was 30–44% on GLM, could be worse on the full low-N tail) | Moderate | Policy from P2 already handles this; flag it in supplementary table |
| Null model runtime prohibitive on laptop | Moderate | Parallelize across permutations; move to cloud compute if needed |
| Matrix storage size larger than expected | Low | HDF5 compression expected to bring each matrix file to ~10–30 MB |
| *Cherax quadricarinatus* `local_full` still produces unstable fits (pilot had 4/8 GBM failures there) | Moderate | Already flagged; may need second-level fallback to `climate_local` only |
| Projection fails on small subset of (A, B) pairs due to NaN in B's env vectors | Moderate | Document and report; handle by either dropping those cells or imputing |
| One or more high-N species (e.g., *P. clarkii*, N ≈ 10,000) takes much longer than projected on GBM | Moderate | Already seen in pilot; add per-species timeout + flag |

## 11. What Stage 3 does NOT include

To keep scope tight:

- **No traits integration** — that's Stage 4 (per Lucian's plan doc).
- **No clustering or community detection** on the matrix — Stage 4.
- **No network centrality analysis** — Stage 4.
- **No Mantel test** between climate-only and full-GeoFRESH matrices — Stage 4.
- **No manuscript writing** — Stage 6.

Stage 3 ends with the matrices, the null, and the QC report. Stage 4 begins the analytical interpretation.
'''

Path("/home/claude/stage3_plan.md").write_text(plan)
print(f"Plan written: {len(plan):,} chars, {plan.count(chr(10))} lines")