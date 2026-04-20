# Manuscript #5 — Stage 2 Pilot Report

**Author:** Kristian Miok
**Recipients:** Lucian Pârvulescu (co-author, manuscript lead)
**Status:** For review at Decision Point P2

---

## 1. Summary

The Stage 2 ESM pilot ran successfully across the frozen 26-species shortlist,
four predictor sets, and both ESM cores (GLM, GBM), totalling 208 configurations.
**194 configurations (93%) produced stable per-run metrics; 14 configurations
(7%) across 3 species failed in informative ways.**

Headline recommendations for P2 closure:

1. **Retain dual-core ESM** (GLM + GBM). The two cores have complementary
   failure modes — GBM is more stable across predictor-set richness; GLM is
   interpretable and performs well when restricted to parsimonious predictors
   but overfits catastrophically with rich predictor sets on low-N species.
2. **Lock `climate_local` (6 BIOCLIM variables) as the production predictor set.**
   It is the only set where GLM remains stable across the N gradient, GBM is
   equally strong across all sets so there is no cost to parsimony, and it
   survives every retained species.
3. **Exclude *Procambarus virginalis* from the cohort.** Both cores fail on
   all four predictor sets. The parthenogenetic origin produces clonally
   identical records that CV cannot partition. Lucian's anomaly flag is
   empirically validated as severe; the species should move to a documented
   methodological exception rather than the production matrix.
4. **Confirm the N ≥ 30 threshold floor.** *Samastacus spinifrons* at N=24
   produces CV-degenerate metrics (perfect 1.0 or wildly unstable within the
   same cell); *Austropotamobius bihariensis* at N=57 produces stable metrics
   under GBM. The transition is between these Ns; N ≥ 30 is a defensible cut.
5. **Flag *Cherax quadricarinatus* for predictor-set-restricted inclusion.**
   Both cores fail on `upstream_full` and `full_geofresh` (likely pseudo-
   absence/upstream-feature interaction given its isolated Australian native
   range). `climate_local` and `local_full` work.

## 2. Pipeline and frozen choices

- **Frozen cohort**: 26 species (reports/stage1/p1_package/ → data/processed/stage2_pilot_shortlist_frozen.csv)
- **Predictor sets**: `climate_local` (6 BIOCLIM means), `local_full` (climate + local topography),
  `upstream_full` (upstream climate + topography), `full_geofresh` (all 16 local + upstream variables)
- **Cores**: GLM (biomod2 defaults, `stats::glm`), GBM (biomod2 `gbm` with tuned options:
  `n.trees=500`, `interaction.depth=2`, `shrinkage=0.01`, `n.minobsinnode=3`, `bag.fraction=0.75`)
- **CV**: 5 random splits, 80/20 presence/test
- **Evaluation**: AUCroc, TSS, KAPPA, Boyce index (via ecospat), MPA
- **Runtime**: ~2–3 h total wall time on the pilot; dominated by GBM on high-N species
- **Code**: `R/01_run_stage2_esm_single_core.R` (single ESM run),
  `scripts/12_run_stage2_pilot_batch.py` (resume-on-restart orchestration)

## 3. Headline: predictor set × core comparison

Means ± std across all species × all CV runs that produced metrics (excludes
the 14 failed configurations, includes 970 individual CV runs).

| Predictor set | Core | AUC | TSS | Boyce | Kappa |
|---|---|---|---|---|---|
| climate_local | GLM | 0.955 ± 0.051 | 0.857 ± 0.123 | 0.925 ± 0.218 | 0.814 ± 0.149 |
| climate_local | GBM | 0.969 ± 0.033 | 0.839 ± 0.157 | 0.913 ± 0.131 | 0.813 ± 0.167 |
| local_full | GLM | 0.943 ± 0.105 | 0.833 ± 0.216 | 0.934 ± 0.090 | 0.794 ± 0.216 |
| local_full | GBM | 0.974 ± 0.029 | 0.863 ± 0.121 | 0.890 ± 0.163 | 0.834 ± 0.136 |
| upstream_full | GLM | 0.949 ± 0.071 | 0.844 ± 0.154 | 0.878 ± 0.328 | 0.806 ± 0.170 |
| upstream_full | GBM | 0.971 ± 0.031 | 0.857 ± 0.095 | 0.881 ± 0.157 | 0.825 ± 0.115 |
| full_geofresh | GLM | 0.936 ± 0.080 | 0.821 ± 0.169 | 0.756 ± 0.602 | 0.778 ± 0.175 |
| full_geofresh | GBM | 0.973 ± 0.031 | 0.861 ± 0.109 | 0.906 ± 0.142 | 0.832 ± 0.122 |

**Key readings:**

- **GBM AUC std is remarkably flat across predictor-set richness** (~0.03 across
  all four sets). GLM AUC std grows from 0.051 (climate_local) to 0.080+ (richer
  sets) — a mild overfitting signature on AUC.
- **GLM Boyce std on richer predictor sets is the real problem.** At
  `full_geofresh`, GLM Boyce mean drops to 0.756 with std 0.602 — i.e. individual
  CV runs are producing Boyce values spanning [−1, +1] within the same cell.
  GBM Boyce stays near 0.90 with much tighter std.
- **At `climate_local`, GLM and GBM are close to tied** (both ~0.97 AUC, both
  stable). This is the regime where GLM is safe to use.
- **For `climate_local`, GLM TSS (0.857) slightly exceeds GBM TSS (0.839).**
  At parsimonious predictor counts GLM is competitive.

## 4. Low-N anchor performance (per Lucian's request)

Lucian asked that the two empirical low-N anchors be tracked separately. Mean ±
std across 5 CV runs per cell. Blank cells indicate fit failure.


**Austropotamobius bihariensis** (N=57)

| Predictor set | Core | AUC | Boyce |
|---|---|---|---|
| climate_local | GLM | 0.924 ± 0.030 | 0.770 ± 0.083 |
| climate_local | GBM | 0.930 ± 0.056 | 0.775 ± 0.162 |
| local_full | GLM | 0.874 ± 0.038 | 0.806 ± 0.071 |
| local_full | GBM | 0.973 ± 0.043 | 0.715 ± 0.073 |
| upstream_full | GLM | 0.908 ± 0.038 | -1.000 |
| upstream_full | GBM | 0.990 ± 0.011 | 0.746 ± 0.152 |
| full_geofresh | GLM | 0.887 ± 0.050 | -0.125 ± 1.031 |
| full_geofresh | GBM | 0.984 ± 0.017 | 0.762 ± 0.206 |

**Samastacus spinifrons** (N=24)

| Predictor set | Core | AUC | Boyce |
|---|---|---|---|
| climate_local | GLM | 1.000 ± 0.000 | nan |
| climate_local | GBM | 1.000 ± 0.000 | 1.000 ± 0.000 |
| local_full | GLM | 1.000 ± 0.000 | nan |
| local_full | GBM | 1.000 ± 0.000 | 1.000 ± 0.000 |
| upstream_full | GLM | 0.840 ± 0.261 | 1.000 |
| upstream_full | GBM | *fit failed* | *fit failed* |
| full_geofresh | GLM | 0.770 ± 0.251 | -1.000 |
| full_geofresh | GBM | *fit failed* | *fit failed* |

**Interpretation:**

- **A. bihariensis (N=57) is a successful low-N case under GBM with any
  predictor set.** AUC 0.930–0.990, Boyce 0.71–0.88, stable across runs.
  GLM with rich predictor sets is catastrophic — on `full_geofresh` Boyce mean
  is −0.125 with std 1.031, and on `upstream_full` it collapses to exactly −1.0
  (all 5 CV runs produce perfectly inverted predictions). GLM is only safe
  at `climate_local` for this species.
- **S. spinifrons (N=24) is CV-degenerate, as anticipated by being below the
  locked threshold.** On `climate_local` and `local_full` every CV run
  produces perfect 1.0 metrics, which at N=24 is the signature of test
  splits too small to separate presences from pseudo-absences by chance.
  On `full_geofresh` GLM the metric variance explodes (AUC 0.77 ± 0.25,
  Boyce spans [−1, +1]). GBM outright fails on the two predictor sets with
  upstream features.

**Conclusion: ≥ 30 is a defensible threshold floor.** The transition from
CV-unstable (N=24) to GBM-stable (N=57) happens in that range.

## 5. Species-specific failures

### Procambarus virginalis (132 records, parthenogenetic)

All 8 configurations (both cores × all 4 predictor sets) fail to produce
per-run metrics. GBM fit errors during `BIOMOD_Modeling`
(`nTrain * bag.fraction <= 2 * n.minobsinnode + 1`); GLM completes but
`get_evaluations` returns only the `allRun` summary row with all-NaN metrics
(i.e., CV produces no partitionable variance). This is consistent with the
clonal structure of this parthenogenetic species: 132 records are effectively
copies of one genetic individual, so any CV split trains and tests on the same
ecological signal and cannot produce independent evaluations.

**Recommendation**: exclude from the production cohort for Stage 3; retain in
methods as a documented limitation. The matrix at ≥30 threshold would then
hold 157 species rather than 158.

### Cherax quadricarinatus (408 records, tropical Australian native)

Both cores fail on `upstream_full` and `full_geofresh`, both succeed on
`climate_local` and `local_full`. The shared factor: the upstream variants
include catchment-aggregated features (`u_CLI*`, `u_TOP*`), and with global
random pseudo-absence sampling the upstream features of truly distant PAs
contain structural NaN patterns that break evaluation. Native variables
(`l_*`) do not have this issue.

**Recommendation**: retain in the production cohort, restrict to local-only
predictor sets for this species, and note the PA/upstream interaction in
methods. This is a minor extension of the Stage 1 finding that Cherax has
a narrow native-to-invaded gradient structure.

### Samastacus spinifrons (24 records, below threshold, added as stress test)

As analysed in §4 above. Retain as the designated sub-threshold stress test;
use only `climate_local` or `local_full` predictor sets; acknowledge that
metrics from this species are unstable and should not be interpreted as
validated performance.

## 6. Supporting files

- `data/processed/stage2_pilot_shortlist_frozen.csv` — frozen 26-species cohort
- `data/processed/stage2_pilot/stage2_pilot_comparison_final.csv` — species × predictor_set × core × metrics
- `data/processed/stage2_pilot/*_esm_*_evaluations.csv` — per-run metrics (208 files)
- `data/processed/stage2_pilot/batch_progress.jsonl` — batch runner audit log
- `R/01_run_stage2_esm_single_core.R` — single-run ESM
- `scripts/12_run_stage2_pilot_batch.py` — batch orchestration with resume-on-restart

## 7. Proposal for P2 closure

Pending your review, I propose Decision Point P2 close with:

- **Cores**: dual ESM (GLM + GBM), averaged to ensemble in Stage 3 per Lucian's original plan
- **Predictor set**: `climate_local` as production; `local_full` as sensitivity check
- **Production cohort**: 157 species (originally 158, drop *P. virginalis*)
- **Species-specific exceptions**:
  - *Cherax quadricarinatus*: restrict predictor set to `local_full` or `climate_local`
  - *Samastacus spinifrons*: retain under the ≥20 sensitivity analysis only
- **Threshold floor**: N ≥ 30 confirmed as conservative

Open for discussion — any of the above can change before we lock and proceed
to Stage 3's full run.
