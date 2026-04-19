# Etapa 1 — Stage 1 Report (Data Audit and Species Selection)

**Author:** Kristian Miok
**Date:** _to fill_
**Status:** _draft / for review / closed at P1_

---

## 1. Summary

_2–3 sentences: how much data we have, how many species are modellable
under each threshold, the recommended threshold, and the main data-quality
issues encountered._

## 2. Inventory results (Task 1.1)

- Total species in WoC export: _N_
- After high-accuracy filter: _N_
- After ≤200 m snap distance: _N_
- After segment deduplication: _N_

Reference: `data/processed/species_inventory.csv`

## 3. Threshold sensitivity (Task 1.2)

| Threshold | Species retained | Continents covered | Families covered |
|-----------|------------------|--------------------|------------------|
| ≥ 80      |                  |                    |                  |
| ≥ 200     |                  |                    |                  |
| ≥ 500     |                  |                    |                  |

Plot: `reports/stage1/threshold_sensitivity_plot.png`

**Recommendation (Kristian):** _which threshold and why._
Considerations: taxonomic coverage, geographic balance, downstream SDM
sample-size requirements (Stage 2 pilot needs ≥20–30 species with reliable
fits at the chosen threshold).

## 4. Environmental coverage (Task 1.3)

Summary statistics across the candidate cohort:

- Median completeness ratio per variable (across species): _fill_
- Species in the bottom decile of completeness for ≥3 variables: _list_

Reference: `data/processed/env_coverage.csv`,
per-species plots in `reports/stage1/completeness/`.

**Flag for P1:** species with consistently low environmental coverage —
these may model technically but extrapolate badly when used as the *target*
in the directional matrix. Consider a "modelled but flagged" tier.

## 5. Basin overlap (Task 1.4)

- Number of geographically isolated species (no shared basins with any
  other in the cohort): _N_
- Median Jaccard similarity in cohort: _value_
- Mention any surprising overlaps or surprising independences.

Reference: `data/processed/basin_overlap_matrix.csv` and
`basin_jaccard_matrix.csv`.

## 6. Invasion contamination audit (Task 1.5)

- Records flagged: _N_ across _N_ species
- Breakdown by reason: _outside_native_continent / outside_native_basin_post_year
  / source_marked_invasive / unknown_native_range_

Reference: `data/processed/invasion_contamination_flagged.csv`.

**Action items:**

- [ ] Lucian, Mihaela, Dave — manual review of flagged records
- [ ] Decide on default treatment for `unknown_native_range`
  (provisionally: keep but mark in metadata)

## 7. Data-quality observations and surprises

_Any schema mismatches in the WoC export, GeoFRESH coverage gaps, low
temporal density species, suspicious clusters, etc._

## 8. Recommendation for P1

- **Final threshold:** _value_
- **Locked cohort:** _N species_ (full list in Appendix A)
- **Repository:** _new repo `directional-compatibility-crayfish` confirmed_
  / _continue in HybridSDM_

---

### Appendix A — Locked species list (after P1)

_To be appended after the P1 meeting with Lucian._
