# Manuscript #5 — Work Plan (reference copy)

> This is a verbatim copy of Lucian's plan for in-repo reference. The
> canonical version lives in the original `Manuscript5_Plan_Kristian.docx`.
> If the plan is updated, replace this file too.

## 1. Overview

Macroecological framework paper targeting **Ecography**. Introduces a
**directional ecological compatibility framework** for freshwater crayfish
and applies it systematically across all species with sufficient native-range
data.

**Core concept.** For each species A, train an SDM on native records only.
Project A onto the native environmental space of every other species B,
producing an asymmetric compatibility matrix. Analyse the matrix as a
macroecological object via asymmetry, clustering, and network metrics.
Supplement with (i) climate-only vs. full GeoFRESH comparison and
(ii) integration with GeoTraits ecological traits.

**Scope constraint.** Projections are restricted to the
**crayfish-occupied freshwater network envelope**. No projection onto the
full GeoFRESH global network in this paper.

**Validation.** Prospective and native-only by design. Internal validation
only: null models for asymmetry, bootstrap stability, algorithm comparison,
provenance-stratified robustness. Invasion-outcome validation belongs to the
companion paper (Miok, Pârvulescu et al., currently at *Journal of
Biogeography*).

## 2. Working principles

- Each stage has explicit tasks, expected outputs, and a decision point.
- Kristian delivers, Lucian reviews, decision point closed jointly.
- No advance to stage N+1 without closing decision point at stage N.
- All outputs on GitHub.
- Intermediate data: CSV (tables), GeoJSON (geometries), HDF5 (large arrays).
- Pipeline frozen after Stage 2.

## 3. Stage 1 — Data Audit and Species Selection (~2 weeks)

| Task | Deliverable                                | Notes                                                |
|------|--------------------------------------------|------------------------------------------------------|
| 1.1  | `species_inventory.csv`                    | One row per species, all required count columns      |
| 1.2  | `threshold_sensitivity.csv` + plot         | Thresholds: 80 (relaxed), 200 (default), 500 (strict)|
| 1.3  | `env_coverage.csv` + per-species plots     | 6 key variables vs. global crayfish-occupied envelope|
| 1.4  | `basin_overlap_matrix.csv`                 | Sanity check on geographic independence              |
| 1.5  | `invasion_contamination_flagged.csv`       | For manual review by Lucian/Mihaela/Dave             |

**Stage deliverable:** `etapa1_raport.md` (2–3 pages) summarizing findings
and Kristian's recommended threshold.

**Decision Point P1:** Final threshold and **locked species list**.

## 4. Downstream stages (overview only — detail issued after each P-point)

| Stage | Content                                             | Duration   | Cumulative   |
|-------|-----------------------------------------------------|------------|--------------|
| 1     | Data audit & species selection                      | 2 weeks    | 2 weeks      |
| 2     | SDM pipeline calibration (pilot, 20–30 species)     | 3 weeks    | 5 weeks      |
| 3     | Full SDM run + directional matrix construction      | 4–6 weeks  | 9–11 weeks   |
| 4     | Analytical layer (asymmetry, network, traits)       | 3–4 weeks  | 12–15 weeks  |
| 5     | Figures & supplementary                             | 2 weeks    | 14–17 weeks  |
| 6     | Manuscript writing (Lucian lead)                    | 4–6 weeks  | 18–23 weeks  |

Total: ~4.5–6 months. See the source `.docx` for Stage 2–6 previews.
