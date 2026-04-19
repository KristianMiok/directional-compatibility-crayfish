# Stage 1 Report — Manuscript #5  
**Directional ecological compatibility framework for freshwater crayfish**  
**Prepared by Kristian Miok**

## 1. Objective of Stage 1

Stage 1 had the role of auditing the available occurrence dataset, defining the candidate modellable cohort, and identifying any issues that need to be resolved before SDM calibration. The work followed the Stage 1 structure from Lucian’s plan: species inventory, threshold sensitivity, environmental coverage, basin overlap, and invasion contamination audit.

The analyses were run on the real master dataset `combined_data_true_master.csv`, which already contains crayfish occurrence records linked to the GeoFRESH river-network environmental layers.

---

## 2. Dataset used

The master dataset contained **115,191 records** and **456 species** after loading into the Stage 1 pipeline.

The available information included:
- species identity (`Crayfish_scientific_name`)
- record status (`Native`, `Alien`, `Introduced`, `Type locality`)
- record year
- original coordinates
- snapping distance to the river network
- basin and subcatchment identifiers
- large sets of local and upstream GeoFRESH-linked environmental variables

For the present Stage 1 workflow, the native filter was defined as:

- **native** = `Native` + `Type locality`
- **non-native** = `Alien` + `Introduced`

This choice was also used for the native-only environmental coverage calculations.

---

## 3. Species inventory

The inventory was built for all species in the dataset. For each species, the following quantities were calculated:

- total records
- high-accuracy records
- records snapped within 200 m
- deduplicated segment-level records
- number of basins
- convex hull extent
- temporal span
- native-only flag

The inventory showed a highly uneven distribution of record density across species. A small number of widespread taxa dominated the dataset, especially well-known invasive or broadly distributed taxa such as *Procambarus clarkii*, *Pacifastacus leniusculus*, *Faxonius limosus*, *Austropotamobius pallipes*, and *Astacus astacus*.

This confirms that threshold selection is a meaningful decision point rather than a formality.

---

## 4. Threshold sensitivity

Threshold sensitivity was assessed using `records_deduplicated_segment` as the benchmark count per species.

The number of retained species at each threshold was:

- **≥80 records**: 82 species
- **≥200 records**: 44 species
- **≥500 records**: 24 species

These values are based on the current filtering pipeline:
- high-accuracy records only
- snapping distance ≤200 m
- one record per species per subcatchment / segment

### Interpretation

The relaxed threshold (≥80) maximizes taxonomic breadth, but likely includes more species with narrower environmental representation and greater model instability.

The strict threshold (≥500) gives a very conservative cohort, but reduces the study to only 24 species and would remove many taxa that appear biologically informative and otherwise usable.

The intermediate threshold (≥200) appears to provide the best balance between model robustness and cohort diversity. It yields a cohort of 44 species, which is still large enough for a genuinely macroecological comparison while remaining restrictive enough to avoid the weakest data-poor species.

### Kristian’s recommendation for P1

At this stage, my recommendation is to use **≥200 deduplicated segment-level records** as the default threshold for Manuscript #5.

This threshold keeps the cohort broad enough for interspecific compatibility analysis, while avoiding the large drop in data quality expected under a more permissive threshold.

---

## 5. Environmental coverage

Environmental coverage was calculated for the species passing the ≥200 threshold, using **native-only records**. Eight representative GeoFRESH-linked variables were used:

1. annual mean temperature  
2. minimum temperature of coldest month  
3. maximum temperature of warmest month  
4. temperature seasonality  
5. annual precipitation  
6. precipitation seasonality  
7. stream gradient  
8. mean elevation  

This produced **352 rows** in `env_coverage.csv` (44 species × 8 variables).

### General pattern

The environmental coverage results appear biologically sensible.

Examples:
- mountain-stream taxa such as *Austropotamobius* species and *Astacus astacus* occupied the highest stream-gradient ranges
- lowland / swamp taxa such as *Cherax quadricarinatus* and several *Procambarus* species occupied flatter systems

This is an important sanity check because it indicates that the environmental variables are behaving in ecologically meaningful ways.

### Species with narrow coverage

Several species in the ≥200 cohort showed low median completeness across the selected variables, indicating relatively narrow environmental representation. These species should be treated carefully during Stage 2 because even if they pass the count threshold, they may still provide a narrower calibration envelope.

### Species with broad coverage

A few species showed broad environmental coverage and appear to span a large fraction of the crayfish-occupied environmental envelope. These taxa may become important anchors in later directional compatibility analyses.

---

## 6. Basin overlap preliminary matrix

The basin overlap matrix and Jaccard similarity matrix were generated successfully for all **456 species**.

This step was intended as a geographic sanity check rather than a niche analysis. It confirms that the dataset contains both:
- highly geographically isolated taxa
- very widespread taxa spanning many basins

At the cohort level, the number of basins per species ranged widely, from very localized taxa to extremely widespread species such as *Procambarus clarkii*.

This result will be useful later when interpreting directional compatibility asymmetries, because broad basin occupancy may reflect both niche breadth and spread history.

---

## 7. Invasion contamination audit

The contamination audit flagged all records outside the native set.

Results:
- **55,029 flagged records**
- **30 flagged species**

This confirms that invasion-related contamination is concentrated in a minority of species rather than spread uniformly across the dataset.

Within the ≥200 cohort, **18 of the 44 candidate species** have mixed native / non-native records. This means that the native-only filtering rule will have a real impact during Stage 2 calibration.

This is not a problem in itself, but it is an important methodological point: these species should be trained only on native-range data, even when total record counts appear very large.

---

## 8. Taxonomic and geographic balance

The current ≥200 cohort contains **44 species** across **10 genera**. The largest genera in the cohort are:

- *Cambarus*
- *Faxonius*
- *Procambarus*

This suggests a strong representation of North American cambarids, with smaller but still important representation of European taxa such as *Astacus*, *Austropotamobius*, *Pontastacus*, and *Pacifastacus*.

At the moment, continent and family fields are not directly available in the master dataset, so the current Stage 1 outputs do not yet include a complete continent-level or family-level stratification. If needed, these can be added later via an external taxonomy / geography enrichment step.

---

## 9. Data-quality issues and surprises

The main issues encountered in Stage 1 were:

1. **Missing direct family / country / continent fields**  
   These were requested in the original Stage 1 plan, but are not directly present in the current master dataset.

2. **Mixed native / non-native records in important species**  
   Several high-profile and high-record taxa have substantial non-native components, so strict native-only training is essential for Stage 2.

3. **Variable scaling issue in the topography codebook**  
   The stream-gradient variable required rescaling according to the GeoFRESH codebook. This has now been handled correctly in the pipeline.

4. **A small number of missing-status records**  
   These should be checked, but they do not currently affect the overall Stage 1 conclusions.

5. **Very old occurrence years in a few records**  
   A few historical outliers exist and may be worth reviewing later, but they do not currently justify additional filtering.

---

## 9.1 Native-only threshold stability check

Because the manuscript framework is explicitly native-only, I ran an additional sensitivity check within the current ≥200 cohort using the same Stage 1 quality filters (high accuracy, snap distance ≤200 m, deduplication by segment), but restricting the records to `Native` + `Type locality`.

This check showed that the current threshold recommendation is highly stable under native-only filtering:

- current cohort at ≥200: **44 species**
- still ≥200 after native-only restriction: **42 species**
- drop below threshold after native-only restriction: **2 species**

The two species that fall below the threshold are:

- **Cherax quadricarinatus**: 408 → 4 native deduplicated records
- **Pacifastacus leniusculus**: 4285 → 115 native deduplicated records

Several other widespread species remain above threshold but lose a large proportion of records when non-native occurrences are removed. The strongest examples are:

- **Faxonius limosus**: 4174 → 363
- **Procambarus clarkii**: 9606 → 1915
- **Austropotamobius fulcisianus**: 744 → 205
- **Pontastacus leptodactylus**: 1496 → 750

This result supports two conclusions:

1. the **≥200 threshold remains a robust default recommendation**
2. strict **native-only training is essential** in Stage 2, because some apparently data-rich species are strongly inflated by non-native records

## 9.1 Native-only threshold stability check

Because the framework is explicitly native-only, I ran an additional sensitivity check within the current ≥200 cohort using the same Stage 1 quality filters (high accuracy, snap distance ≤200 m, deduplication by segment), but restricting the records to `Native` + `Type locality`.

This check showed that the current threshold recommendation is highly stable under native-only filtering:

- current cohort at ≥200: **44 species**
- still ≥200 after native-only restriction: **42 species**
- drop below threshold after native-only restriction: **2 species**

The two species that fall below the threshold are:

- **Cherax quadricarinatus**: 408 → 4 native deduplicated records
- **Pacifastacus leniusculus**: 4285 → 115 native deduplicated records

Several other widespread species remain above threshold but lose a large proportion of records when non-native occurrences are removed. The strongest examples are:

- **Faxonius limosus**: 4174 → 363
- **Procambarus clarkii**: 9606 → 1915
- **Austropotamobius fulcisianus**: 744 → 205
- **Pontastacus leptodactylus**: 1496 → 750

This result supports two conclusions:

1. the **≥200 threshold remains a robust default recommendation**
2. strict **native-only training is essential** in Stage 2, because some apparently data-rich species are strongly inflated by non-native records

## 10. Conclusion for Decision Point P1

Stage 1 is now complete on the real master dataset.

Main conclusions:
- the Stage 1 workflow runs successfully end-to-end on the full real dataset
- the master dataset is suitable for the planned framework
- the candidate cohort changes strongly depending on the threshold
- the **≥200 threshold** currently appears to be the most reasonable compromise
- native-only filtering is essential for a substantial subset of candidate species
- environmental coverage outputs are biologically interpretable and useful for cohort evaluation

## 11. Proposed decision at P1

I recommend that we:

1. adopt **≥200 deduplicated segment-level records** as the cohort threshold  
2. freeze the candidate species list at this threshold after review  
3. proceed to Stage 2 SDM calibration using native-only records for training

