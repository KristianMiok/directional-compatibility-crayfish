#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
STAGE4A = DATA / "stage3" / "stage4a"
OUTDIR = STAGE4A

CLUSTER_FP = STAGE4A / "table_stage4a_cluster_membership_k6.csv"
INVENTORY_FP = DATA / "species_inventory.csv"


def safe_metric_pair(cluster_labels, external_labels):
    ari = adjusted_rand_score(external_labels, cluster_labels)
    nmi = normalized_mutual_info_score(external_labels, cluster_labels)
    return float(ari), float(nmi)


def summarize_label_distribution(df: pd.DataFrame, label_col: str) -> pd.DataFrame:
    out = (
        df.groupby(["cluster", label_col], dropna=False)
        .size()
        .reset_index(name="n_species")
        .sort_values(["cluster", "n_species", label_col], ascending=[True, False, True])
        .reset_index(drop=True)
    )
    return out


def main() -> None:
    clusters = pd.read_csv(CLUSTER_FP)
    inventory = pd.read_csv(INVENTORY_FP)

    merged = clusters.merge(
        inventory,
        left_on="species",
        right_on="species_name",
        how="left",
        validate="one_to_one",
    )

    missing = merged["species_name"].isna().sum()
    print(f"Merged rows: {len(merged)}")
    print(f"Missing inventory matches: {missing}")

    if missing > 0:
        print("\nSpecies missing from inventory:")
        print(merged.loc[merged["species_name"].isna(), "species"].to_string(index=False))

    # Keep only rows with external labels present for each metric
    metrics_rows = []

    comparisons = [
        ("family", "family"),
        ("genus", "genus"),
        ("modal_continent", "modal_continent"),
    ]

    for label_name, col in comparisons:
        sub = merged.dropna(subset=[col]).copy()
        ari, nmi = safe_metric_pair(sub["cluster"], sub[col])
        metrics_rows.append({
            "comparison": f"cluster_vs_{label_name}",
            "n_species_used": int(len(sub)),
            "n_unique_clusters": int(sub["cluster"].nunique()),
            "n_unique_external_labels": int(sub[col].nunique()),
            "adjusted_rand_index": ari,
            "normalized_mutual_info": nmi,
        })

    metrics_df = pd.DataFrame(metrics_rows)

    family_breakdown = summarize_label_distribution(
        merged.dropna(subset=["family"]).copy(),
        "family",
    )
    genus_breakdown = summarize_label_distribution(
        merged.dropna(subset=["genus"]).copy(),
        "genus",
    )
    continent_breakdown = summarize_label_distribution(
        merged.dropna(subset=["modal_continent"]).copy(),
        "modal_continent",
    )

    # For readability: dominant label per cluster
    def dominant_label_table(breakdown: pd.DataFrame, label_col: str) -> pd.DataFrame:
        dom = (
            breakdown.sort_values(["cluster", "n_species", label_col], ascending=[True, False, True])
            .groupby("cluster", as_index=False)
            .first()
            .rename(columns={label_col: f"dominant_{label_col}", "n_species": f"dominant_{label_col}_count"})
        )
        return dom

    dominant_family = dominant_label_table(family_breakdown, "family")
    dominant_genus = dominant_label_table(genus_breakdown, "genus")
    dominant_continent = dominant_label_table(continent_breakdown, "modal_continent")

    cluster_summary = (
        merged.groupby("cluster", as_index=False)
        .size()
        .rename(columns={"size": "cluster_size"})
        .merge(dominant_family, on="cluster", how="left")
        .merge(dominant_genus, on="cluster", how="left")
        .merge(dominant_continent, on="cluster", how="left")
        .sort_values("cluster")
        .reset_index(drop=True)
    )

    metrics_df.to_csv(OUTDIR / "table_stage4a_cluster_concordance_metrics.csv", index=False)
    cluster_summary.to_csv(OUTDIR / "table_stage4a_cluster_summary_k6.csv", index=False)
    family_breakdown.to_csv(OUTDIR / "table_stage4a_cluster_family_breakdown.csv", index=False)
    genus_breakdown.to_csv(OUTDIR / "table_stage4a_cluster_genus_breakdown.csv", index=False)
    continent_breakdown.to_csv(OUTDIR / "table_stage4a_cluster_continent_breakdown.csv", index=False)

    summary = {
        "cluster_file": str(CLUSTER_FP),
        "inventory_file": str(INVENTORY_FP),
        "n_species_clusters": int(len(clusters)),
        "n_species_inventory_matched": int(len(merged) - missing),
        "n_species_inventory_missing": int(missing),
        "comparisons_run": [r["comparison"] for r in metrics_rows],
    }
    (OUTDIR / "stage4a_concordance_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {OUTDIR / 'table_stage4a_cluster_concordance_metrics.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_cluster_summary_k6.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_cluster_family_breakdown.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_cluster_genus_breakdown.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_cluster_continent_breakdown.csv'}")
    print(f"  {OUTDIR / 'stage4a_concordance_summary.json'}")

    print("\nConcordance metrics:")
    print(metrics_df.to_string(index=False))

    print("\nCluster summary:")
    print(cluster_summary.to_string(index=False))


if __name__ == "__main__":
    main()