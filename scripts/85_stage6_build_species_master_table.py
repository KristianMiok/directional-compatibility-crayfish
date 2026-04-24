#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
OUTDIR = DATA / "stage6_invasion_followup"
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIMARY_MATRIX_FP = DATA / "stage3" / "matrices_full_geofresh_gbm.h5"

OUT_CSV = OUTDIR / "species_master_table_base.csv"
OUT_JSON = OUTDIR / "species_master_table_base_summary.json"


def read_csv_if_exists(fp: Path) -> pd.DataFrame | None:
    if fp.exists():
        return pd.read_csv(fp)
    return None


def ensure_species_col(df: pd.DataFrame, candidates: list[str]) -> pd.DataFrame:
    for c in candidates:
        if c in df.columns:
            if c != "species":
                df = df.rename(columns={c: "species"})
            return df
    raise KeyError(f"Could not find species column among {candidates}. Found columns: {list(df.columns)}")


def load_primary_species(h5_fp: Path) -> list[str]:
    with h5py.File(h5_fp, "r") as h5:
        return [
            s.decode("utf-8") if isinstance(s, (bytes, bytes)) else str(s)
            for s in h5["species"][:]
        ]


def main() -> None:
    primary_species = load_primary_species(PRIMARY_MATRIX_FP)
    primary_set = set(primary_species)

    inventory_fp = DATA / "species_inventory.csv"
    donors_fp = DATA / "stage3" / "tables" / "table_stage3_top20_donors_acceptors.csv"
    network_fp = DATA / "stage3" / "stage4b" / "table_stage4b_species_network_metrics.csv"
    cluster_fp = DATA / "stage3" / "stage4a" / "table_stage4a_cluster_membership_k6.csv"
    community_fp = DATA / "stage3" / "stage4b" / "table_stage4b_network_communities.csv"
    reposition_fp = DATA / "stage3" / "stage4c" / "table_stage4c_top50_repositioners.csv"
    residual_fp = DATA / "stage3" / "stage4c" / "table_stage4c_top50_species_pairwise_residuals.csv"
    bootstrap_fp = DATA / "stage3" / "stage4d" / "summary" / "table_stage4d_species_top20_frequencies.csv"

    inventory = read_csv_if_exists(inventory_fp)
    if inventory is None:
        raise FileNotFoundError(f"Missing required file: {inventory_fp}")

    inventory = ensure_species_col(inventory, ["species", "species_name"])
    inventory = inventory[inventory["species"].isin(primary_set)].copy()

    keep_cols = [c for c in [
        "species",
        "family",
        "genus",
        "total_records",
        "records_high_accuracy",
        "records_snapped_le_200m",
        "records_deduplicated_segment",
        "records_after_thinning",
        "n_basins",
        "n_countries",
        "n_continents",
        "modal_continent",
        "range_bbox_area_km2",
        "native_only_flag",
        "temporal_span_years",
        "temporal_density",
        "ecological_strategy",
    ] if c in inventory.columns]
    master = inventory[keep_cols].copy()

    donors = read_csv_if_exists(donors_fp)
    if donors is not None:
        donors = ensure_species_col(donors, ["species"])
        if "role" in donors.columns:
            donor_rows = donors[donors["role"].astype(str).str.lower() == "donor"].copy()
            acceptor_rows = donors[donors["role"].astype(str).str.lower() == "acceptor"].copy()

            donor_merge_cols = ["species"]
            acceptor_merge_cols = ["species"]

            if "score" in donor_rows.columns:
                donor_rows = donor_rows.rename(columns={"score": "donor_score"})
                donor_merge_cols.append("donor_score")
            if "rank" in donor_rows.columns:
                donor_rows = donor_rows.rename(columns={"rank": "donor_rank"})
                donor_merge_cols.append("donor_rank")

            if "score" in acceptor_rows.columns:
                acceptor_rows = acceptor_rows.rename(columns={"score": "acceptor_score"})
                acceptor_merge_cols.append("acceptor_score")
            if "rank" in acceptor_rows.columns:
                acceptor_rows = acceptor_rows.rename(columns={"rank": "acceptor_rank"})
                acceptor_merge_cols.append("acceptor_rank")

            donor_rows = donor_rows[donor_merge_cols].drop_duplicates(subset=["species"])
            acceptor_rows = acceptor_rows[acceptor_merge_cols].drop_duplicates(subset=["species"])

            master = master.merge(donor_rows, on="species", how="left")
            master = master.merge(acceptor_rows, on="species", how="left")

    network = read_csv_if_exists(network_fp)
    if network is not None:
        network = ensure_species_col(network, ["species"])
        network = network[network["species"].isin(primary_set)].copy()
        keep = [c for c in [
            "species",
            "out_strength",
            "in_strength",
            "out_degree",
            "in_degree",
            "out_degree_centrality",
            "in_degree_centrality",
            "betweenness_centrality",
            "eigenvector_centrality",
            "mean_pairwise_asymmetry",
            "cluster",
            "network_community",
        ] if c in network.columns]
        master = master.merge(network[keep], on="species", how="left")

    cluster = read_csv_if_exists(cluster_fp)
    if cluster is not None:
        cluster = ensure_species_col(cluster, ["species", "species_name"])
        cluster = cluster[cluster["species"].isin(primary_set)].copy()
        rename_map = {}
        if "cluster" in cluster.columns and "cluster" in master.columns:
            rename_map["cluster"] = "cluster_k6"
        if "clustered order" in cluster.columns:
            rename_map["clustered order"] = "clustered_order"
        cluster = cluster.rename(columns=rename_map)
        keep = [c for c in ["species", "cluster", "cluster_k6", "clustered_order"] if c in cluster.columns]
        cluster = cluster[keep].drop_duplicates(subset=["species"])
        master = master.merge(cluster, on="species", how="left")

    community = read_csv_if_exists(community_fp)
    if community is not None:
        community = ensure_species_col(community, ["species", "species_name"])
        community = community[community["species"].isin(primary_set)].copy()
        rename_map = {}
        if "network_community" in community.columns and "network_community" in master.columns:
            rename_map["network_community"] = "network_community_table"
        elif "community" in community.columns:
            rename_map["community"] = "network_community_table"
        community = community.rename(columns=rename_map)
        keep = [c for c in ["species", "network_community", "network_community_table"] if c in community.columns]
        community = community[keep].drop_duplicates(subset=["species"])
        master = master.merge(community, on="species", how="left")

    reposition = read_csv_if_exists(reposition_fp)
    if reposition is not None:
        reposition = ensure_species_col(reposition, ["species"])
        reposition = reposition[reposition["species"].isin(primary_set)].copy()
        keep = [c for c in [
            "species",
            "full_geofresh_outgoing_mean",
            "climate_local_outgoing_mean",
            "shift_full_minus_climate",
            "abs_shift",
        ] if c in reposition.columns]
        reposition = reposition[keep].drop_duplicates(subset=["species"])
        master = master.merge(reposition, on="species", how="left")

    residual = read_csv_if_exists(residual_fp)
    if residual is not None:
        residual = ensure_species_col(residual, ["species"])
        residual = residual[residual["species"].isin(primary_set)].copy()
        keep = [c for c in ["species", "mean_abs_pairwise_residual"] if c in residual.columns]
        residual = residual[keep].drop_duplicates(subset=["species"])
        master = master.merge(residual, on="species", how="left")

    bootstrap = read_csv_if_exists(bootstrap_fp)
    if bootstrap is not None:
        bootstrap = ensure_species_col(bootstrap, ["species"])
        bootstrap = bootstrap[bootstrap["species"].isin(primary_set)].copy()
        keep = [c for c in [
            "species",
            "donor_top20_frequency",
            "acceptor_top20_frequency",
        ] if c in bootstrap.columns]
        bootstrap = bootstrap[keep].drop_duplicates(subset=["species"])
        master = master.merge(bootstrap, on="species", how="left")

    if "invasion_status" not in master.columns:
        master["invasion_status"] = pd.NA
    if "is_invasive" not in master.columns:
        master["is_invasive"] = pd.NA

    # preserve primary matrix species order
    master["species"] = pd.Categorical(master["species"], categories=primary_species, ordered=True)
    master = master.sort_values("species").reset_index(drop=True)
    master["species"] = master["species"].astype(str)

    preferred_order = [
        "species",
        "family",
        "genus",
        "modal_continent",
        "n_countries",
        "n_continents",
        "range_bbox_area_km2",
        "n_basins",
        "total_records",
        "records_after_thinning",
        "native_only_flag",
        "ecological_strategy",
        "invasion_status",
        "is_invasive",
        "donor_score",
        "acceptor_score",
        "donor_rank",
        "acceptor_rank",
        "out_strength",
        "in_strength",
        "out_degree",
        "in_degree",
        "out_degree_centrality",
        "in_degree_centrality",
        "betweenness_centrality",
        "eigenvector_centrality",
        "mean_pairwise_asymmetry",
        "cluster",
        "cluster_k6",
        "clustered_order",
        "network_community",
        "network_community_table",
        "full_geofresh_outgoing_mean",
        "climate_local_outgoing_mean",
        "shift_full_minus_climate",
        "abs_shift",
        "mean_abs_pairwise_residual",
        "donor_top20_frequency",
        "acceptor_top20_frequency",
        "temporal_span_years",
        "temporal_density",
        "records_high_accuracy",
        "records_snapped_le_200m",
        "records_deduplicated_segment",
    ]
    cols = [c for c in preferred_order if c in master.columns] + [c for c in master.columns if c not in preferred_order]
    master = master[cols]

    master.to_csv(OUT_CSV, index=False)

    summary = {
        "n_species_primary_matrix": int(len(primary_species)),
        "n_species_table": int(len(master)),
        "all_primary_species_matched": bool(len(master) == len(primary_species)),
        "columns": list(master.columns),
        "n_missing_donor_score": int(master["donor_score"].isna().sum()) if "donor_score" in master.columns else None,
        "n_missing_acceptor_score": int(master["acceptor_score"].isna().sum()) if "acceptor_score" in master.columns else None,
        "n_missing_network_metrics": int(master["out_strength"].isna().sum()) if "out_strength" in master.columns else None,
        "n_missing_bootstrap_freq": int(master["donor_top20_frequency"].isna().sum()) if "donor_top20_frequency" in master.columns else None,
        "output_csv": str(OUT_CSV),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_CSV}")
    print(f"  {OUT_JSON}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nFirst 10 rows:")
    print(master.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
