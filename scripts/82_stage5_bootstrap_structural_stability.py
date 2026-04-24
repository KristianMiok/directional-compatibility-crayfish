#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
import networkx as nx
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
BOOTROOT = STAGE3 / "stage4d"
OUTDIR = ROOT / "data" / "processed" / "stage5" / "supplementary_analyses"
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIMARY_MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"
PRIMARY_CLUSTER_FP = STAGE3 / "stage4a" / "table_stage4a_cluster_membership_k6.csv"
PRIMARY_COMMUNITY_FP = STAGE3 / "stage4b" / "table_stage4b_network_communities.csv"

EDGE_THRESHOLD = 0.15
K_CLUSTERS = 6


def load_species_and_matrix(h5_fp: Path, dataset_name: str = "mean_suitability") -> tuple[list[str], np.ndarray]:
    with h5py.File(h5_fp, "r") as h5:
        species = [
            s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)
            for s in h5["species"][:]
        ]
        if dataset_name not in h5:
            raise KeyError(f"Dataset '{dataset_name}' not found in {h5_fp}")
        mat = np.array(h5[dataset_name][:], dtype=float)
    return species, mat


def compute_cluster_labels(species: list[str], mat: np.ndarray, k: int = 6) -> pd.Series:
    sym = (mat + mat.T) / 2.0
    dist = 1.0 - sym
    dist = np.clip(dist, 0.0, None)
    np.fill_diagonal(dist, 0.0)

    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="ward")
    labels = fcluster(Z, t=k, criterion="maxclust")

    return pd.Series(labels.astype(int), index=species, name="cluster")


def build_directed_graph(species: list[str], mat: np.ndarray, edge_threshold: float = 0.15) -> nx.DiGraph:
    G = nx.DiGraph()
    G.add_nodes_from(species)

    n = len(species)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            w = float(mat[i, j])
            if np.isfinite(w) and w >= edge_threshold:
                G.add_edge(species[i], species[j], weight=w)
    return G


def compute_communities(species: list[str], mat: np.ndarray, edge_threshold: float = 0.15) -> pd.Series:
    G = build_directed_graph(species, mat, edge_threshold=edge_threshold)

    # Use undirected weighted graph for Louvain-like community detection,
    # aggregating reciprocal directed weights.
    UG = nx.Graph()
    UG.add_nodes_from(species)

    for u, v, data in G.edges(data=True):
        w = float(data["weight"])
        if UG.has_edge(u, v):
            UG[u][v]["weight"] += w
        else:
            UG.add_edge(u, v, weight=w)

    communities = nx.community.louvain_communities(UG, weight="weight", seed=42)

    membership = {}
    for cid, comm in enumerate(communities, start=1):
        for sp in comm:
            membership[sp] = cid

    return pd.Series([membership[sp] for sp in species], index=species, name="community")


def load_primary_cluster_partition() -> pd.Series:
    df = pd.read_csv(PRIMARY_CLUSTER_FP)
    species_col = "species" if "species" in df.columns else "species_name"
    cluster_col = "cluster" if "cluster" in df.columns else "cluster_k6"
    return pd.Series(df[cluster_col].astype(int).to_numpy(), index=df[species_col].astype(str), name="cluster")


def load_primary_community_partition() -> pd.Series:
    df = pd.read_csv(PRIMARY_COMMUNITY_FP)
    species_col = "species" if "species" in df.columns else "species_name"
    community_col = "network_community" if "network_community" in df.columns else "community"
    return pd.Series(df[community_col].astype(int).to_numpy(), index=df[species_col].astype(str), name="community")


def main() -> None:
    primary_species, primary_mat = load_species_and_matrix(PRIMARY_MATRIX_FP, dataset_name="mean_suitability")
    primary_cluster = load_primary_cluster_partition().loc[primary_species]
    primary_community = load_primary_community_partition().loc[primary_species]

    rows = []

    for i in range(1, 101):
        iter_name = f"iter_{i:03d}"
        h5_fp = BOOTROOT / iter_name / "matrices_full_geofresh_gbm.h5"
        if not h5_fp.exists():
            print(f"Skipping {iter_name}: missing matrix")
            continue

        species, mat = load_species_and_matrix(h5_fp, dataset_name="mean_suitability")
        if species != primary_species:
            raise ValueError(f"Species order mismatch in {h5_fp}")

        cluster_boot = compute_cluster_labels(species, mat, k=K_CLUSTERS)
        community_boot = compute_communities(species, mat, edge_threshold=EDGE_THRESHOLD)

        cluster_ari = adjusted_rand_score(primary_cluster.to_numpy(), cluster_boot.to_numpy())
        cluster_nmi = normalized_mutual_info_score(primary_cluster.to_numpy(), cluster_boot.to_numpy())
        community_ari = adjusted_rand_score(primary_community.to_numpy(), community_boot.to_numpy())
        community_nmi = normalized_mutual_info_score(primary_community.to_numpy(), community_boot.to_numpy())
        community_count = int(community_boot.nunique())

        rows.append(
            {
                "iteration": iter_name,
                "cluster_ARI_vs_primary": float(cluster_ari),
                "cluster_NMI_vs_primary": float(cluster_nmi),
                "community_ARI_vs_primary": float(community_ari),
                "community_NMI_vs_primary": float(community_nmi),
                "community_count": community_count,
            }
        )

        print(
            f"{iter_name}: "
            f"cluster ARI={cluster_ari:.3f}, cluster NMI={cluster_nmi:.3f}, "
            f"community ARI={community_ari:.3f}, community NMI={community_nmi:.3f}, "
            f"community_count={community_count}"
        )

    if not rows:
        raise RuntimeError("No bootstrap matrices found.")

    out_df = pd.DataFrame(rows).sort_values("iteration").reset_index(drop=True)

    csv_fp = OUTDIR / "bootstrap_structural_stability.csv"
    json_fp = OUTDIR / "bootstrap_structural_stability_summary.json"

    out_df.to_csv(csv_fp, index=False)

    summary = {
        "n_iterations_used": int(len(out_df)),
        "cluster_ARI_mean": float(out_df["cluster_ARI_vs_primary"].mean()),
        "cluster_ARI_sd": float(out_df["cluster_ARI_vs_primary"].std(ddof=1)),
        "cluster_NMI_mean": float(out_df["cluster_NMI_vs_primary"].mean()),
        "cluster_NMI_sd": float(out_df["cluster_NMI_vs_primary"].std(ddof=1)),
        "community_ARI_mean": float(out_df["community_ARI_vs_primary"].mean()),
        "community_ARI_sd": float(out_df["community_ARI_vs_primary"].std(ddof=1)),
        "community_NMI_mean": float(out_df["community_NMI_vs_primary"].mean()),
        "community_NMI_sd": float(out_df["community_NMI_vs_primary"].std(ddof=1)),
        "community_count_mean": float(out_df["community_count"].mean()),
        "community_count_counts": {str(k): int(v) for k, v in out_df["community_count"].value_counts().sort_index().items()},
        "edge_threshold": float(EDGE_THRESHOLD),
        "k_clusters": int(K_CLUSTERS),
    }
    json_fp.write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {csv_fp}")
    print(f"  {json_fp}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
