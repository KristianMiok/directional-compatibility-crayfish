#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community import louvain_communities


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
STAGE4A = STAGE3 / "stage4a"
OUTDIR = STAGE3 / "stage4b"
OUTDIR.mkdir(parents=True, exist_ok=True)

MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"
CLUSTER_FP = STAGE4A / "table_stage4a_cluster_membership_k6.csv"


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def reciprocity_per_species(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        diffs = np.abs(np.delete(mat[i, :] - mat[:, i], i))
        out[i] = np.mean(diffs)
    return out


def build_directed_graph(
    species: list[str],
    mat: np.ndarray,
    threshold: float = 0.15,
) -> nx.DiGraph:
    G = nx.DiGraph()
    for sp in species:
        G.add_node(sp)

    n = len(species)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            w = float(mat[i, j])
            if w >= threshold:
                G.add_edge(species[i], species[j], weight=w)
    return G


def make_metrics_table(species: list[str], mat: np.ndarray, G: nx.DiGraph) -> pd.DataFrame:
    out_strength = mat.sum(axis=1) - np.diag(mat)
    in_strength = mat.sum(axis=0) - np.diag(mat)

    out_degree = pd.Series(dict(G.out_degree()), name="out_degree")
    in_degree = pd.Series(dict(G.in_degree()), name="in_degree")
    out_degree_c = pd.Series(nx.out_degree_centrality(G), name="out_degree_centrality")
    in_degree_c = pd.Series(nx.in_degree_centrality(G), name="in_degree_centrality")

    # Weighted betweenness on inverse weights so stronger edges = shorter paths
    H = G.copy()
    for u, v, d in H.edges(data=True):
        d["inv_weight"] = 1.0 / max(d["weight"], 1e-9)

    betweenness = pd.Series(
        nx.betweenness_centrality(H, weight="inv_weight", normalized=True),
        name="betweenness_centrality",
    )

    # Eigenvector centrality on weighted directed graph
    try:
        eigen = pd.Series(
            nx.eigenvector_centrality_numpy(G, weight="weight"),
            name="eigenvector_centrality",
        )
    except Exception:
        eigen = pd.Series({sp: np.nan for sp in species}, name="eigenvector_centrality")

    reciprocity = pd.Series(
        reciprocity_per_species(mat),
        index=species,
        name="mean_pairwise_asymmetry",
    )

    df = pd.DataFrame({"species": species})
    df["out_strength"] = [out_strength[species.index(sp)] for sp in species]
    df["in_strength"] = [in_strength[species.index(sp)] for sp in species]

    for s in [out_degree, in_degree, out_degree_c, in_degree_c, betweenness, eigen, reciprocity]:
        df = df.merge(s.rename_axis("species").reset_index(), on="species", how="left")

    return df


def community_detection(G: nx.DiGraph) -> pd.DataFrame:
    # Louvain in networkx is for undirected graphs, so use symmetrized weighted graph
    UG = nx.Graph()
    for u, v, d in G.edges(data=True):
        w = d["weight"]
        if UG.has_edge(u, v):
            UG[u][v]["weight"] += w
        else:
            UG.add_edge(u, v, weight=w)

    comms = louvain_communities(UG, weight="weight", seed=42, resolution=1.0)
    rows = []
    for cid, members in enumerate(comms, start=1):
        for sp in members:
            rows.append({"species": sp, "network_community": cid})
    return pd.DataFrame(rows).sort_values(["network_community", "species"]).reset_index(drop=True)


def agreement_rate(a: pd.Series, b: pd.Series) -> float:
    merged = pd.DataFrame({"a": a, "b": b}).dropna()
    if len(merged) == 0:
        return float("nan")
    # crude same-label agreement only if labels share numbering; mostly informative placeholder
    return float((merged["a"] == merged["b"]).mean())


def community_level_graph(G: nx.DiGraph, community_df: pd.DataFrame) -> nx.DiGraph:
    cmap = dict(zip(community_df["species"], community_df["network_community"]))
    CG = nx.DiGraph()

    for c in sorted(community_df["network_community"].unique()):
        CG.add_node(c)

    for u, v, d in G.edges(data=True):
        cu = cmap[u]
        cv = cmap[v]
        if CG.has_edge(cu, cv):
            CG[cu][cv]["weight"] += d["weight"]
            CG[cu][cv]["n_edges"] += 1
        else:
            CG.add_edge(cu, cv, weight=d["weight"], n_edges=1)
    return CG


def plot_community_graph(CG: nx.DiGraph, outpath: Path) -> None:
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(CG, seed=42, weight="weight")

    weights = np.array([d["weight"] for _, _, d in CG.edges(data=True)], dtype=float)
    if len(weights) == 0:
        weights = np.array([1.0])

    ptp_val = np.ptp(weights)
    widths = 1.0 + 5.0 * (weights - weights.min()) / (ptp_val if ptp_val > 0 else 1.0)

    nx.draw_networkx_nodes(CG, pos, node_size=1600)
    nx.draw_networkx_labels(CG, pos, labels={n: f"C{n}" for n in CG.nodes()}, font_size=11)
    nx.draw_networkx_edges(
        CG,
        pos,
        width=widths,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=18,
        connectionstyle="arc3,rad=0.08",
    )

    edge_labels = {(u, v): f"{d['weight']:.1f}" for u, v, d in CG.edges(data=True)}
    nx.draw_networkx_edge_labels(CG, pos, edge_labels=edge_labels, font_size=8)

    plt.title("Community-level directed compatibility graph")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    species, mat = load_matrix(MATRIX_FP)
    clusters = pd.read_csv(CLUSTER_FP)

    print("Building directed graph...")
    G = build_directed_graph(species, mat, threshold=0.15)
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges retained (weight >= 0.15): {G.number_of_edges()}")

    print("Computing per-species network metrics...")
    metrics = make_metrics_table(species, mat, G)

    print("Running community detection...")
    comm_df = community_detection(G)

    out = metrics.merge(clusters[["species", "cluster"]], on="species", how="left")
    out = out.merge(comm_df, on="species", how="left")

    out.to_csv(OUTDIR / "table_stage4b_species_network_metrics.csv", index=False)
    comm_df.to_csv(OUTDIR / "table_stage4b_network_communities.csv", index=False)

    # Compare network communities to stage4a clusters
    merged_compare = out[["species", "cluster", "network_community"]].dropna()
    community_cluster_crosstab = pd.crosstab(
        merged_compare["cluster"],
        merged_compare["network_community"],
        rownames=["cluster_k6"],
        colnames=["network_community"],
    )
    community_cluster_crosstab.to_csv(OUTDIR / "table_stage4b_cluster_vs_community_crosstab.csv")

    agreement = agreement_rate(merged_compare["cluster"], merged_compare["network_community"])

    CG = community_level_graph(G, comm_df)
    plot_community_graph(CG, OUTDIR / "fig_stage4b_community_graph.png")

    summary = {
        "n_species": len(species),
        "graph_threshold": 0.15,
        "n_edges_retained": int(G.number_of_edges()),
        "n_network_communities": int(comm_df["network_community"].nunique()),
        "cluster_vs_community_exact_label_agreement": agreement,
        "outputs": [
            "table_stage4b_species_network_metrics.csv",
            "table_stage4b_network_communities.csv",
            "table_stage4b_cluster_vs_community_crosstab.csv",
            "fig_stage4b_community_graph.png",
        ],
    }
    (OUTDIR / "stage4b_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {OUTDIR / 'table_stage4b_species_network_metrics.csv'}")
    print(f"  {OUTDIR / 'table_stage4b_network_communities.csv'}")
    print(f"  {OUTDIR / 'table_stage4b_cluster_vs_community_crosstab.csv'}")
    print(f"  {OUTDIR / 'fig_stage4b_community_graph.png'}")
    print(f"  {OUTDIR / 'stage4b_summary.json'}")

    print("\nTop 10 by out_strength:")
    print(out.sort_values("out_strength", ascending=False).head(10)[
        ["species", "out_strength", "in_strength", "betweenness_centrality", "eigenvector_centrality", "mean_pairwise_asymmetry", "cluster", "network_community"]
    ].to_string(index=False))

    print("\nTop 10 by in_strength:")
    print(out.sort_values("in_strength", ascending=False).head(10)[
        ["species", "out_strength", "in_strength", "betweenness_centrality", "eigenvector_centrality", "mean_pairwise_asymmetry", "cluster", "network_community"]
    ].to_string(index=False))

    print("\nCluster vs network community crosstab:")
    print(community_cluster_crosstab.to_string())


if __name__ == "__main__":
    main()