#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, dendrogram
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
OUTDIR = STAGE3 / "stage4a"
OUTDIR.mkdir(parents=True, exist_ok=True)

MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def build_symmetric_similarity(mat: np.ndarray) -> np.ndarray:
    # For clustering we need a symmetric species-by-species similarity.
    # We average the two directions so the clustering reflects overall pair affinity.
    sim = 0.5 * (mat + mat.T)
    np.fill_diagonal(sim, 1.0)
    return sim


def build_distance_from_similarity(sim: np.ndarray) -> np.ndarray:
    dist = 1.0 - sim
    dist = np.clip(dist, 0.0, None)
    dist = 0.5 * (dist + dist.T)
    np.fill_diagonal(dist, 0.0)
    return dist


def cluster_species(dist: np.ndarray, method: str = "ward") -> tuple[np.ndarray, np.ndarray]:
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method=method)
    order = leaves_list(Z)
    return Z, order


def choose_cluster_counts(
    Z: np.ndarray,
    species: list[str],
    candidate_k: list[int] = [4, 5, 6, 7, 8]
) -> pd.DataFrame:
    rows = []
    n = len(species)
    for k in candidate_k:
        labels = fcluster(Z, t=k, criterion="maxclust")
        sizes = pd.Series(labels).value_counts().sort_index().tolist()
        rows.append({
            "k": k,
            "n_clusters_found": len(set(labels)),
            "min_cluster_size": int(min(sizes)),
            "max_cluster_size": int(max(sizes)),
            "cluster_sizes": ",".join(map(str, sizes)),
        })
    return pd.DataFrame(rows)


def make_cluster_membership_table(
    species: list[str],
    order: np.ndarray,
    labels: np.ndarray,
) -> pd.DataFrame:
    df = pd.DataFrame({
        "species": species,
        "cluster": labels,
        "clustered_order": np.arange(len(species)),
    })
    ordered_species = [species[i] for i in order]
    pos = {sp: i for i, sp in enumerate(ordered_species)}
    df["clustered_order"] = df["species"].map(pos)
    return df.sort_values(["cluster", "clustered_order"]).reset_index(drop=True)


def plot_clustered_heatmap(
    ordered_mat: np.ndarray,
    ordered_species: list[str],
    Z: np.ndarray,
    outpath: Path,
) -> None:
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[0.22, 1.0],
        height_ratios=[0.22, 1.0],
        wspace=0.02,
        hspace=0.02,
    )

    ax_empty = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_heat = fig.add_subplot(gs[1, 1])

    ax_empty.axis("off")

    dendrogram(
        Z,
        ax=ax_top,
        no_labels=True,
        color_threshold=None,
    )
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    for spine in ax_top.spines.values():
        spine.set_visible(False)

    dendrogram(
        Z,
        ax=ax_left,
        orientation="left",
        no_labels=True,
        color_threshold=None,
    )
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    for spine in ax_left.spines.values():
        spine.set_visible(False)

    im = ax_heat.imshow(ordered_mat, aspect="auto", interpolation="nearest")
    step = max(1, len(ordered_species) // 30)
    tick_idx = np.arange(0, len(ordered_species), step)

    ax_heat.set_xticks(tick_idx)
    ax_heat.set_yticks(tick_idx)
    ax_heat.set_xticklabels([ordered_species[i] for i in tick_idx], rotation=90, fontsize=7)
    ax_heat.set_yticklabels([ordered_species[i] for i in tick_idx], fontsize=7)

    ax_heat.set_xlabel("Target species (clustered order)")
    ax_heat.set_ylabel("Source species (clustered order)")

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02)
    cbar.set_label("Mean suitability")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_simple_heatmap(
    ordered_mat: np.ndarray,
    ordered_species: list[str],
    outpath: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(ordered_mat, aspect="auto", interpolation="nearest")

    step = max(1, len(ordered_species) // 30)
    tick_idx = np.arange(0, len(ordered_species), step)

    ax.set_xticks(tick_idx)
    ax.set_yticks(tick_idx)
    ax.set_xticklabels([ordered_species[i] for i in tick_idx], rotation=90, fontsize=7)
    ax.set_yticklabels([ordered_species[i] for i in tick_idx], fontsize=7)

    ax.set_xlabel("Target species (clustered order)")
    ax.set_ylabel("Source species (clustered order)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean suitability")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    species, mat = load_matrix(MATRIX_FP)

    sim = build_symmetric_similarity(mat)
    dist = build_distance_from_similarity(sim)

    print("Clustering species...")
    Z, order = cluster_species(dist, method="ward")
    ordered_species = [species[i] for i in order]
    ordered_mat = mat[np.ix_(order, order)]

    print("Evaluating candidate cluster counts...")
    cluster_count_df = choose_cluster_counts(Z, species, candidate_k=[4, 5, 6, 7, 8])
    cluster_count_df.to_csv(OUTDIR / "table_stage4a_candidate_cluster_counts.csv", index=False)

    # Choose one default cut for a first pass.
    default_k = 6
    labels = fcluster(Z, t=default_k, criterion="maxclust")
    membership = make_cluster_membership_table(species, order, labels)
    membership.to_csv(OUTDIR / "table_stage4a_cluster_membership_k6.csv", index=False)

    pd.DataFrame({
        "clustered_order": np.arange(len(ordered_species)),
        "species": ordered_species,
    }).to_csv(OUTDIR / "table_stage4a_species_order.csv", index=False)

    plot_clustered_heatmap(
        ordered_mat,
        ordered_species,
        Z,
        OUTDIR / "fig_stage4a_clustered_heatmap_with_dendrogram.png",
    )
    plot_simple_heatmap(
        ordered_mat,
        ordered_species,
        OUTDIR / "fig_stage4a_clustered_heatmap.png",
    )

    summary = {
        "n_species": len(species),
        "default_k": default_k,
        "distance_definition": "1 - 0.5*(M + M^T)",
        "linkage_method": "ward",
        "outputs": [
            "table_stage4a_candidate_cluster_counts.csv",
            "table_stage4a_cluster_membership_k6.csv",
            "table_stage4a_species_order.csv",
            "fig_stage4a_clustered_heatmap_with_dendrogram.png",
            "fig_stage4a_clustered_heatmap.png",
        ],
    }
    (OUTDIR / "stage4a_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    print(f"  {OUTDIR / 'table_stage4a_candidate_cluster_counts.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_cluster_membership_k6.csv'}")
    print(f"  {OUTDIR / 'table_stage4a_species_order.csv'}")
    print(f"  {OUTDIR / 'fig_stage4a_clustered_heatmap_with_dendrogram.png'}")
    print(f"  {OUTDIR / 'fig_stage4a_clustered_heatmap.png'}")
    print(f"  {OUTDIR / 'stage4a_summary.json'}")

    print("\nCandidate cluster counts:")
    print(cluster_count_df.to_string(index=False))

    print("\nFirst 20 species in clustered order:")
    for i, sp in enumerate(ordered_species[:20], start=1):
        print(f"  {i:2d}. {sp}")


if __name__ == "__main__":
    main()
