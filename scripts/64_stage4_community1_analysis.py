#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
OUTDIR = STAGE3 / "stage4_community1"
OUTDIR.mkdir(parents=True, exist_ok=True)

MATRIX_FP = STAGE3 / "matrices_full_geofresh_gbm.h5"
NET_FP = STAGE3 / "stage4b" / "table_stage4b_species_network_metrics.csv"
CLUSTER_FP = STAGE3 / "stage4a" / "table_stage4a_cluster_membership_k6.csv"
INV_FP = ROOT / "data" / "processed" / "species_inventory.csv"
ASYM_FP = STAGE3 / "tables" / "table_stage3_top20_asymmetric_pairs.csv"


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def plot_submatrix(species: list[str], mat: np.ndarray, out_fp: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, aspect="auto")

    ax.set_xticks(np.arange(len(species)))
    ax.set_yticks(np.arange(len(species)))
    ax.set_xticklabels(species, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(species, fontsize=8)
    ax.set_title("Community 1 within-community mean compatibility")

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=7)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mean suitability")
    fig.tight_layout()
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_reach_into_clusters(df: pd.DataFrame, cluster_cols: list[str], out_fp: Path) -> None:
    species = df["species"].tolist()
    vals = df[cluster_cols].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(len(species))

    for col in cluster_cols:
        ax.bar(species, df[col].to_numpy(), bottom=bottom, label=col)
        bottom += df[col].to_numpy()

    ax.set_ylabel("Mean outgoing compatibility")
    ax.set_title("Community 1 non-Cambaridae reach into Cambaridae clusters")
    ax.set_xticklabels(species, rotation=45, ha="right")
    ax.legend(frameon=False, title="Target cluster")
    fig.tight_layout()
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    species_order, mat = load_matrix(MATRIX_FP)

    net = pd.read_csv(NET_FP)
    cl = pd.read_csv(CLUSTER_FP)
    inv = pd.read_csv(INV_FP)

    df = (
        net.merge(cl[["species", "cluster"]], on="species", how="left", suffixes=("", "_clusterfile"))
           .merge(inv[["species_name", "family", "genus", "modal_continent"]], left_on="species", right_on="species_name", how="left")
    )

    # Community 1 members
    community1 = df[df["network_community"] == 1].copy().sort_values("species")
    community1_species = community1["species"].tolist()

    # Non-Cambaridae subset inside Community 1
    non_cambaridae = community1[community1["family"] != "Cambaridae"].copy().sort_values("species")
    non_c_species = non_cambaridae["species"].tolist()

    # 1) within-community 10x10 matrix
    idx = [species_order.index(sp) for sp in community1_species]
    submat = mat[np.ix_(idx, idx)]
    submat_df = pd.DataFrame(submat, index=community1_species, columns=community1_species)
    submat_df.to_csv(OUTDIR / "table_stage4_community1_within_matrix.csv")

    plot_submatrix(
        community1_species,
        submat,
        OUTDIR / "fig_stage4_community1_within_matrix.png",
    )

    # 2) reach into Cambaridae clusters
    camb = df[df["family"] == "Cambaridae"].copy()
    camb_clusters = sorted(camb["cluster"].dropna().unique().tolist())

    rows = []
    for sp in non_c_species:
        i = species_order.index(sp)
        row = {"species": sp}
        for k in camb_clusters:
            target_species = camb[camb["cluster"] == k]["species"].tolist()
            target_idx = [species_order.index(x) for x in target_species]
            row[f"cluster_{int(k)}"] = float(np.nanmean(mat[i, target_idx]))
        rows.append(row)

    reach_df = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    reach_df.to_csv(OUTDIR / "table_stage4_community1_reach_into_cambaridae_clusters.csv", index=False)

    cluster_cols = [c for c in reach_df.columns if c.startswith("cluster_")]
    plot_reach_into_clusters(
        reach_df,
        cluster_cols,
        OUTDIR / "fig_stage4_community1_reach_into_clusters.png",
    )

    # 3) cross-continental asymmetric cases
    asym = pd.read_csv(ASYM_FP)

    meta = inv[["species_name", "family", "genus", "modal_continent"]].rename(columns={"species_name": "species"})
    asym = (
        asym.merge(meta.add_suffix("_a"), left_on="species_a", right_on="species_a", how="left")
            .merge(meta.add_suffix("_b"), left_on="species_b", right_on="species_b", how="left")
    )

    cross_cont = asym[asym["modal_continent_a"] != asym["modal_continent_b"]].copy()
    cross_cont = cross_cont.sort_values("abs_diff", ascending=False).reset_index(drop=True)
    cross_cont.insert(0, "rank_cross_continental", np.arange(1, len(cross_cont) + 1))
    cross_cont.to_csv(OUTDIR / "table_stage4_community1_cross_continental_asym_pairs.csv", index=False)

    summary = {
        "community1_n_total": int(len(community1)),
        "community1_n_non_cambaridae": int(len(non_cambaridae)),
        "community1_species": community1_species,
        "community1_non_cambaridae_species": non_c_species,
        "cambaridae_clusters_present": [int(x) for x in camb_clusters],
        "n_cross_continental_pairs_in_top20": int(len(cross_cont)),
    }
    (OUTDIR / "stage4_community1_summary.json").write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUTDIR / 'table_stage4_community1_within_matrix.csv'}")
    print(f"  {OUTDIR / 'table_stage4_community1_reach_into_cambaridae_clusters.csv'}")
    print(f"  {OUTDIR / 'table_stage4_community1_cross_continental_asym_pairs.csv'}")
    print(f"  {OUTDIR / 'fig_stage4_community1_within_matrix.png'}")
    print(f"  {OUTDIR / 'fig_stage4_community1_reach_into_clusters.png'}")
    print(f"  {OUTDIR / 'stage4_community1_summary.json'}")

    print("\nCommunity 1 species:")
    for sp in community1_species:
        print(" ", sp)

    print("\nNon-Cambaridae species:")
    for sp in non_c_species:
        print(" ", sp)

    print("\nCross-continental asymmetric pairs found in top20:")
    if len(cross_cont) == 0:
        print("  none in current top20 file")
    else:
        print(cross_cont[["species_a", "species_b", "modal_continent_a", "modal_continent_b", "abs_diff"]].to_string(index=False))


if __name__ == "__main__":
    main()