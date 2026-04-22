#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "data" / "processed" / "stage3"
TABLES = STAGE3 / "tables"
FIGURES = STAGE3 / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


PRIMARY_MATRIX = STAGE3 / "matrices_full_geofresh_gbm.h5"
CLIMATE_COMMON155 = STAGE3 / "matrices_climate_local_gbm_common155.h5"


def load_matrix_bundle(path: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(path, "r") as h5:
        species = [s.decode() if isinstance(s, bytes) else str(s) for s in h5["species"][:]]
        mat = h5["mean_suitability"][:].astype(np.float64)
    return species, mat


def offdiag_row_col_means(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = mat.shape[0]
    row_means = np.empty(n, dtype=np.float64)
    col_means = np.empty(n, dtype=np.float64)

    for i in range(n):
        row_means[i] = np.nanmean(np.delete(mat[i, :], i))
        col_means[i] = np.nanmean(np.delete(mat[:, i], i))

    return row_means, col_means


def top_donors_acceptors(species: list[str], mat: np.ndarray, top_n: int = 20) -> tuple[pd.DataFrame, pd.DataFrame]:
    row_means, col_means = offdiag_row_col_means(mat)

    donor_order = np.argsort(-row_means)
    acceptor_order = np.argsort(-col_means)

    donors = pd.DataFrame({
        "rank": np.arange(1, top_n + 1),
        "species": np.array(species)[donor_order][:top_n],
        "score": row_means[donor_order][:top_n],
        "role": "donor",
    })

    acceptors = pd.DataFrame({
        "rank": np.arange(1, top_n + 1),
        "species": np.array(species)[acceptor_order][:top_n],
        "score": col_means[acceptor_order][:top_n],
        "role": "acceptor",
    })

    return donors, acceptors


def top_asymmetric_pairs(species: list[str], mat: np.ndarray, top_n: int = 20) -> pd.DataFrame:
    n = len(species)
    rows = []

    for i in range(n):
        for j in range(i + 1, n):
            a_to_b = float(mat[i, j])
            b_to_a = float(mat[j, i])
            rows.append({
                "species_a": species[i],
                "species_b": species[j],
                "a_to_b": a_to_b,
                "b_to_a": b_to_a,
                "abs_diff": abs(a_to_b - b_to_a),
                "signed_diff_a_minus_b": a_to_b - b_to_a,
            })

    df = (
        pd.DataFrame(rows)
        .sort_values("abs_diff", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def top_repositioning(
    species_fg: list[str],
    mat_fg: np.ndarray,
    species_cl: list[str],
    mat_cl: np.ndarray,
    top_n: int = 20,
) -> pd.DataFrame:
    if species_fg != species_cl:
        raise ValueError("Species rosters do not match between full_geofresh and climate_local common155 matrices.")

    fg_out, _ = offdiag_row_col_means(mat_fg)
    cl_out, _ = offdiag_row_col_means(mat_cl)

    df = pd.DataFrame({
        "species": species_fg,
        "full_geofresh_outgoing_mean": fg_out,
        "climate_local_outgoing_mean": cl_out,
        "shift_full_minus_climate": fg_out - cl_out,
        "abs_shift": np.abs(fg_out - cl_out),
    }).sort_values("abs_shift", ascending=False).head(top_n).reset_index(drop=True)

    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def write_workbook_tables(
    donors: pd.DataFrame,
    acceptors: pd.DataFrame,
    asym_pairs: pd.DataFrame,
    reposition: pd.DataFrame,
) -> None:
    donors_acceptors = pd.concat([donors, acceptors], ignore_index=True)
    donors_acceptors.to_csv(TABLES / "table_stage3_top20_donors_acceptors.csv", index=False)
    asym_pairs.to_csv(TABLES / "table_stage3_top20_asymmetric_pairs.csv", index=False)
    reposition.to_csv(TABLES / "table_stage3_top20_repositioning.csv", index=False)

    with pd.ExcelWriter(TABLES / "stage3_compact_workbook.xlsx", engine="openpyxl") as writer:
        donors.to_excel(writer, sheet_name="top20_donors", index=False)
        acceptors.to_excel(writer, sheet_name="top20_acceptors", index=False)
        asym_pairs.to_excel(writer, sheet_name="top20_asymmetric_pairs", index=False)
        reposition.to_excel(writer, sheet_name="top20_repositioning", index=False)


def short_name(name: str, max_len: int = 28) -> str:
    return name if len(name) <= max_len else name[: max_len - 1] + "…"


def plot_top10_asymmetric_pairs(asym_pairs: pd.DataFrame, outpath: Path) -> None:
    top10 = asym_pairs.head(10).copy()

    fig, ax = plt.subplots(figsize=(12, 7))
    y_positions = np.arange(len(top10))[::-1]
    x_max = float(np.nanmax(top10[["a_to_b", "b_to_a"]].to_numpy())) * 1.05

    for y, (_, row) in zip(y_positions, top10.iterrows()):
        a_to_b = row["a_to_b"]
        b_to_a = row["b_to_a"]

        left = min(a_to_b, b_to_a)
        right = max(a_to_b, b_to_a)

        ax.hlines(y=y, xmin=left, xmax=right, linewidth=2)

        ax.plot(a_to_b, y, marker="o", markersize=8)
        ax.plot(b_to_a, y, marker="s", markersize=8)

        label = (
            f"{row['rank']}. {short_name(row['species_a'], 22)} → "
            f"{short_name(row['species_b'], 22)}   "
            f"{a_to_b:.3f} vs {b_to_a:.3f}   Δ={row['abs_diff']:.3f}"
        )
        ax.text(0.01, y + 0.18, label, fontsize=9)

    ax.set_xlim(0, x_max)
    ax.set_ylim(-1, len(top10))
    ax.set_yticks([])
    ax.set_xlabel("Mean suitability")
    ax.set_title("Top 10 most asymmetric species pairs")
    ax.grid(True, axis="x", alpha=0.3)

    ax.plot([], [], marker="o", linestyle="None", label="A → B")
    ax.plot([], [], marker="s", linestyle="None", label="B → A")
    ax.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

def extract_submatrix(species: list[str], mat: np.ndarray, wanted: Iterable[str]) -> tuple[list[str], np.ndarray]:
    wanted = list(wanted)
    indices = [species.index(sp) for sp in wanted]
    sub = mat[np.ix_(indices, indices)]
    return wanted, sub


def plot_zoomed_matrix(labels: list[str], submat: np.ndarray, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(submat, interpolation="nearest", aspect="auto")

    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    for i in range(submat.shape[0]):
        for j in range(submat.shape[1]):
            ax.text(
                j,
                i,
                f"{submat[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=9,
            )

    ax.set_title("Procambarus fallax complex subcluster\n(full_geofresh / GBM mean suitability)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mean suitability")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    fg_species, fg_mat = load_matrix_bundle(PRIMARY_MATRIX)
    cl_species, cl_mat = load_matrix_bundle(CLIMATE_COMMON155)

    donors, acceptors = top_donors_acceptors(fg_species, fg_mat, top_n=20)
    asym_pairs = top_asymmetric_pairs(fg_species, fg_mat, top_n=20)
    reposition = top_repositioning(fg_species, fg_mat, cl_species, cl_mat, top_n=20)

    write_workbook_tables(donors, acceptors, asym_pairs, reposition)

    plot_top10_asymmetric_pairs(
        asym_pairs,
        FIGURES / "fig_stage3_top10_asymmetric_pairs.png",
    )

    subcluster_species = [
        "Procambarus fallax",
        "Procambarus pearsei",
        "Procambarus alleni",
        "Procambarus plumimanus",
    ]
    labels, submat = extract_submatrix(fg_species, fg_mat, subcluster_species)
    plot_zoomed_matrix(
        labels,
        submat,
        FIGURES / "fig_stage3_procambarus_fallax_complex_zoom.png",
    )

    print("Wrote tables:")
    print(f"  {TABLES / 'table_stage3_top20_donors_acceptors.csv'}")
    print(f"  {TABLES / 'table_stage3_top20_asymmetric_pairs.csv'}")
    print(f"  {TABLES / 'table_stage3_top20_repositioning.csv'}")
    print(f"  {TABLES / 'stage3_compact_workbook.xlsx'}")

    print("\nWrote figures:")
    print(f"  {FIGURES / 'fig_stage3_top10_asymmetric_pairs.png'}")
    print(f"  {FIGURES / 'fig_stage3_procambarus_fallax_complex_zoom.png'}")

    print("\nTop 5 donors:")
    print(donors.head(5).to_string(index=False))
    print("\nTop 5 acceptors:")
    print(acceptors.head(5).to_string(index=False))
    print("\nTop 5 asymmetric pairs:")
    print(asym_pairs.head(5).to_string(index=False))
    print("\nTop 5 repositioning species:")
    print(reposition.head(5).to_string(index=False))


if __name__ == "__main__":
    main()