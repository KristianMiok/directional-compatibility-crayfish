#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FP = ROOT / "data" / "processed" / "stage6_invasion_followup" / "species_master_table_labeled_fullcompat_broad.csv"
OUT_FP = ROOT / "data" / "processed" / "stage6_invasion_followup" / "fig_invasion_broad_outgoing_vs_incoming_warrenI.png"


def main() -> None:
    df = pd.read_csv(FP)
    df["is_invasive"] = pd.to_numeric(df["is_invasive"], errors="coerce")

    xcol = "mean_outgoing_warren_I"
    ycol = "mean_incoming_warren_I"

    plot_df = df[["species", "is_invasive", xcol, ycol]].copy()
    plot_df[xcol] = pd.to_numeric(plot_df[xcol], errors="coerce")
    plot_df[ycol] = pd.to_numeric(plot_df[ycol], errors="coerce")
    plot_df = plot_df.dropna(subset=[xcol, ycol, "is_invasive"]).copy()

    non = plot_df[plot_df["is_invasive"] == 0]
    inv = plot_df[plot_df["is_invasive"] == 1]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        non[xcol],
        non[ycol],
        alpha=0.7,
        s=35,
        label=f"Non-invasive (n={len(non)})",
    )

    ax.scatter(
        inv[xcol],
        inv[ycol],
        alpha=0.9,
        s=55,
        marker="D",
        label=f"Invasive (n={len(inv)})",
    )

    # annotate invasive species
    for _, row in inv.iterrows():
        ax.annotate(
            row["species"],
            (row[xcol], row[ycol]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    # reference diagonal
    min_val = min(plot_df[xcol].min(), plot_df[ycol].min())
    max_val = max(plot_df[xcol].max(), plot_df[ycol].max())
    ax.plot([min_val, max_val], [min_val, max_val], linestyle="--", linewidth=1)

    ax.set_xlabel("Mean outgoing Warren's I")
    ax.set_ylabel("Mean incoming Warren's I")
    ax.legend(frameon=True)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_FP, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(f"  {OUT_FP}")
    print("\nInvasive species coordinates:")
    print(inv[["species", xcol, ycol]].sort_values(xcol, ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
