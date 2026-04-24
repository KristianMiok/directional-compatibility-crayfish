#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "data" / "processed" / "stage6_invasion_followup"
OUTDIR = INDIR
OUTDIR.mkdir(parents=True, exist_ok=True)

MASTER_FP = INDIR / "species_master_table_labeled_fullcompat_broad.csv"
OUT_FP = OUTDIR / "fig_invasion_broad_key_features.png"


PANELS = [
    ("mean_outgoing_schoener_D", "Outgoing Schoener's D"),
    ("mean_outgoing_warren_I", "Outgoing Warren's I"),
    ("mean_incoming_warren_I", "Incoming Warren's I"),
    ("n_continents", "Number of continents"),
]


def prep_groups(df: pd.DataFrame, col: str):
    inv = pd.to_numeric(df.loc[df["is_invasive"] == 1, col], errors="coerce").dropna()
    non = pd.to_numeric(df.loc[df["is_invasive"] == 0, col], errors="coerce").dropna()
    return inv, non


def main() -> None:
    df = pd.read_csv(MASTER_FP)
    df["is_invasive"] = pd.to_numeric(df["is_invasive"], errors="coerce")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    for ax, (col, title) in zip(axes, PANELS):
        inv, non = prep_groups(df, col)

        ax.boxplot(
            [non.values, inv.values],
            labels=["Non-invasive", "Invasive"],
            widths=0.6,
        )
        ax.set_title(title)
        ax.set_ylabel(col)

        # add points with jitter
        for xpos, series in zip([1, 2], [non, inv]):
            if len(series) > 0:
                jitter_x = [xpos] * len(series)
                ax.scatter(jitter_x, series.values, alpha=0.6, s=20)

        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_FP, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(f"  {OUT_FP}")


if __name__ == "__main__":
    main()
