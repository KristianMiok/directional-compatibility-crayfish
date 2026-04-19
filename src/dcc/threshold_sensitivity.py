"""Task 1.2 — threshold sensitivity.

Counts retained species at each minimum-records threshold, stratified by
continent and family. Plot drives Decision Point P1.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def threshold_sweep(
    inventory: pd.DataFrame,
    thresholds: list[int],
    woc: pd.DataFrame,
) -> pd.DataFrame:
    """For each threshold, count species retained, stratified by continent and family.

    Continent is taken as the *modal* continent for the species' records, which
    is robust to the few cross-continental introductions still present in WoC.
    """
    # Modal continent per species
    modal_cont = (
        woc.groupby("species_name")["continent"]
        .agg(lambda s: s.mode().iat[0] if not s.mode().empty else None)
        .rename("modal_continent")
    )
    inv = inventory.merge(modal_cont, left_on="species_name", right_index=True, how="left")

    rows = []
    for t in thresholds:
        kept = inv[inv["records_deduplicated_segment"] >= t]
        rows.append({
            "threshold": t,
            "n_species_total": len(kept),
            "by_continent": kept["modal_continent"].value_counts().to_dict(),
            "by_family": kept["family"].value_counts().to_dict(),
        })
    return pd.DataFrame(rows)


def plot_threshold_sweep(sweep: pd.DataFrame, out_path: str | Path, default: int = 200) -> None:
    """Bar chart: species retained by continent at each threshold."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Long-form for plotting
    records = []
    for _, row in sweep.iterrows():
        for cont, n in row["by_continent"].items():
            records.append({"threshold": row["threshold"], "continent": cont, "n_species": n})
    long = pd.DataFrame(records)

    if long.empty:
        # Nothing to plot — write an empty placeholder so the pipeline doesn't break.
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No species retained at any threshold", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    pivot = long.pivot(index="threshold", columns="continent", values="n_species").fillna(0)
    fig, ax = plt.subplots(figsize=(8, 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, edgecolor="white")
    ax.set_xlabel("Minimum records (deduplicated by segment)")
    ax.set_ylabel("Species retained")
    ax.set_title("Threshold sensitivity — species retained by continent")
    ax.axvline(
        list(pivot.index).index(default) if default in pivot.index else -1,
        color="red", linestyle="--", linewidth=1, label=f"default = {default}",
    )
    ax.legend(title="Continent", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
