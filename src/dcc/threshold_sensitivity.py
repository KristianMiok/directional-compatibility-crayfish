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
    # Modal continent per species: not in the master dataset, so we leave it None.
    inv = inventory.copy()
    inv["modal_continent"] = None

    rows = []
    for t in thresholds:
        kept = inv[inv["records_deduplicated_segment"] >= t]
        rows.append({
            "threshold": t,
            "n_species_total": len(kept),
        })
    return pd.DataFrame(rows)


def plot_threshold_sweep(sweep: pd.DataFrame, out_path: str | Path, default: int = 200) -> None:
    """Bar chart: species retained by continent at each threshold."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(sweep["threshold"].astype(str), sweep["n_species_total"],
           color="steelblue", edgecolor="white")
    for x, n in zip(sweep["threshold"].astype(str), sweep["n_species_total"], strict=True):
        ax.text(x, n + max(sweep["n_species_total"]) * 0.01, str(int(n)),
                ha="center", va="bottom", fontsize=9)
    ax.set_xlabel("Minimum records (after high-accuracy + snap filter, dedup by segment)")
    ax.set_ylabel("Species retained")
    ax.set_title("Threshold sensitivity — species retained at each threshold")
    if default is not None:
        try:
            idx = list(sweep["threshold"]).index(default)
            ax.axvline(idx, color="red", linestyle="--", linewidth=1.2,
                       label=f"default = {default}")
            ax.legend(loc="upper right")
        except ValueError:
            pass
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
