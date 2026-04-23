#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


ROOT = Path(__file__).resolve().parents[1]
STAGE5 = ROOT / "data" / "processed" / "stage5"
FIGDIR = STAGE5 / "figures"


LEFT_FP = FIGDIR / "fig_stage4_community1_within_matrix.png"
RIGHT_FP = FIGDIR / "fig_stage4_community1_reach_into_clusters.png"
OUT_FP = FIGDIR / "Fig5_community1_panel.png"


def main() -> None:
    if not LEFT_FP.exists():
        raise FileNotFoundError(f"Missing left panel: {LEFT_FP}")
    if not RIGHT_FP.exists():
        raise FileNotFoundError(f"Missing right panel: {RIGHT_FP}")

    left_img = mpimg.imread(LEFT_FP)
    right_img = mpimg.imread(RIGHT_FP)

    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.15], wspace=0.08)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    ax1.imshow(left_img)
    ax1.axis("off")
    ax1.set_title("A", loc="left", fontweight="bold", fontsize=16, pad=8)

    ax2.imshow(right_img)
    ax2.axis("off")
    ax2.set_title("B", loc="left", fontweight="bold", fontsize=16, pad=8)


    fig.savefig(OUT_FP, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(f"  {OUT_FP}")


if __name__ == "__main__":
    main()