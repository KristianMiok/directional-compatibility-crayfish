#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "data" / "processed" / "stage6_invasion_followup"

GROUP_FP = INDIR / "invasion_fullcompat_group_summary_broad.csv"
MASTER_FP = INDIR / "species_master_table_labeled_fullcompat_broad.csv"

OUT_CSV = INDIR / "invasion_broad_key_features_summary.csv"
OUT_JSON = INDIR / "invasion_broad_key_features_summary.json"


KEY_FEATURES = [
    "n_countries",
    "n_continents",
    "range_bbox_area_km2",
    "n_basins",
    "records_after_thinning",
    "temporal_density",
    "mean_outgoing_schoener_D",
    "mean_outgoing_warren_I",
    "mean_incoming_schoener_D",
    "mean_incoming_warren_I",
    "mean_incoming_fraction_above",
    "mean_incoming_q90_suitability",
]


def main() -> None:
    group = pd.read_csv(GROUP_FP)
    master = pd.read_csv(MASTER_FP)

    out = group[group["feature"].isin(KEY_FEATURES)].copy()

    # Add ranks by p-value for readability
    out["perm_p_rank"] = out["perm_p_mean_diff"].rank(method="dense", ascending=True)

    # Keep ordered like KEY_FEATURES
    out["feature"] = pd.Categorical(out["feature"], categories=KEY_FEATURES, ordered=True)
    out = out.sort_values("feature").reset_index(drop=True)

    out.to_csv(OUT_CSV, index=False)

    inv_n = int((pd.to_numeric(master["is_invasive"], errors="coerce") == 1).sum())
    non_n = int((pd.to_numeric(master["is_invasive"], errors="coerce") == 0).sum())

    summary = {
        "scheme": "broad",
        "n_invasive": inv_n,
        "n_noninvasive": non_n,
        "n_features": int(len(out)),
        "features": [str(x) for x in out["feature"].astype(str).tolist()],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_CSV}")
    print(f"  {OUT_JSON}")

    print("\nTable preview:")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
