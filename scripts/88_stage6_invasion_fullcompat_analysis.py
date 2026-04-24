#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "data" / "processed" / "stage6_invasion_followup"


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = 0
    lt = 0
    for xi in x:
        gt += np.sum(xi > y)
        lt += np.sum(xi < y)
    return float((gt - lt) / (len(x) * len(y)))


def permutation_pvalue_mean_diff(x: np.ndarray, y: np.ndarray, n_perm: int = 999, seed: int = 42) -> float:
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan

    obs = abs(np.mean(x) - np.mean(y))
    combined = np.concatenate([x, y])
    n_x = len(x)

    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(combined)
        diff = abs(np.mean(perm[:n_x]) - np.mean(perm[n_x:]))
        count += diff >= obs
    return float((count + 1) / (n_perm + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", choices=["conservative", "broad"], required=True)
    args = ap.parse_args()

    MASTER_FP = INDIR / f"species_master_table_labeled_fullcompat_{args.scheme}.csv"
    OUT_GROUP = INDIR / f"invasion_fullcompat_group_summary_{args.scheme}.csv"
    OUT_LOGIT = INDIR / f"invasion_fullcompat_logistic_coefficients_{args.scheme}.csv"
    OUT_JSON = INDIR / f"invasion_fullcompat_summary_{args.scheme}.json"

    df = pd.read_csv(MASTER_FP)
    df["is_invasive"] = pd.to_numeric(df["is_invasive"], errors="coerce")

    inv = df[df["is_invasive"] == 1].copy()
    non = df[df["is_invasive"] == 0].copy()

    group_features = [
        "mean_outgoing_mean_suitability",
        "mean_incoming_mean_suitability",
        "rank_outgoing_mean_suitability",
        "rank_incoming_mean_suitability",
        "mean_outgoing_fraction_above",
        "mean_incoming_fraction_above",
        "rank_outgoing_fraction_above",
        "rank_incoming_fraction_above",
        "mean_outgoing_schoener_D",
        "mean_incoming_schoener_D",
        "rank_outgoing_schoener_D",
        "rank_incoming_schoener_D",
        "mean_outgoing_warren_I",
        "mean_incoming_warren_I",
        "rank_outgoing_warren_I",
        "rank_incoming_warren_I",
        "mean_outgoing_q90_suitability",
        "mean_incoming_q90_suitability",
        "rank_outgoing_q90_suitability",
        "rank_incoming_q90_suitability",
        "out_strength",
        "in_strength",
        "betweenness_centrality",
        "eigenvector_centrality",
        "mean_pairwise_asymmetry",
        "shift_full_minus_climate",
        "abs_shift",
        "mean_abs_pairwise_residual",
        "donor_top20_frequency",
        "acceptor_top20_frequency",
        "n_countries",
        "n_continents",
        "range_bbox_area_km2",
        "n_basins",
        "records_after_thinning",
        "temporal_span_years",
        "temporal_density",
    ]
    group_features = [c for c in group_features if c in df.columns]

    rows = []
    for col in group_features:
        x = pd.to_numeric(inv[col], errors="coerce").to_numpy()
        y = pd.to_numeric(non[col], errors="coerce").to_numpy()
        rows.append({
            "feature": col,
            "n_invasive_nonmissing": int(np.isfinite(x).sum()),
            "n_noninvasive_nonmissing": int(np.isfinite(y).sum()),
            "mean_invasive": float(np.nanmean(x)) if np.isfinite(x).sum() else np.nan,
            "mean_noninvasive": float(np.nanmean(y)) if np.isfinite(y).sum() else np.nan,
            "median_invasive": float(np.nanmedian(x)) if np.isfinite(x).sum() else np.nan,
            "median_noninvasive": float(np.nanmedian(y)) if np.isfinite(y).sum() else np.nan,
            "mean_diff_invasive_minus_noninvasive": (
                float(np.nanmean(x) - np.nanmean(y))
                if np.isfinite(x).sum() and np.isfinite(y).sum()
                else np.nan
            ),
            "cliffs_delta": cliffs_delta(x, y),
            "perm_p_mean_diff": permutation_pvalue_mean_diff(x, y, n_perm=999, seed=42),
        })

    group_df = pd.DataFrame(rows).sort_values("perm_p_mean_diff", na_position="last")
    group_df.to_csv(OUT_GROUP, index=False)

    model_features = [
        "mean_outgoing_mean_suitability",
        "mean_incoming_mean_suitability",
        "mean_outgoing_fraction_above",
        "mean_incoming_fraction_above",
        "mean_outgoing_schoener_D",
        "mean_incoming_schoener_D",
        "mean_outgoing_warren_I",
        "mean_incoming_warren_I",
        "mean_outgoing_q90_suitability",
        "mean_incoming_q90_suitability",
        "out_strength",
        "in_strength",
        "betweenness_centrality",
        "eigenvector_centrality",
        "mean_pairwise_asymmetry",
        "shift_full_minus_climate",
        "abs_shift",
        "mean_abs_pairwise_residual",
        "donor_top20_frequency",
        "acceptor_top20_frequency",
        "n_countries",
        "n_continents",
        "range_bbox_area_km2",
        "n_basins",
        "records_after_thinning",
        "temporal_span_years",
        "temporal_density",
    ]
    model_features = [c for c in model_features if c in df.columns]

    model_df = df[["species", "is_invasive"] + model_features].dropna(subset=["is_invasive"]).copy()
    X = model_df[model_features].apply(pd.to_numeric, errors="coerce")
    y = model_df["is_invasive"].astype(int)

    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="l2",
            solver="liblinear",
            class_weight="balanced",
            random_state=42,
            max_iter=5000,
        )),
    ])
    pipe.fit(X, y)

    coefs = pipe.named_steps["clf"].coef_[0]
    coef_df = pd.DataFrame({
        "feature": model_features,
        "coefficient": coefs,
        "abs_coefficient": np.abs(coefs),
    }).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(OUT_LOGIT, index=False)

    summary = {
        "scheme": args.scheme,
        "n_species_total": int(len(df)),
        "n_invasive": int((df["is_invasive"] == 1).sum()),
        "n_noninvasive": int((df["is_invasive"] == 0).sum()),
        "top_group_differences": group_df.head(15)["feature"].tolist(),
        "top_logistic_features": coef_df.head(15)["feature"].tolist(),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_GROUP}")
    print(f"  {OUT_LOGIT}")
    print(f"  {OUT_JSON}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nTop group differences:")
    print(group_df.head(20).to_string(index=False))

    print("\nTop logistic coefficients:")
    print(coef_df.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
