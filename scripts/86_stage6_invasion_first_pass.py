#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "data" / "processed" / "stage6_invasion_followup"
OUTDIR = INDIR

MASTER_FP = INDIR / "species_master_table_base.csv"
LABELS_FP = INDIR / "species_invasion_labels_template.csv"

OUT_MASTER_LABELED = OUTDIR / "species_master_table_labeled.csv"
OUT_GROUP_STATS = OUTDIR / "invasion_group_summary.csv"
OUT_LOGIT_COEFS = OUTDIR / "invasion_logistic_regression_coefficients.csv"
OUT_SUMMARY_JSON = OUTDIR / "invasion_first_pass_summary.json"


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


def mann_whitney_u_stat(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    rx = ranks[:len(x)]
    u = rx.sum() - len(x) * (len(x) + 1) / 2
    return float(u)


def permutation_pvalue_mean_diff(x: np.ndarray, y: np.ndarray, n_perm: int = 999, seed: int = 42) -> float:
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan

    observed = abs(np.mean(x) - np.mean(y))
    combined = np.concatenate([x, y])
    n_x = len(x)

    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(combined)
        diff = abs(np.mean(perm[:n_x]) - np.mean(perm[n_x:]))
        count += diff >= observed

    return float((count + 1) / (n_perm + 1))


def main() -> None:
    master = pd.read_csv(MASTER_FP)
    labels = pd.read_csv(LABELS_FP)

    labeled = master.drop(columns=["invasion_status", "is_invasive"], errors="ignore").merge(
        labels[["species", "invasion_status", "is_invasive"]],
        on="species",
        how="left",
    )

    labeled["is_invasive"] = pd.to_numeric(labeled["is_invasive"], errors="coerce")

    labeled.to_csv(OUT_MASTER_LABELED, index=False)

    features = [
        "donor_score",
        "acceptor_score",
        "donor_rank",
        "acceptor_rank",
        "out_strength",
        "in_strength",
        "out_degree",
        "in_degree",
        "out_degree_centrality",
        "in_degree_centrality",
        "betweenness_centrality",
        "eigenvector_centrality",
        "mean_pairwise_asymmetry",
        "cluster",
        "cluster_k6",
        "network_community",
        "network_community_table",
        "full_geofresh_outgoing_mean",
        "climate_local_outgoing_mean",
        "shift_full_minus_climate",
        "abs_shift",
        "mean_abs_pairwise_residual",
        "donor_top20_frequency",
        "acceptor_top20_frequency",
        "n_countries",
        "n_continents",
        "range_bbox_area_km2",
        "n_basins",
        "total_records",
        "records_after_thinning",
        "temporal_span_years",
        "temporal_density",
    ]
    features = [c for c in features if c in labeled.columns]

    inv = labeled[labeled["is_invasive"] == 1].copy()
    non = labeled[labeled["is_invasive"] == 0].copy()

    rows = []
    for col in features:
        x = pd.to_numeric(inv[col], errors="coerce").to_numpy()
        y = pd.to_numeric(non[col], errors="coerce").to_numpy()

        rows.append(
            {
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
                "mann_whitney_u": mann_whitney_u_stat(x, y),
                "perm_p_mean_diff": permutation_pvalue_mean_diff(x, y, n_perm=999, seed=42),
            }
        )

    group_df = pd.DataFrame(rows).sort_values("perm_p_mean_diff", na_position="last")
    group_df.to_csv(OUT_GROUP_STATS, index=False)

    # Logistic regression first pass
    model_features = [
        c for c in [
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
        ] if c in labeled.columns
    ]

    model_df = labeled[["species", "is_invasive"] + model_features].copy()
    model_df = model_df.dropna(subset=["is_invasive"]).copy()

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

    coef_df.to_csv(OUT_LOGIT_COEFS, index=False)

    summary = {
        "n_species_total": int(len(labeled)),
        "n_invasive": int((labeled["is_invasive"] == 1).sum()),
        "n_noninvasive": int((labeled["is_invasive"] == 0).sum()),
        "n_unlabeled": int(labeled["is_invasive"].isna().sum()),
        "top_group_differences_by_perm_p": group_df.head(10)["feature"].tolist(),
        "top_logistic_regression_features": coef_df.head(10)["feature"].tolist(),
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2))

    print("Wrote:")
    print(f"  {OUT_MASTER_LABELED}")
    print(f"  {OUT_GROUP_STATS}")
    print(f"  {OUT_LOGIT_COEFS}")
    print(f"  {OUT_SUMMARY_JSON}")

    print("\nSummary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nTop group differences:")
    print(group_df.head(15).to_string(index=False))

    print("\nTop logistic regression coefficients:")
    print(coef_df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
