from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config.stage2_predictor_sets import PREDICTOR_SETS


NATIVE_VALUES = {"Native", "Type locality"}


def slugify_species(name: str) -> str:
    return name.lower().replace(" ", "_").replace("/", "_")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True, help="Exact species name")
    parser.add_argument("--predictor-set", required=True, choices=sorted(PREDICTOR_SETS.keys()))
    parser.add_argument("--input", default="data/raw/combined_data_true_master.csv")
    parser.add_argument("--output-dir", default="data/processed/stage2_pilot")
    parser.add_argument("--bg-mult", type=int, default=3)
    args = parser.parse_args()

    species = args.species
    predictor_set = args.predictor_set
    predictors = PREDICTOR_SETS[predictor_set]
    slug = slugify_species(species)

    df = pd.read_csv(args.input, low_memory=False)

    df = df[
        (df["Accuracy"].astype(str).str.lower() == "high") &
        (df["distance_m"] <= 200) &
        (df["Status"].isin(NATIVE_VALUES))
    ].copy()

    df = df.drop_duplicates(subset=["Crayfish_scientific_name", "subc_id"]).copy()

    need = ["Crayfish_scientific_name", "subc_id", "basin_id", "long_or", "lat_or"] + predictors
    df = df[need].dropna(subset=predictors).copy()

    # scale gradient fields if present
    for col in ["l_TOP15", "u_TOP15"]:
        if col in df.columns:
            df[col] = df[col] * 1.0e-6

    pres = df[df["Crayfish_scientific_name"] == species].copy()
    if pres.empty:
        raise ValueError(f"No rows found for species: {species!r} after filtering")

    pres["resp"] = 1

    focal_segments = set(pres["subc_id"])
    bg = df[(df["Crayfish_scientific_name"] != species) & (~df["subc_id"].isin(focal_segments))].copy()
    bg = bg.drop_duplicates(subset=["subc_id"]).copy()

    n_bg = min(len(bg), len(pres) * args.bg_mult)
    if n_bg == 0:
        raise ValueError("No background rows available after filtering")

    bg = bg.sample(n=n_bg, random_state=42).copy()
    bg["resp"] = 0

    out = pd.concat([pres, bg], ignore_index=True)
    out = out[["resp", "long_or", "lat_or"] + predictors + ["subc_id", "basin_id", "Crayfish_scientific_name"]]

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outfile = out_dir / f"{slug}__{predictor_set}_biomod_pilot.csv"
    out.to_csv(outfile, index=False)

    print(f"Wrote {outfile}")
    print(f"Species: {species}")
    print(f"Predictor set: {predictor_set}")
    print(f"Presence rows:   {len(pres)}")
    print(f"Background rows: {len(bg)}")
    print(f"Total rows:      {len(out)}")
    print("\nResponse counts:")
    print(out['resp'].value_counts().sort_index().to_string())
    print("\nPredictors:")
    print(predictors)


if __name__ == "__main__":
    main()
