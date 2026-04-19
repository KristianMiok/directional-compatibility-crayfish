from __future__ import annotations

from pathlib import Path
import pandas as pd

pilot_dir = Path("data/processed/stage2_pilot")
files = sorted(pilot_dir.glob("*_esm_*_evaluations.csv"))

rows = []
for f in files:
    df = pd.read_csv(f)
    stem = f.stem  # e.g. astacus_astacus_esm_glm_evaluations
    if "_esm_" not in stem:
        continue
    species_slug, rest = stem.split("_esm_", 1)
    model_core = rest.replace("_evaluations", "").upper()

    for _, row in df.iterrows():
        rows.append({
            "species_slug": species_slug,
            "model_core": model_core,
            "RUN": row.get("RUN"),
            "AUC": row.get("AUC"),
            "TSS": row.get("TSS"),
            "Boyce": row.get("Boyce"),
            "Kappa": row.get("Kappa"),
            "MPA": row.get("MPA"),
        })

out = pd.DataFrame(rows)

outfile = pilot_dir / "pilot_species_core_comparison.csv"
out.to_csv(outfile, index=False)

print(f"Wrote {outfile}")

if out.empty:
    print("\nNo evaluation files found.")
else:
    print("\nMean metrics by species x core:")
    print(
        out.groupby(["species_slug", "model_core"])[["AUC", "TSS", "Boyce", "Kappa", "MPA"]]
        .mean()
        .round(4)
        .to_string()
    )
