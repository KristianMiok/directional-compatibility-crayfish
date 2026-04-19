"""Task 1.3 — environmental coverage per species.

For each species (passing the candidate threshold), compute the range of
observed values for each key environmental variable, and the completeness
ratio relative to the full crayfish-occupied envelope.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def env_coverage(
    occurrences_with_env: pd.DataFrame,
    variables: list[dict],
    *,
    native_only: bool = False,
    native_status_column: str = "Status",
    native_values: list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-species environmental coverage.

    Parameters
    ----------
    occurrences_with_env
        Occurrence records joined to the environmental values for their segment.
        Must contain `species_name` plus one column per variable.
    variables
        List of variable dicts from the config (each has `name`, `column`,
        optional `categorical: true`).

    Returns
    -------
    DataFrame with columns:
        species_name, variable, vmin, vmax, range_observed, env_min, env_max,
        env_range, completeness_ratio
    For categorical variables, range_observed is the count of unique classes
    and completeness_ratio is the share of the global class set covered.

    If native_only=True, both the per-species range AND the global envelope
    are computed from native records only. This matches what Stage 2 SDM
    training will see, and is the relevant denominator for cohort decisions.
    """
    if native_only:
        if native_status_column not in occurrences_with_env.columns:
            raise ValueError(
                f"native_only=True but column '{native_status_column}' is missing. "
                f"Either add it upstream or pass native_only=False."
            )
        nv = {v.lower() for v in (native_values or ["Native", "Type locality"])}
        occurrences_with_env = occurrences_with_env[
            occurrences_with_env[native_status_column].astype(str).str.lower().isin(nv)
        ].copy()

    sp_col = "Crayfish_scientific_name"
    rows: list[dict] = []
    for var in variables:
        col = var["column"]
        if col not in occurrences_with_env.columns:
            continue

        if var.get("categorical"):
            global_classes = set(occurrences_with_env[col].dropna().unique())
            n_global = max(len(global_classes), 1)
            for sp, group in occurrences_with_env.groupby(sp_col):
                sp_classes = set(group[col].dropna().unique())
                rows.append({
                    "species_name": sp,
                    "variable": var["name"],
                    "vmin": np.nan,
                    "vmax": np.nan,
                    "range_observed": float(len(sp_classes)),
                    "env_min": np.nan,
                    "env_max": np.nan,
                    "env_range": float(n_global),
                    "completeness_ratio": len(sp_classes) / n_global,
                })
        else:
            scale = float(var.get("scale_factor", 1.0))
            if scale != 1.0:
                occurrences_with_env = occurrences_with_env.copy()
                occurrences_with_env[col] = occurrences_with_env[col] * scale
            global_min = occurrences_with_env[col].min()
            global_max = occurrences_with_env[col].max()
            global_range = float(global_max - global_min) or 1.0
            for sp, group in occurrences_with_env.groupby(sp_col):
                vmin, vmax = group[col].min(), group[col].max()
                rng = float(vmax - vmin) if pd.notna(vmin) and pd.notna(vmax) else 0.0
                rows.append({
                    "species_name": sp,
                    "variable": var["name"],
                    "vmin": vmin,
                    "vmax": vmax,
                    "range_observed": rng,
                    "env_min": global_min,
                    "env_max": global_max,
                    "env_range": global_range,
                    "completeness_ratio": rng / global_range,
                })

    return pd.DataFrame(rows)


def plot_completeness_per_species(coverage: pd.DataFrame, out_dir: str | Path) -> None:
    """One small heatmap per species: variable × completeness_ratio."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if coverage.empty or "species_name" not in coverage.columns:
        # Nothing to plot — happens when env values haven't been joined to occurrences yet.
        return

    for sp, group in coverage.groupby("species_name"):
        fig, ax = plt.subplots(figsize=(5, 0.4 * len(group) + 1))
        sorted_g = group.sort_values("variable")
        ax.barh(sorted_g["variable"], sorted_g["completeness_ratio"], color="steelblue")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Completeness vs. global crayfish-occupied envelope")
        ax.set_title(sp, fontsize=10)
        fig.tight_layout()
        # filesystem-safe name
        safe = sp.replace(" ", "_").replace("/", "_")
        fig.savefig(out_dir / f"completeness_{safe}.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
