"""Task 1.4 — basin overlap preliminary matrix.

Sanity check: which species share basins, which are geographically
independent. Produces a species × species matrix of shared-basin counts and
a parallel matrix of Jaccard similarities.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def basin_overlap_matrices(
    woc: pd.DataFrame,
    basin_id_column: str = "basin_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (shared_basins, jaccard) matrices, both species × species.

    Self-overlap on the diagonal equals the species' own basin count
    (shared) and 1.0 (Jaccard).
    """
    sp_basins: dict[str, set] = (
        woc.groupby("Crayfish_scientific_name")[basin_id_column]
        .apply(lambda s: set(s.dropna())).to_dict()
    )
    species = sorted(sp_basins)
    n = len(species)
    shared = np.zeros((n, n), dtype=int)
    jacc = np.zeros((n, n), dtype=float)

    for i, a in enumerate(species):
        sa = sp_basins[a]
        for j in range(i, n):
            b = species[j]
            sb = sp_basins[b]
            inter = len(sa & sb)
            union = len(sa | sb) or 1
            shared[i, j] = shared[j, i] = inter
            jacc[i, j] = jacc[j, i] = inter / union

    shared_df = pd.DataFrame(shared, index=species, columns=species)
    jacc_df = pd.DataFrame(jacc, index=species, columns=species)
    return shared_df, jacc_df
