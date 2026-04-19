"""Stage 1, Task 1.4 — basin overlap matrix."""

from __future__ import annotations

from pathlib import Path

import click

from dcc.basin_overlap import basin_overlap_matrices
from dcc.config import ensure_dirs, load_config
from dcc.data import load_woc
from dcc.logging_utils import setup_logger


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.4", cfg["logging"]["level"], cfg["logging"]["file"])

    woc = load_woc(cfg["paths"]["woc_occurrences"])
    shared, jacc = basin_overlap_matrices(
        woc, basin_id_column=cfg["basin_overlap"]["basin_id_column"]
    )

    out_dir = Path(cfg["paths"]["processed_dir"])
    shared.to_csv(out_dir / "basin_overlap_matrix.csv")
    jacc.to_csv(out_dir / "basin_jaccard_matrix.csv")
    log.info("Wrote %s and %s (%d×%d)",
             out_dir / "basin_overlap_matrix.csv",
             out_dir / "basin_jaccard_matrix.csv",
             *shared.shape)


if __name__ == "__main__":
    main()
