"""Stage 1, Task 1.1 — write species_inventory.csv."""

from __future__ import annotations

from pathlib import Path

import click

from dcc.config import ensure_dirs, load_config
from dcc.data import load_woc
from dcc.inventory import build_inventory
from dcc.logging_utils import setup_logger


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.1", cfg["logging"]["level"], cfg["logging"]["file"])

    woc_path = cfg["paths"]["woc_occurrences"]
    log.info("Loading WoC occurrences from %s", woc_path)
    woc = load_woc(woc_path)
    log.info("Loaded %d raw records across %d species", len(woc), woc["Crayfish_scientific_name"].nunique())

    inv = build_inventory(woc, max_snap_m=cfg["filtering"]["max_snap_distance_m"])
    out = Path(cfg["paths"]["processed_dir"]) / "species_inventory.csv"
    inv.to_csv(out, index=False)
    log.info("Wrote %s (%d species)", out, len(inv))


if __name__ == "__main__":
    main()
