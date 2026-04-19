"""Stage 1, Task 1.2 — threshold sensitivity table + plot."""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from dcc.config import ensure_dirs, load_config
from dcc.data import load_woc
from dcc.logging_utils import setup_logger
from dcc.threshold_sensitivity import plot_threshold_sweep, threshold_sweep


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.2", cfg["logging"]["level"], cfg["logging"]["file"])

    inv_path = Path(cfg["paths"]["processed_dir"]) / "species_inventory.csv"
    inventory = pd.read_csv(inv_path)
    woc = load_woc(cfg["paths"]["woc_occurrences"])

    sweep = threshold_sweep(inventory, cfg["thresholds"]["values"], woc)
    out_csv = Path(cfg["paths"]["processed_dir"]) / "threshold_sensitivity.csv"
    sweep.to_csv(out_csv, index=False)
    log.info("Wrote %s", out_csv)

    out_png = Path(cfg["paths"]["reports_dir"]) / "threshold_sensitivity_plot.png"
    plot_threshold_sweep(sweep, out_png, default=cfg["thresholds"]["default"])
    log.info("Wrote %s", out_png)


if __name__ == "__main__":
    main()
