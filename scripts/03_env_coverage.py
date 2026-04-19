"""Stage 1, Task 1.3 — environmental coverage per species.

Note: this expects occurrences to already carry environmental values per
segment. If they don't, do that join in src/dcc/data.py first (TODO once
the GeoFRESH attribute table is settled).
"""

from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from dcc.config import ensure_dirs, load_config
from dcc.data import deduplicate_by_segment, filter_records, load_woc
from dcc.env_coverage import env_coverage, plot_completeness_per_species
from dcc.logging_utils import setup_logger


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.3", cfg["logging"]["level"], cfg["logging"]["file"])

    woc = load_woc(cfg["paths"]["woc_occurrences"])
    woc = filter_records(
        woc,
        high_accuracy_only=cfg["filtering"]["high_accuracy_only"],
        max_snap_distance_m=cfg["filtering"]["max_snap_distance_m"],
    )
    woc = deduplicate_by_segment(woc)

    # Restrict to species at or above the candidate threshold.
    inv = pd.read_csv(Path(cfg["paths"]["processed_dir"]) / "species_inventory.csv")
    keep = set(inv.loc[inv["records_deduplicated_segment"] >= cfg["thresholds"]["default"],
                       "species_name"])
    log.info("Coverage restricted to %d species at threshold %d",
             len(keep), cfg["thresholds"]["default"])
    woc = woc[woc["species_name"].isin(keep)]

    ec_cfg = cfg.get("env_coverage", {})
    coverage = env_coverage(
        woc,
        cfg["env_variables"],
        native_only=ec_cfg.get("native_only", False),
        native_status_column=ec_cfg.get("native_status_column", "native_status"),
        native_value=ec_cfg.get("native_value", "native"),
    )
    if ec_cfg.get("native_only"):
        log.info("env_coverage computed on native-only records")
    out_csv = Path(cfg["paths"]["processed_dir"]) / "env_coverage.csv"
    coverage.to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(coverage))

    plot_completeness_per_species(coverage, Path(cfg["paths"]["reports_dir"]) / "completeness")
    log.info("Per-species completeness plots written to %s/completeness/",
             cfg["paths"]["reports_dir"])


if __name__ == "__main__":
    main()
