"""Stage 1, Task 1.5 — invasion contamination audit."""

from __future__ import annotations

from pathlib import Path

import click

from dcc.config import ensure_dirs, load_config
from dcc.data import load_woc
from dcc.invasion_audit import flag_records
from dcc.logging_utils import setup_logger


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.5", cfg["logging"]["level"], cfg["logging"]["file"])

    woc = load_woc(cfg["paths"]["woc_occurrences"])
    nv = cfg.get("invasion_audit", {}).get("native_values", ["Native", "Type locality"])
    flagged = flag_records(woc, native_values=nv)

    out = Path(cfg["paths"]["processed_dir"]) / "invasion_contamination_flagged.csv"
    flagged.to_csv(out, index=False)
    log.info(
        "Wrote %s — %d records flagged across %d species (for manual review)",
        out, len(flagged),
        flagged["Crayfish_scientific_name"].nunique() if not flagged.empty else 0,
    )


if __name__ == "__main__":
    main()
