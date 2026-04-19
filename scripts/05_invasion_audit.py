"""Stage 1, Task 1.5 — invasion contamination audit.

Output goes to manual review by Lucian / Mihaela / Dave. Nothing here
auto-removes records — that happens after their validation.
"""

from __future__ import annotations

from pathlib import Path

import click

from dcc.config import ensure_dirs, load_config
from dcc.data import load_woc
from dcc.invasion_audit import flag_records, load_native_ranges
from dcc.logging_utils import setup_logger


@click.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def main(config_path: str) -> None:
    cfg = load_config(config_path)
    ensure_dirs(cfg)
    log = setup_logger("task1.5", cfg["logging"]["level"], cfg["logging"]["file"])

    woc = load_woc(cfg["paths"]["woc_occurrences"])

    nat_path = cfg["invasion_audit"]["native_range_source"]
    if not Path(nat_path).exists():
        log.warning(
            "Native ranges file not found at %s — every record will be flagged "
            "as 'unknown_native_range'. Populate this CSV before P1.",
            nat_path,
        )
        native_ranges: dict = {}
    else:
        native_ranges = load_native_ranges(nat_path)
        log.info("Loaded native-range info for %d species", len(native_ranges))

    flagged = flag_records(
        woc,
        native_ranges,
        flag_outside_continent=cfg["invasion_audit"]["flag_outside_continent"],
        post_year=cfg["invasion_audit"]["flag_post_year"],
    )

    out = Path(cfg["paths"]["processed_dir"]) / "invasion_contamination_flagged.csv"
    flagged.to_csv(out, index=False)
    log.info(
        "Wrote %s — %d records flagged across %d species (for manual review)",
        out,
        len(flagged),
        flagged["species_name"].nunique() if not flagged.empty else 0,
    )


if __name__ == "__main__":
    main()
