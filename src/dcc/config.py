"""Config loading utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML config file into a dict."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ensure_dirs(cfg: dict[str, Any]) -> None:
    """Create output directories declared in the config if they don't exist."""
    for key in ("interim_dir", "processed_dir", "reports_dir"):
        d = cfg.get("paths", {}).get(key)
        if d:
            Path(d).mkdir(parents=True, exist_ok=True)
