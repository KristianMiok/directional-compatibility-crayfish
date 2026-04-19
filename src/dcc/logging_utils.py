"""Centralized logging setup."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def setup_logger(name: str = "dcc", level: str = "INFO", file: str | None = None) -> logging.Logger:
    """Return a configured logger that writes to stdout and optionally a file."""
    logger = logging.getLogger(name)
    if logger.handlers:  # already configured
        return logger

    logger.setLevel(level)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(name)s | %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if file:
        Path(file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(file, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger
