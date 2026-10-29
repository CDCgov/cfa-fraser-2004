"""Reproduce deterministic outbreak-control thresholds from Fraser et al."""

import argparse
from pathlib import Path

from fraser.model import (
    exponential_threshold,
    gamma_parameters,
    isolation_threshold,
    tracing_threshold,
)

__all__ = [
    "exponential_threshold",
    "gamma_parameters",
    "isolation_threshold",
    "tracing_threshold",
]


def main() -> None:
    """Write the deterministic figures as SVG, without notebooks or HTML."""
    from fraser.plots import deterministic_chart, exponential_chart

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, chart in [
        ("deterministic", deterministic_chart()),
        ("exponential", exponential_chart()),
    ]:
        path = args.output_dir / f"{name}.svg"
        chart.save(path)
        print(path)
