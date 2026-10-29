"""Reproduce deterministic outbreak-control thresholds from Fraser et al."""

import argparse
from itertools import product
from pathlib import Path

import numpy as np

from fraser.model import (
    exponential_threshold,
    gamma_parameters,
    isolation_threshold,
)
from fraser.simulation import (
    CriticalResult,
    SimulationResult,
    binned_incidence,
    estimate_growth,
    find_critical_r0,
    presymptomatic_fraction,
    simulate,
)

__all__ = [
    "exponential_threshold",
    "gamma_parameters",
    "isolation_threshold",
    "SimulationResult",
    "presymptomatic_fraction",
    "simulate",
    "binned_incidence",
    "estimate_growth",
    "CriticalResult",
    "find_critical_r0",
]

SIMULATION_SETTINGS: dict = {
    "theta": [0.05, 0.2, 0.4, 0.6, 0.8, 0.95],
    "shape": [2, 10],
    "isolation": [1, 0.9, 0.75],
    "tracing": 1,
    "initial_infections": 300,
    "timestep": 0.02,
    "duration": 10,
    "max_infections": 1_000_000,
    "seed": 2004,
    "bin_width": 1,
}

SEARCH_SETTINGS: dict = {
    "lower": 1,
    "upper": 15,
    "tolerance": 0.5,
    "repeats": 3,
}


def main() -> None:
    """Search simulated thresholds and write one SVG with analytical overlays."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _critical_figures(args.output_dir)


def _critical_figures(output_dir: Path) -> None:
    from fraser.plots import critical_chart

    parameters = SIMULATION_SETTINGS
    settings = list(
        product(parameters["isolation"], parameters["shape"], parameters["theta"])
    )
    streams = np.random.SeedSequence(parameters["seed"]).spawn(len(settings))
    results = []
    for (isolation, shape, theta), stream in zip(settings, streams):
        _, _, incubation_scale = gamma_parameters(
            theta, variance_ratio=1 / shape, generation_time=1
        )
        symptom_mean = shape * incubation_scale
        result = find_critical_r0(
            symptom_mean,
            shape,
            isolation=isolation,
            tracing=parameters["tracing"],
            initial_infections=parameters["initial_infections"],
            timestep=parameters["timestep"],
            duration=parameters["duration"],
            max_infections=parameters["max_infections"],
            bin_width=parameters["bin_width"],
            seed=int(stream.generate_state(1)[0]),
            **SEARCH_SETTINGS,
        )
        results.append(result)
        detail = (
            f"R0={result.estimate:.2f} [{result.lower:.2f}, {result.upper:.2f}]"
            if result.status == "bracketed"
            else result.status
        )
        print(
            f"isolation={isolation:.0%} k={shape:g} theta={theta:g} Ts={symptom_mean:.3g}: {detail}",
            flush=True,
        )
    path = output_dir / "control.svg"
    critical_chart(results).save(path)
    print(path)
