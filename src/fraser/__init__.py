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
from fraser.simulation import SimulationResult, presymptomatic_fraction, simulate

__all__ = [
    "exponential_threshold",
    "gamma_parameters",
    "isolation_threshold",
    "SimulationResult",
    "presymptomatic_fraction",
    "simulate",
]

SIMULATION_SETTINGS: dict = {
    "r0": [2, 3],
    "symptom_mean": [0.5, 1],
    "shape": [2, 10],
    "isolation": 0.9,
    "tracing": 1,
    "initial_infections": 1_000,
    "timestep": 0.01,
    "duration": 10,
    "max_infections": 1_000_000,
    "seed": 2004,
    "bin_width": 1,
}


def main() -> None:
    """Write deterministic or exploratory simulation figures as SVG."""
    from fraser.plots import deterministic_chart

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures"))
    commands = parser.add_subparsers(dest="command")
    simulation = commands.add_parser(
        "simulate", help="plot timestep simulation timecourses"
    )
    simulation.add_argument("--output-dir", type=Path, default=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.command == "simulate":
        _simulation_figures(args.output_dir)
        return
    path = args.output_dir / "deterministic.svg"
    deterministic_chart().save(path)
    print(path)


def _simulation_figures(output_dir: Path) -> None:
    from fraser.plots import incidence_chart

    parameters = SIMULATION_SETTINGS
    tracing_values = [0] if parameters["tracing"] == 0 else [0, parameters["tracing"]]
    settings = list(
        product(
            parameters["shape"],
            parameters["symptom_mean"],
            parameters["r0"],
            tracing_values,
        )
    )
    streams = np.random.SeedSequence(parameters["seed"]).spawn(len(settings))
    results = []
    for (shape, symptom_mean, r0, tracing), stream in zip(settings, streams):
        result = simulate(
            r0,
            symptom_mean,
            shape,
            isolation=parameters["isolation"],
            tracing=tracing,
            initial_infections=parameters["initial_infections"],
            timestep=parameters["timestep"],
            duration=parameters["duration"],
            max_infections=parameters["max_infections"],
            seed=int(stream.generate_state(1)[0]),
        )
        results.append(result)
        status = "complete" if result.completed else "POPULATION LIMIT"
        print(
            f"k={shape:g} Ts={symptom_mean:g} R0={r0:g} tracing={tracing:g}: "
            f"{result.total_infections} infections, t={result.times[-1]:g}, {status}"
        )
    for label, logarithmic in [("linear", False), ("log", True)]:
        path = output_dir / f"simulation-{label}.svg"
        incidence_chart(
            results, logarithmic=logarithmic, bin_width=parameters["bin_width"]
        ).save(path)
        print(path)
