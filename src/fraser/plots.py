"""Minimal Altair plots of deterministic thresholds, not simulated Figure 3."""

import math

import altair as alt
import numpy as np

from fraser.model import (
    exponential_threshold,
    isolation_threshold,
)
from fraser.simulation import SimulationResult

DISEASE_RANGES = [
    {"disease": "SARS", "theta_min": 0, "theta_max": 0.11, "r0_min": 2, "r0_max": 4},
    {
        "disease": "Smallpox",
        "theta_min": 0,
        "theta_max": 0.2,
        "r0_min": 4,
        "r0_max": 10,
    },
    {
        "disease": "Influenza",
        "theta_min": 0.3,
        "theta_max": 0.5,
        "r0_min": 1,
        "r0_max": 21,
    },
    {"disease": "HIV", "theta_min": 0.8, "theta_max": 1, "r0_min": 1, "r0_max": 5},
]


def deterministic_chart() -> alt.VConcatChart:
    """Plot Equations 6 and 18 with historical disease-range rectangles."""
    methods = ["Isolation only", "Tracing: exponential (Eq. 18)"]
    panels = []
    for isolation in [1.0, 0.9, 0.75]:
        records = []
        for theta in np.linspace(0.01, 0.99, 99):
            records.append(
                {
                    "theta": theta,
                    "r0": isolation_threshold(theta, isolation),
                    "method": methods[0],
                }
            )
            records.append(
                {
                    "theta": theta,
                    "r0": exponential_threshold(theta, isolation, tracing=1),
                    "method": methods[1],
                }
            )
        records.extend({"theta": 1, "r0": 1, "method": method} for method in methods)
        lines = (
            alt.Chart(alt.Data(values=records))
            .mark_line(clip=True)
            .encode(  # ty: ignore[unresolved-attribute]
                x=alt.X(
                    "theta:Q",
                    title="Transmission before symptoms or asymptomatic (theta)",
                    scale=alt.Scale(domain=[0, 1]),
                    axis=alt.Axis(format="%"),
                ),
                y=alt.Y(
                    "r0:Q",
                    title="Critical R0",
                    scale=alt.Scale(domain=[0, 15], nice=False),
                    axis=alt.Axis(values=[0, 3, 6, 9, 12, 15]),
                ),
                color=alt.Color(
                    "method:N",
                    title=None,
                    scale=alt.Scale(domain=methods, range=["black", "#2864b4"]),
                    legend=alt.Legend(orient="bottom"),
                ),
                strokeDash=alt.StrokeDash(
                    "method:N",
                    title=None,
                    scale=alt.Scale(domain=methods, range=[[], [6, 3]]),
                    legend=alt.Legend(orient="bottom"),
                ),
            )
        )
        diseases = alt.Chart(alt.Data(values=DISEASE_RANGES))
        tiles = diseases.mark_rect(
            color="gray", opacity=0.12, stroke="gray", clip=True
        ).encode(  # ty: ignore[unresolved-attribute]
            x="theta_min:Q",
            x2="theta_max:Q",
            y=alt.Y("r0_min:Q", scale=alt.Scale(domain=[0, 15], nice=False)),
            y2="r0_max:Q",
        )
        labels = (
            diseases.transform_calculate(
                theta="(datum.theta_min + datum.theta_max) / 2",
                r0="(datum.r0_min + min(datum.r0_max, 15)) / 2",
            )
            .mark_text(fontSize=10, color="#555555")
            .encode(  # ty: ignore[unresolved-attribute]
                x="theta:Q",
                y=alt.Y("r0:Q", scale=alt.Scale(domain=[0, 15], nice=False)),
                text="disease:N",
            )
        )
        panels.append(
            (tiles + lines + labels).properties(
                width=480, height=210, title=f"Isolation efficacy: {isolation:.0%}"
            )
        )
    return (
        alt.vconcat(*panels)
        .resolve_scale(color="shared", strokeDash="shared")
        .properties(
            title=alt.TitleParams(
                text="Deterministic control thresholds",
                subtitle="Solid: Eq. 6; dashed: exponential Eq. 18, not simulated Figure 3",
            )
        )
    )


def incidence_records(
    result: SimulationResult,
    bin_width: float = 1.0,
) -> list[dict]:
    """Aggregate new infections over bins, excluding the seed infections."""
    if not math.isfinite(bin_width) or bin_width <= 0:
        raise ValueError("bin_width must be positive and finite")
    per_bin = round(bin_width / result.timestep)
    if per_bin < 1 or not math.isclose(per_bin * result.timestep, bin_width):
        raise ValueError("bin_width must be an integer number of timesteps")
    records = []
    for start in range(0, len(result.incidence) - 1, per_bin):
        end = min(start + per_bin, len(result.incidence) - 1)
        records.append(
            {
                "time": (result.times[start] + result.times[end]) / 2,
                "incidence": int(result.incidence[start + 1 : end + 1].sum())
                / (result.times[end] - result.times[start]),
                "r0": f"R0 = {result.r0:g}",
                "scenario": f"k={result.shape:g}, Ts={result.symptom_mean:g}, theta={result.theta:.1%}",
                "policy": f"Isolation {result.isolation:.0%}; tracing {result.tracing:.0%}",
                "truncated": not result.completed,
                "last": end == len(result.incidence) - 1,
            }
        )
    if not records:
        records.append(
            {
                "time": 0,
                "incidence": None,
                "r0": f"R0 = {result.r0:g}",
                "scenario": f"k={result.shape:g}, Ts={result.symptom_mean:g}, theta={result.theta:.1%}",
                "policy": f"Isolation {result.isolation:.0%}; tracing {result.tracing:.0%}",
                "truncated": not result.completed,
                "last": True,
            }
        )
    return records


def incidence_chart(
    results: list[SimulationResult],
    *,
    logarithmic: bool = False,
    bin_width: float = 1.0,
) -> alt.FacetChart:
    """Plot timecourses by parameter pair, marking population-limit endpoints."""
    if not results:
        raise ValueError("at least one simulation result is required")
    initial_counts = {result.initial_infections for result in results}
    initial_label = (
        f"{next(iter(initial_counts)):,} initial infections per run"
        if len(initial_counts) == 1
        else "Initial infections vary by run"
    )
    records = [
        record for result in results for record in incidence_records(result, bin_width)
    ]
    if logarithmic:
        for record in records:
            if record["incidence"] == 0:
                record["incidence"] = None
    base = alt.Chart(alt.Data(values=records)).encode(
        x=alt.X("time:Q", title="Time / Tg"),
        y=alt.Y(
            "incidence:Q",
            title="New infections / Tg",
            scale=alt.Scale(type="log" if logarithmic else "linear"),
        ),
        color=alt.Color("policy:N", title=None, legend=alt.Legend(orient="bottom")),
    )
    lines = base.mark_line()
    stops = base.transform_filter("datum.truncated && datum.last").mark_point(
        shape="cross", size=60
    )
    return (
        (lines + stops)  # ty: ignore[unsupported-operator]
        .properties(width=240, height=150)
        .facet(
            row=alt.Row("scenario:N", title=None, header=alt.Header(labelAngle=0)),
            column=alt.Column("r0:N", title=None),
        )
        .resolve_scale(y="independent")
        .properties(
            title=alt.TitleParams(
                text=f"Simulated incidence ({'log' if logarithmic else 'linear'} scale)",
                subtitle=[
                    f"{initial_label}; plotted incidence excludes seeds",
                    "Crosses mark population-limit truncation",
                ],
            )
        )
    )
