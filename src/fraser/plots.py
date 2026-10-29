"""Minimal Altair plots of deterministic thresholds, not simulated Figure 3."""

import math

import altair as alt
import numpy as np

from fraser.model import (
    exponential_threshold,
    isolation_threshold,
)

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
        at_zero = isolation_threshold(0, isolation)
        if math.isfinite(at_zero):
            records.append({"theta": 0, "r0": at_zero, "method": methods[0]})
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
