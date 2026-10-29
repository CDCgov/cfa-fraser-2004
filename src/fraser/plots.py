"""Minimal Altair plots of deterministic thresholds, not simulated Figure 3."""

import math

import altair as alt
import numpy as np

from fraser.model import (
    exponential_threshold,
    isolation_threshold,
    tracing_threshold,
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
    """Plot nine threshold curves and historical disease-range rectangles."""
    methods = ["Isolation only", "Tracing: low variance", "Tracing: high variance"]
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
            for variance_ratio, method in zip([0.1, 0.5], methods[1:]):
                records.append(
                    {
                        "theta": theta,
                        "r0": tracing_threshold(
                            theta, isolation, variance_ratio=variance_ratio
                        ),
                        "method": method,
                    }
                )
        records.extend({"theta": 1, "r0": 1, "method": method} for method in methods)
        at_zero = isolation_threshold(0, isolation)
        if math.isfinite(at_zero):
            records.append({"theta": 0, "r0": at_zero, "method": methods[0]})
        lines = (
            alt.Chart(alt.Data(values=records))
            .mark_line(clip=True)
            .encode(
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
                    scale=alt.Scale(
                        domain=methods, range=["black", "#2864b4", "#d34836"]
                    ),
                    legend=alt.Legend(orient="bottom"),
                ),
                strokeDash=alt.StrokeDash(
                    "method:N",
                    title=None,
                    scale=alt.Scale(domain=methods, range=[[], [6, 3], [6, 3]]),
                    legend=alt.Legend(orient="bottom"),
                ),
            )
        )
        diseases = alt.Chart(alt.Data(values=DISEASE_RANGES))
        tiles = diseases.mark_rect(
            color="gray", opacity=0.12, stroke="gray", clip=True
        ).encode(
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
            .encode(
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
                subtitle="Dashed: deterministic upper bounds (Eq. 16), not simulated Figure 3",
            )
        )
    )


def exponential_chart() -> alt.Chart:
    """Compare Equation 16 quadrature to corrected Equation 18 for exponentials."""
    records = []
    for isolation in [1.0, 0.9, 0.75]:
        for theta in np.linspace(0.05, 0.95, 19):
            for method, value in [
                ("Corrected Eq. 18", exponential_threshold(theta, isolation)),
                (
                    "Numerical Eq. 16",
                    tracing_threshold(theta, isolation, variance_ratio=1),
                ),
            ]:
                records.append(
                    {
                        "theta": theta,
                        "r0": value,
                        "isolation": f"{isolation:.0%}",
                        "method": method,
                    }
                )
    base = alt.Chart(alt.Data(values=records)).encode(
        x=alt.X(
            "theta:Q",
            title="theta",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(format="%"),
        ),
        y=alt.Y(
            "r0:Q",
            title="Critical R0",
            scale=alt.Scale(domain=[0, 15], nice=False),
            axis=alt.Axis(values=[0, 3, 6, 9, 12, 15]),
        ),
        color=alt.Color("isolation:N", title="Isolation efficacy"),
    )
    lines = base.transform_filter(alt.datum.method == "Corrected Eq. 18").mark_line(
        clip=True
    )
    points = base.transform_filter(alt.datum.method == "Numerical Eq. 16").mark_point(
        clip=True
    )
    return (lines + points).properties(
        width=480,
        height=270,
        title=alt.TitleParams(
            text="Exponential special case",
            subtitle="Lines: corrected Eq. 18; points: numerical Eq. 16; tracing efficacy 100%",
        ),
    )
