"""Simulation-led control thresholds with inexpensive analytical references."""

import math
from typing import cast

import altair as alt
import numpy as np

from fraser.model import exponential_threshold, isolation_threshold
from fraser.simulation import CriticalResult

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


def critical_chart(results: list[CriticalResult]) -> alt.FacetChart:
    """Show simulated thresholds, analytical references, and historical ranges."""
    shapes = sorted({result.shape for result in results}, reverse=True)
    series = (
        ["Isolation only (Eq. 6)"]
        + [f"Tracing (k={shape:g})" for shape in shapes]
        + ["Exponential tracing (Eq. 18)"]
    )
    colors = (
        ["black"]
        + ["#2864b4" if shape == 10 else "#d34836" for shape in shapes]
        + ["#777777"]
    )
    records = []
    for result in results:
        value = result.estimate
        if result.status == "above_range":
            value = result.upper
        elif result.status == "below_range":
            value = result.lower
        records.append(
            {
                "theta": result.theta,
                "r0": value if math.isfinite(value) else None,
                "lower": result.lower,
                "upper": result.upper,
                "efficacy": f"Isolation {result.isolation:.0%}",
                "series": f"Tracing (k={result.shape:g})",
                "status": result.status,
                "source": "Simulation",
            }
        )
    for isolation, tracing in sorted(
        {(result.isolation, result.tracing) for result in results}
    ):
        for theta in np.linspace(0.01, 1, 100):
            for name, value in [
                ("Isolation only (Eq. 6)", isolation_threshold(theta, isolation)),
                (
                    "Exponential tracing (Eq. 18)",
                    exponential_threshold(theta, isolation, tracing),
                ),
            ]:
                records.append(
                    {
                        "theta": theta,
                        "r0": value,
                        "efficacy": f"Isolation {isolation:.0%}",
                        "series": name,
                        "status": "reference",
                        "source": "Deterministic",
                    }
                )
    for isolation in sorted({result.isolation for result in results}):
        for disease in DISEASE_RANGES:
            records.append(
                {
                    **disease,
                    "efficacy": f"Isolation {isolation:.0%}",
                    "source": "Disease",
                }
            )
    y_scale = alt.Scale(domain=[0, 15], nice=False)
    base = alt.Chart(alt.Data(values=records)).encode(
        x=alt.X(
            "theta:Q",
            title="Transmission before symptoms (theta)",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(format="%"),
        ),
        y=alt.Y(
            "r0:Q",
            title="Critical R0",
            scale=y_scale,
            axis=alt.Axis(values=[0, 3, 6, 9, 12, 15]),
        ),
    )
    curves = base.encode(
        color=alt.Color(
            "series:N",
            title=None,
            scale=alt.Scale(domain=series, range=colors),
            legend=alt.Legend(orient="bottom", columns=2),
        ),
    )
    lines = cast(
        alt.Chart,
        curves.transform_filter(
            "datum.source == 'Deterministic' || datum.status == 'bracketed'"
        )
        .encode(
            strokeDash=alt.StrokeDash(
                "series:N",
                title=None,
                scale=alt.Scale(
                    domain=series,
                    range=[[] for _ in series[:-1]] + [[6, 3]],
                ),
                legend=None,
            ),
            strokeWidth=alt.condition(
                "datum.source == 'Simulation'", alt.value(2.5), alt.value(1)
            ),
            opacity=alt.condition(
                "datum.series == 'Exponential tracing (Eq. 18)'",
                alt.value(0.5),
                alt.value(1),
            ),
        )
        .mark_line(clip=True),
    )
    bars = cast(
        alt.Chart,
        curves.transform_filter("datum.status == 'bracketed'")
        .encode(
            y=alt.Y("lower:Q", scale=y_scale),
            y2="upper:Q",
        )
        .mark_rule(strokeWidth=2),
    )
    points = cast(
        alt.Chart,
        curves.transform_filter("datum.source == 'Simulation'")
        .encode(
            shape=alt.Shape(
                "status:N",
                title=None,
                legend=None,
                scale=alt.Scale(
                    domain=["bracketed", "above_range", "below_range"],
                    range=["circle", "triangle-up", "triangle-down"],
                ),
            ),
        )
        .mark_point(filled=True, size=70),
    )
    diseases = base.transform_filter("datum.source == 'Disease'")
    tiles = cast(
        alt.Chart,
        diseases.encode(
            x=alt.X("theta_min:Q", scale=alt.Scale(domain=[0, 1])),
            x2="theta_max:Q",
            y=alt.Y("r0_min:Q", scale=y_scale),
            y2="r0_max:Q",
        ).mark_rect(color="gray", opacity=0.07, clip=True),
    )
    labels = cast(
        alt.Chart,
        diseases.transform_calculate(
            theta="(datum.theta_min + datum.theta_max) / 2",
            r0="(datum.r0_min + min(datum.r0_max, 15)) / 2",
        )
        .encode(text="disease:N")
        .mark_text(
            color="gray",
            fontSize=9,
        ),
    )
    unresolved = sum(result.status == "unresolved" for result in results)
    return (
        alt.layer(tiles, labels, lines, bars, points)
        .properties(width=480, height=220)
        .facet(
            row=alt.Row(
                "efficacy:N",
                title=None,
                sort=[
                    f"Isolation {value:.0%}"
                    for value in sorted(
                        {result.isolation for result in results}, reverse=True
                    )
                ],
                header=alt.Header(labelAngle=0),
            ),
        )
        .properties(
            title=alt.TitleParams(
                text="Simulated outbreak-control thresholds",
                subtitle=[
                    "Bars: search brackets, not confidence intervals; triangles: outside search range",
                    f"Black: Eq. 6; blue/red: simulated gamma shapes; dashed gray: Eq. 18; unresolved: {unresolved}",
                ],
            )
        )
    )
