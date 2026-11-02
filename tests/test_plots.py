import xml.etree.ElementTree as ET

import pytest

import fraser
from fraser import main
from fraser.plots import (
    critical_chart,
)
from fraser.simulation import CriticalResult, SearchParameters


def search_parameters(theta, isolation=0.9, tracing=1, shape=2):
    _, _, scale = fraser.gamma_parameters(
        theta, variance_ratio=1 / shape, generation_time=1
    )
    return SearchParameters(
        symptom_mean=shape * scale,
        shape=shape,
        isolation=isolation,
        tracing=tracing,
    )


def test_critical_svg_and_out_of_range_markers(tmp_path):
    results = [
        CriticalResult(
            parameters=search_parameters(0.5),
            lower=3.5,
            upper=3.7,
            status="bracketed",
        ),
        CriticalResult(
            parameters=search_parameters(0.2),
            lower=1,
            upper=15,
            status="above_range",
        ),
        CriticalResult(
            parameters=search_parameters(0.8),
            lower=1,
            upper=15,
            status="unresolved",
        ),
    ]
    chart = critical_chart(results)
    path = tmp_path / "critical.svg"
    chart.save(path)
    text = " ".join(ET.parse(path).getroot().itertext())
    assert "Simulated outbreak-control thresholds" in text
    assert "unresolved: 1" in text
    assert "Smallpox" in text
    records = chart.to_dict()["data"]["values"]
    assert records[0]["r0"] == 3.6
    assert records[1]["r0"] == 15
    assert records[2]["r0"] is None
    references = [record for record in records if record["source"] == "Deterministic"]
    assert len(references) == 200
    assert references[-1]["r0"] == 1
    lines = chart.to_dict()["spec"]["layer"][2]
    assert lines["encoding"]["strokeDash"]["scale"]["range"] == [[], [], [6, 3]]
    assert lines["encoding"]["strokeWidth"]["condition"]["value"] == 2.5


def test_combined_isolation_reference():
    result = CriticalResult(
        parameters=search_parameters(0.5),
        lower=1.7,
        upper=1.9,
        status="bracketed",
    )
    records = critical_chart([result]).to_dict()["data"]["values"]
    reference = next(
        record
        for record in records
        if record.get("series") == "Isolation only (Eq. 6)" and record["theta"] == 0.5
    )
    assert reference["r0"] == fraser.isolation_threshold(0.5, 0.9)


def test_efficacy_panels_have_both_gamma_shapes_and_both_references():
    results = [
        CriticalResult(
            parameters=search_parameters(0.5, isolation=isolation, shape=shape),
            lower=3,
            upper=3.5,
            status="bracketed",
        )
        for isolation in [1, 0.9, 0.75]
        for shape in [2, 10]
    ]
    spec = critical_chart(results).to_dict()
    assert spec["facet"]["row"]["field"] == "efficacy"
    color_scale = spec["spec"]["layer"][2]["encoding"]["color"]["scale"]
    palette = dict(zip(color_scale["domain"], color_scale["range"]))
    assert palette["Tracing (k=10)"] == "#2864b4"
    assert palette["Tracing (k=2)"] == "#d34836"
    assert spec["facet"]["row"]["sort"] == [
        "Isolation 100%",
        "Isolation 90%",
        "Isolation 75%",
    ]
    for efficacy in spec["facet"]["row"]["sort"]:
        curves = {
            record["series"]
            for record in spec["data"]["values"]
            if record["efficacy"] == efficacy and "series" in record
        }
        assert curves == {
            "Isolation only (Eq. 6)",
            "Tracing (k=2)",
            "Tracing (k=10)",
            "Exponential tracing (Eq. 18)",
        }


@pytest.mark.parametrize("theta", [0.05, 0.5, 0.95])
@pytest.mark.parametrize("shape", [2, 10])
def test_critical_command(tmp_path, monkeypatch, theta, shape):
    monkeypatch.setattr(
        fraser,
        "SIMULATION_SETTINGS",
        {
            **fraser.SIMULATION_SETTINGS,
            "shape": [shape],
            "theta": [theta],
            "isolation": [0.9],
        },
    )
    calls = []

    def fake_search(base_parameters, **kwargs):
        calls.append(base_parameters)
        return CriticalResult(
            parameters=base_parameters,
            lower=3.5,
            upper=3.7,
            status="bracketed",
        )

    monkeypatch.setattr(
        fraser,
        "find_critical_r0",
        fake_search,
    )
    monkeypatch.setattr("sys.argv", ["fraser", "--output-dir", str(tmp_path)])
    main()
    assert [path.name for path in tmp_path.iterdir()] == ["control.svg"]
    assert len(calls) == 1
    base_parameters = calls[0]
    assert (
        base_parameters.shape,
        base_parameters.isolation,
        base_parameters.tracing,
    ) == (
        shape,
        0.9,
        1,
    )
    assert base_parameters.theta == pytest.approx(theta)
