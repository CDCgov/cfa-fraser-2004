import xml.etree.ElementTree as ET

import pytest

import fraser
from fraser import main
from fraser.plots import deterministic_chart, incidence_chart, incidence_records
from fraser.simulation import simulate


def test_svg_export(tmp_path):
    path = tmp_path / "deterministic.svg"
    deterministic_chart().save(path)
    root = ET.parse(path).getroot()
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    text = " ".join(root.itertext())
    assert "Deterministic control thresholds" in text
    assert "exponential Eq. 18" in text
    assert "Eq. 16" not in text
    assert "100%" in text


def test_deterministic_bounds_and_axis(tmp_path):
    chart = deterministic_chart()
    spec = chart.to_dict()
    for panel in spec["vconcat"]:
        for layer in panel["layer"]:
            assert layer["encoding"]["y"]["scale"] == {
                "domain": [0, 15],
                "nice": False,
            }
        assert panel["layer"][1]["encoding"]["strokeDash"]["scale"]["range"] == [
            [],
            [6, 3],
        ]
        assert panel["layer"][1]["encoding"]["y"]["axis"]["values"] == [
            0,
            3,
            6,
            9,
            12,
            15,
        ]
    path = tmp_path / "deterministic.svg"
    chart.save(path)
    root = ET.parse(path).getroot()
    for color, pattern in [("black", ""), ("#2864b4", "6,3")]:
        curves = [
            element
            for element in root.iter()
            if element.get("stroke") == color and element.tag.endswith("path")
        ]
        assert len(curves) == 4
        assert all(curve.get("stroke-dasharray", "") == pattern for curve in curves)


@pytest.mark.parametrize("logarithmic", [False, True])
def test_incidence_svg(tmp_path, logarithmic):
    results = [
        simulate(1.5, 1, 2, initial_infections=20, duration=2, tracing=tracing)
        for tracing in [0, 1]
    ]
    path = tmp_path / "incidence.svg"
    incidence_chart(results, logarithmic=logarithmic).save(path)
    root = ET.parse(path).getroot()
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    text = " ".join(root.itertext())
    assert "Simulated incidence" in text
    assert "theta=50.0%" in text
    assert "20 initial infections per run" in text


def test_incidence_bins_preserve_counts_and_partial_bin():
    result = simulate(1.5, 1, 2, initial_infections=100, duration=2.5)
    records = incidence_records(result)
    assert len(records) == 3
    assert (
        sum(record["incidence"] for record in records[:2])
        + records[2]["incidence"] * 0.5
        == result.incidence.sum()
    )


def test_log_plot_keeps_zero_incidence_as_missing():
    chart = incidence_chart([simulate(0, 1, 2, duration=2)], logarithmic=True)
    spec = chart.to_dict()
    assert all(record["incidence"] is None for record in spec["data"]["values"])


def test_population_limit_is_marked():
    result = simulate(
        20, 1, 2, initial_infections=100, max_infections=200, timestep=0.1, duration=2
    )
    records = incidence_records(result)
    assert records[-1]["truncated"]
    assert records[-1]["last"]


@pytest.mark.parametrize("output_after_command", [False, True])
def test_simulation_command(tmp_path, monkeypatch, output_after_command):
    monkeypatch.setattr(
        fraser,
        "SIMULATION_SETTINGS",
        {
            **fraser.SIMULATION_SETTINGS,
            "r0": [2],
            "symptom_mean": [1],
            "shape": [2],
            "initial_infections": 20,
            "timestep": 0.1,
            "duration": 1,
        },
    )
    output_args = ["--output-dir", str(tmp_path)]
    monkeypatch.setattr(
        "sys.argv",
        ["fraser", "simulate", *output_args]
        if output_after_command
        else ["fraser", *output_args, "simulate"],
    )
    main()
    assert (tmp_path / "simulation-linear.svg").is_file()
    assert (tmp_path / "simulation-log.svg").is_file()
