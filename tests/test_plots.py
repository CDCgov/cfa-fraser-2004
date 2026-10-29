import xml.etree.ElementTree as ET

from fraser.plots import deterministic_chart, exponential_chart


def test_svg_export(tmp_path):
    path = tmp_path / "exponential.svg"
    exponential_chart().save(path)
    root = ET.parse(path).getroot()
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    text = " ".join(root.itertext())
    assert "Exponential special case" in text
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
    for color, pattern in [("#2864b4", "6,3"), ("#d34836", "6,3")]:
        curves = [
            element
            for element in root.iter()
            if element.get("stroke") == color and element.tag.endswith("path")
        ]
        assert len(curves) == 4
        assert all(curve.get("stroke-dasharray") == pattern for curve in curves)
