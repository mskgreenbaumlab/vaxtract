import json
import pathlib

import figure_benchmark.model as m
import figure_benchmark.numeric_reader as nr


PANEL = m.Panel(
    panel_id="fig4b", figure_label="Figure 4B", chart_type="bar", page=1,
    bbox_frac=(0.1, 0.1, 0.4, 0.4), truth_file="x", truth_sheet="Figure 4B",
    truth_adapter="fig4b_bars", unit="% all blood T cells",
)


def test_prompt_mentions_panel_and_no_fabrication():
    prompt = nr.build_numeric_prompt(PANEL, image_path="/abs/fig4b.png")
    assert "Figure 4B" in prompt
    assert "/abs/fig4b.png" in prompt
    assert "null" in prompt.lower()           # unreadable → null
    assert "log" in prompt.lower()            # axis-scale awareness
    assert "json" in prompt.lower()


def test_parse_response_flattens_series_to_points():
    payload = {
        "series": [
            {"series_label": "blood_T_pct", "unit": "%", "points": [
                {"key": "Blood (30 weeks)", "value": 1.7},
                {"key": "Liver mass (35 weeks)", "value": 2.5},
                {"key": "Blood (41 weeks)", "value": None},
            ]},
        ]
    }
    text = "```json\n" + json.dumps(payload) + "\n```"
    pts = nr.parse_numeric_response(text)
    assert len(pts) == 3
    assert pts[0] == m.Point("blood_T_pct", "Blood (30 weeks)", 1.7, "%")
    assert pts[2].value is None


def test_parse_response_raises_on_garbage():
    import pytest
    with pytest.raises(ValueError):
        nr.parse_numeric_response("not json at all")


class FakeProvider:
    name = "fake-opus"

    def __init__(self, text):
        self._text = text
        self.calls = 0

    def extract(self, image_path, prompt, effort):
        self.calls += 1
        return self._text, {"wall_time_s": 0.01, "cost_usd": 0.0}


def test_read_panel_parses_and_caches(tmp_path):
    payload = {"series": [{"series_label": "blood_T_pct", "unit": "%",
               "points": [{"key": "Blood (30 weeks)", "value": 1.7}]}]}
    text = json.dumps(payload)
    img = tmp_path / "fig4b.png"
    img.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 64)  # fake but stable bytes
    prov = FakeProvider(text)
    cache = tmp_path / ".cache"

    r1 = nr.read_panel(PANEL, img, prov, cache_dir=cache)
    assert [p.value for p in r1.points] == [1.7]
    assert prov.calls == 1
    r2 = nr.read_panel(PANEL, img, prov, cache_dir=cache)  # cache hit
    assert prov.calls == 1  # provider NOT called again
    assert r2.from_cache is True


def test_parse_coerces_bool_and_string_values():
    # no-fabrication: a JSON bool must NOT become 1.0/0.0; numeric strings parse;
    # non-numeric strings become None (present-but-unreadable), never guessed.
    payload = {"series": [{"series_label": "s", "unit": "u", "points": [
        {"key": "a", "value": True},
        {"key": "b", "value": "1.7"},
        {"key": "c", "value": "bad"},
    ]}]}
    pts = nr.parse_numeric_response(json.dumps(payload))
    assert pts[0].value is None   # bool is not a number
    assert pts[1].value == 1.7    # numeric string parsed
    assert pts[2].value is None   # non-numeric string → None


def test_read_panel_parse_error_is_not_cached(tmp_path):
    # retry-safety invariant: a failed parse stores error + empty points, writes
    # NO cache, so a second call re-hits the provider (the failure is retryable).
    img = tmp_path / "fig.png"
    img.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 64)
    prov = FakeProvider("not json at all")
    cache = tmp_path / ".cache"

    r = nr.read_panel(PANEL, img, prov, cache_dir=cache)
    assert r.points == []
    assert r.error is not None
    assert r.from_cache is False
    r2 = nr.read_panel(PANEL, img, prov, cache_dir=cache)
    assert prov.calls == 2                 # provider called again (not cached)
    assert not list(cache.glob("*.json"))  # nothing written
