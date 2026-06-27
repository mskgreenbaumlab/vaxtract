import json
import pathlib
import agent_core

PKT = pathlib.Path(__file__).resolve().parents[1]
ROJAS_REF = PKT / "reference_records" / "rojas_extracted.json"


def test_schema_version_is_2_14_0():
    assert agent_core.SCHEMA_VERSION == "2.16.0"


def test_validate_accepts_the_rojas_reference_record():
    ok, msg = agent_core.validate_record(ROJAS_REF.read_text())
    assert ok is True
    assert msg == "VALID"


def test_validate_rejects_a_broken_record():
    ok, msg = agent_core.validate_record("{}")
    assert ok is False
    assert msg.startswith("INVALID")


def test_save_writes_a_revalidating_file_for_a_good_record(tmp_path):
    out = tmp_path / "rojas_out.json"
    ok, msg = agent_core.save_record(ROJAS_REF.read_text(), str(out))
    assert ok is True
    assert "SAVED" in msg
    ok2, _ = agent_core.validate_record(out.read_text())
    assert ok2 is True


def test_save_refuses_a_broken_record_and_writes_nothing(tmp_path):
    out = tmp_path / "bad_out.json"
    ok, msg = agent_core.save_record("{}", str(out))
    assert ok is False
    assert "NOT SAVED" in msg
    assert not out.exists()


def test_outer_guard_quarantines_an_invalid_written_file(tmp_path):
    bad = tmp_path / "newpaper_extracted.json"
    bad.write_text("{}")
    ok, msg = agent_core.outer_guard(str(bad))
    assert ok is False
    assert not bad.exists()
    assert (tmp_path / "newpaper_extracted.json.QUARANTINED").exists()
    assert "quarantine" in msg.lower()


def test_outer_guard_passes_a_valid_written_file(tmp_path):
    good = tmp_path / "ok_extracted.json"
    good.write_text(ROJAS_REF.read_text())
    ok, msg = agent_core.outer_guard(str(good))
    assert ok is True
    assert good.exists()


def test_outer_guard_handles_a_missing_file_without_raising(tmp_path):
    missing = tmp_path / "does_not_exist.json"
    ok, msg = agent_core.outer_guard(str(missing))
    assert ok is False
    assert "not found" in msg.lower()
    # no .QUARANTINED is created for a file that never existed
    assert not (tmp_path / "does_not_exist.json.QUARANTINED").exists()


def test_save_creates_missing_parent_dirs(tmp_path):
    out = tmp_path / "nested" / "deeper" / "rojas_out.json"
    ok, msg = agent_core.save_record(ROJAS_REF.read_text(), str(out))
    assert ok is True
    assert out.exists()
    ok2, _ = agent_core.validate_record(out.read_text())
    assert ok2 is True
