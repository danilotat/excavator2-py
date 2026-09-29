"""Keep failed or malformed performance evidence out of the README."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "benchmark_readme", ROOT / "tools/update_benchmark_readme.py"
)
renderer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)


@pytest.fixture
def report():
    return {
        "passed": True,
        "source_commit": "a" * 40,
        "run_url": "https://github.com/danilotat/excavator2-py/actions/runs/123",
        "windows": 4113,
        "stages": {
            "target": {"original_seconds": 10, "port_seconds": 2, "runs": 1},
            "paired_analysis": {"original_seconds": 3, "port_seconds": 2, "runs": 3},
        },
    }


def test_update_preserves_document_and_is_idempotent(report):
    before = "Introduction\n" + renderer.START + "\nstale\n" + renderer.END + "\nInstallation\n"
    after = renderer.update(before, report)
    assert after.startswith("Introduction\n" + renderer.START)
    assert after.endswith(renderer.END + "\nInstallation\n")
    assert "5.00×" in after and "1.50×" in after and "4,113" in after
    assert "stale" not in after
    assert renderer.update(after, report) == after


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf"), True, "1"])
def test_rejects_invalid_timing(report, value):
    report["stages"]["target"]["port_seconds"] = value
    with pytest.raises(ValueError, match="positive finite"):
        renderer.render(report)


def test_rejects_failed_parity_and_insufficient_repetitions(report):
    report["passed"] = False
    with pytest.raises(ValueError, match="scientific parity"):
        renderer.render(report)
    report["passed"] = True
    report["stages"]["paired_analysis"]["runs"] = 1
    with pytest.raises(ValueError, match="repetitions"):
        renderer.render(report)


@pytest.mark.parametrize(
    "text", ["none", renderer.END + renderer.START, renderer.START * 2 + renderer.END]
)
def test_rejects_missing_duplicate_or_reversed_markers(report, text):
    with pytest.raises(ValueError, match="marker"):
        renderer.update(text, report)


def test_cli_refuses_other_revision_without_changing_readme(report, tmp_path):
    source, readme = tmp_path / "report.json", tmp_path / "README.md"
    source.write_text(json.dumps(report))
    readme.write_text(renderer.START + renderer.END)
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "tools/update_benchmark_readme.py"),
            "--report",
            str(source),
            "--readme",
            str(readme),
            "--expected-sha",
            "b" * 40,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "revision does not match" in result.stderr
    assert readme.read_text() == renderer.START + renderer.END


def test_failed_timed_command_does_not_publish_report(tmp_path):
    output = tmp_path / "timing.json"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "tools/oracle/time_command.py"),
            "--report",
            str(output),
            "--",
            sys.executable,
            "-c",
            "raise SystemExit(7)",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize(
    "url",
    [
        "https://example.org/actions/runs/123",
        "https://github.com/o/r/actions/runs/0",
        "https://github.com/o/r/actions/runs/123)bad",
    ],
)
def test_rejects_invalid_report_link(report, url):
    report["run_url"] = url
    with pytest.raises(ValueError, match="CI run URL"):
        renderer.render(report)
