"""Installation and CLI safety checks; these do not establish scientific parity."""

import subprocess
import sys

import pytest

from excavator2 import __version__, _core


def test_compiled_extension_imports():
    assert callable(_core.fastcall_posterior)


def test_module_version():
    result = subprocess.run(
        [sys.executable, "-m", "excavator2", "--version"], capture_output=True, text=True
    )
    assert result.returncode == 0
    assert result.stdout.strip() == f"excavator2 {__version__}"


@pytest.mark.parametrize(
    "command",
    [
        ["target", "--config", "config.yaml"],
        ["prepare", "--samples", "samples.yaml", "--target", "target/", "--threads", "4"],
    ],
)
def test_unimplemented_stage_fails_without_creating_output(tmp_path, command):
    output = tmp_path / "output"
    result = subprocess.run(
        [sys.executable, "-m", "excavator2", *command, "--output", str(output)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "not implemented" in result.stderr
    assert not output.exists()
