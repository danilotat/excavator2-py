"""Time a checked subprocess, excluding its surrounding container/harness setup."""

import argparse
import json
import subprocess
import time
from pathlib import Path


def measure(command, **kwargs):
    start = time.perf_counter()
    subprocess.run(command, check=True, **kwargs)
    return time.perf_counter() - start


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command is required")
    args.report.write_text(json.dumps({"seconds": measure(command)}) + "\n")
