"""Render validated CI timing JSON into a bounded README section."""

import argparse
import json
import math
import re
from pathlib import Path

START = "<!-- ci-benchmark:start -->"
END = "<!-- ci-benchmark:end -->"


def render(report):
    if report.get("passed") is not True:
        raise ValueError("cannot publish a benchmark without scientific parity")
    revision, url = report["source_commit"], report["run_url"]
    if not re.fullmatch(r"[0-9a-f]{40}", revision) or not re.fullmatch(
        r"https://github.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+/actions/runs/[1-9][0-9]*", url
    ):
        raise ValueError("expected an exact revision and a GitHub CI run URL")
    if type(report["windows"]) is not int or report["windows"] < 1:
        raise ValueError("window count must be positive")
    lines = [
        f"Latest passing [CI measurement]({url}) · revision `{revision[:7]}` · "
        f"{report['windows']:,} synthetic windows · one analysis test sample.",
        "",
        "| Stage | Original | Python/C++ | Original / port | Measured runs |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for key, label in [("target", "Target generation"), ("paired_analysis", "Paired analysis")]:
        stage = report["stages"][key]
        old, new = stage["original_seconds"], stage["port_seconds"]
        if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in [old, new]):
            raise ValueError("timings must be positive finite numbers")
        count = stage["runs"]
        if type(count) is not int or count < (3 if key == "paired_analysis" else 1):
            raise ValueError("insufficient benchmark repetitions")
        lines.append(f"| {label} | {old:.3f} s | {new:.3f} s | {old / new:.2f}× | {count} |")
    lines += [
        "",
        "Analysis values are medians after one excluded warm-up. "
        "Ratios above 1 mean the port was faster.",
    ]
    return "\n".join(lines)


def update(text, report):
    if text.count(START) != 1 or text.count(END) != 1:
        raise ValueError("README must contain exactly one benchmark marker pair")
    start, end = text.index(START) + len(START), text.index(END)
    if end < start:
        raise ValueError("benchmark markers are out of order")
    return text[:start] + "\n\n" + render(report) + "\n\n" + text[end:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--readme", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--expected-sha", required=True)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    if report.get("source_commit") != args.expected_sha:
        parser.error("benchmark revision does not match the expected commit")
    rendered = render(report)
    if args.readme:
        args.readme.write_text(update(args.readme.read_text(), report))
    if args.summary:
        with args.summary.open("a") as handle:
            handle.write("## Original versus port benchmark\n\n" + rendered + "\n")
