"""Regenerate pure-R HSLM boundary evidence using the pinned original source."""

import argparse
import subprocess
from pathlib import Path

from oracle import BASELINE, ROOT, container_command, file_record, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare-fixtures", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    # Docker writes root-owned files on Linux. Keep the directory host-owned
    # so the invoking user can add the provenance manifest after the run.
    (output / "edges").mkdir()
    source = output / "LibraryJSLMIn.R"
    source.write_bytes(
        subprocess.check_output(
            ["git", "show", f"{BASELINE['source_commit']}:excavator2/lib/R/LibraryJSLMIn.R"],
            cwd=ROOT,
        )
    )
    script = ROOT / "tools/oracle/characterize_hslm_edges.R"
    with (output / "legacy.log").open("w") as log:
        subprocess.run(
            container_command(output)
            + [
                "Rscript",
                f"/repo/{script.relative_to(ROOT)}",
                "/work/LibraryJSLMIn.R",
                "/work/edges",
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    files = {p.name: file_record(p) for p in sorted((output / "edges").iterdir())}
    if args.compare_fixtures:
        fixtures = ROOT / "tests/fixtures/legacy-hslm/edges"
        for name in files:
            if (output / "edges" / name).read_bytes() != (fixtures / name).read_bytes():
                raise ValueError(f"legacy HSLM boundary fixture differs: {name}")
    save_json(
        output / "edges/manifest.json",
        {
            "baseline": BASELINE,
            "source": file_record(source),
            "harness": file_record(script),
            "files": files,
        },
    )


if __name__ == "__main__":
    main()
