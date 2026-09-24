"""Create small deterministic BAMs from the shared legacy fixture's count profiles."""

from pathlib import Path

import numpy as np
import pysam
import yaml

from excavator2.artifacts import read_export


def make_bams(output):
    output = Path(output)
    folder = output / "bams"
    folder.mkdir()
    samples = {}
    for exported in sorted((output / "exports/prepared").glob("*/RCNorm/*.NRC.RData/MatrixNorm")):
        name = exported.parents[2].name
        matrix = read_export(exported)
        chromosomes = list(dict.fromkeys(matrix[:, 0]))
        path = folder / f"{name}.bam"
        header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": c, "LN": 100000} for c in chromosomes]}
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for i, row in enumerate(matrix):
                count = int(np.rint(float(row[5])))
                for j in range(count + 2):
                    read = pysam.AlignedSegment()
                    read.query_name = f"r{i}_{j}"
                    read.flag = (
                        [0, 256, 2048, 64, 128][j % 5] if j < count else [4, 1024][j - count]
                    )
                    read.reference_id = chromosomes.index(row[0])
                    read.reference_start = int(float(row[2]))
                    read.mapping_quality = [0, 1, 19, 20, 60][j % 5]
                    read.query_sequence = "A"
                    read.cigarstring = "1M"
                    bam.write(read)
        pysam.index(str(path))
        samples[name] = str(path.resolve())
    (output / "bam-samples.yaml").write_text(yaml.safe_dump(samples, sort_keys=False))
    return samples
