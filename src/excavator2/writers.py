"""Tables and VCFs with explicit corrected diploid-equivalent semantics."""

from datetime import datetime
from decimal import Decimal

import numpy as np


def number(value):
    value = float(value)
    text = format(value, ".15g")
    fixed = format(Decimal(text), "f")
    mantissa, exponent = format(value, ".14e").split("e")
    scientific = mantissa.rstrip("0").rstrip(".") + "e" + exponent
    return scientific if len(scientific) < len(fixed) else fixed


def vcf_header(sample, assembly, windows):
    lines = [
        "##fileformat=VCFv4.0",
        datetime.now().strftime("##fileDate=%Y%d%m"),
        "##source=EXCAVATOR2v1.0",
        f"##reference={assembly}",
        "##assembly=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta",
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,"
        'Description="Difference in length between REF and ALT alleles">',
        "##INFO=<ID=END,Number=1,Type=Integer,"
        'Description="End position of the variant described in this record">',
        '##ALT=<ID=CNV,Description="Copy number variable region">',
    ]
    formats = [
        ("GT", "String", "Genotype"),
        ("DECN", "Integer", "Rounded cellularity-corrected diploid-equivalent copy number"),
        ("DECNF", "Float", "Cellularity-corrected diploid-equivalent copy number; not absolute CN"),
        (
            "FCL",
            "Float",
            "FastCall inference Label" if windows else "Label Inferred by FastCall algorithm",
        ),
        (
            "FCP",
            "Float",
            "Probability that the segment level lies in the called class interval, "
            "from segment measurement uncertainty; not calibrated event confidence",
        ),
        (
            "FPA",
            "Float",
            "Probability that the segment level lies outside the normal interval, "
            "from segment measurement uncertainty; not calibrated event confidence",
        ),
    ]
    if windows:
        formats.append(
            ("L2R", "Float", "Observed mixture log2 ratio before cellularity correction")
        )
    for field, kind, description in formats:
        lines.append(f'##FORMAT=<ID={field},Number=1,Type={kind},Description="{description}">')
    lines.append(
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample])
    )
    return lines


def call_summary(probabilities, labels):
    """Probabilities of the called class and of any non-normal class."""
    called = probabilities[np.arange(len(labels)), labels + 2]
    return called, probabilities[:, [0, 1, 3, 4]].sum(axis=1)


def genotypes(corrected, labels, probabilities):
    called, altered = call_summary(probabilities, labels)
    cnf = 2 * 2 ** np.asarray(corrected)
    return [
        f".:{number(np.rint(c))}:{number(c)}:{label}:"
        f"{number(np.round(p, 2))}:{number(np.round(a, 2))}"
        for c, label, p, a in zip(cnf, labels, called, altered, strict=True)
    ]


def write_results(
    folder,
    sample,
    rows,
    starts,
    ends,
    original,
    corrected,
    calls,
    probabilities,
    events,
    target,
    target_folder,
):
    """Segments in tables and the Window VCF; merged events in the Region VCF."""
    header = "Chromosome\tPosition\tStart\tEnd\tLog2R\tSegMean\tClass"
    with (folder / f"HSLMResults_{sample}.txt").open("w") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    significant = np.flatnonzero(calls.labels != 0)
    called, altered = call_summary(probabilities, calls.labels)
    with (folder / f"FastCallResults_{sample}.txt").open("w") as handle:
        handle.write(
            "Chromosome\tStart\tEnd\tMixtureSegment\tCorrectedSegment"
            "\tDECNF\tDECN\tCall\tProbCall\tProbAltered\tEventId\n"
        )
        for index in significant:
            start, end = starts[index], ends[index] - 1
            cnf = 2 * 2 ** corrected[index]
            handle.write(
                "\t".join(
                    [
                        rows[start, 0],
                        rows[start, 2],
                        rows[end, 3],
                        number(original[index]),
                        number(corrected[index]),
                        number(cnf),
                        number(np.rint(cnf)),
                        str(calls.labels[index]),
                        number(called[index]),
                        number(altered[index]),
                        str(events.ids[index]),
                    ]
                )
                + "\n"
            )
    segment_genotypes = genotypes(corrected, calls.labels, probabilities)
    event_genotypes = genotypes(events.corrected, events.labels, events.probabilities)
    for windows in [True, False]:
        records = []
        if windows:
            intervals = [(starts[i], ends[i] - 1, segment_genotypes[i]) for i in significant]
        else:
            intervals = [
                (events.starts[i], events.ends[i] - 1, event_genotypes[i])
                for i in np.flatnonzero(events.labels != 0)
            ]
        for start, end, genotype in intervals:
            selected = range(start, end + 1) if windows else [start]
            for row in selected:
                records.append(
                    (
                        rows[row, 0],
                        rows[row, 2],
                        rows[row if windows else end, 3],
                        genotype + (":" + rows[row, 4] if windows else ""),
                    )
                )
        lines = vcf_header(sample, target["assembly"], windows)
        for chromosome in dict.fromkeys(record[0] for record in records):
            group = [record for record in records if record[0] == chromosome]
            filename = target["references"][chromosome]
            if filename not in target["files"]:
                raise ValueError("reference file lacks a manifest checksum")
            with np.load(target_folder / filename, allow_pickle=False) as archive:
                reference = archive["matrix"]
            # Match numerically: rows may render round positions as e.g. "1e+05".
            wanted = {int(float(record[1])) for record in group}
            # Legacy selects matching FRB rows in FRB order, not by reordering lookup.
            bases = reference[np.isin(reference[:, 0].astype(int), list(wanted)), 1]
            if len(bases) != len(group):
                raise ValueError("legacy reference-base selection does not match call rows")
            for (_, start, end, genotype), base in zip(group, bases, strict=True):
                start, end = str(int(float(start))), str(int(float(end)))
                info = f"IMPRECISE;SVTYPE=CNV;END={end};SVLEN={int(end) - int(start) + 1};"
                fields = "GT:DECN:DECNF:FCL:FCP:FPA" + (":L2R" if windows else "")
                lines.append(
                    "\t".join(
                        [chromosome, start, ".", base, "<CNV>", ".", "PASS", info, fields, genotype]
                    )
                )
        kind = "Window" if windows else "Region"
        (folder / f"EXCAVATOR{kind}Call_{sample}.vcf").write_text("\n".join(lines) + "\n")
