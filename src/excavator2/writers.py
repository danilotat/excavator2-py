"""Legacy table and VCF layout, including historical rounding and grouping."""

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
        ("CN", "Integer", "Copy number genotype for imprecise events"),
        ("CNF", "Float", "Copy number genotype fraction for imprecise events"),
        (
            "FCL",
            "Float",
            "FastCall inference Label" if windows else "Label Inferred by FastCall algorithm",
        ),
        (
            "FCP",
            "Float",
            "FastCall Posterior Probability"
            if windows
            else "Posterior Probability inferred by FastCall algorithm",
        ),
    ]
    if windows:
        formats.append(("L2R", "Float", "Normalized log2Ratio value"))
    for field, kind, description in formats:
        lines.append(f'##FORMAT=<ID={field},Number=1,Type={kind},Description="{description}">')
    lines.append(
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample])
    )
    return lines


def write_results(folder, sample, rows, starts, ends, original, calls, target, target_folder):
    header = "Chromosome\tPosition\tStart\tEnd\tLog2R\tSegMean\tClass"
    with (folder / f"HSLMResults_{sample}.txt").open("w") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    significant = np.flatnonzero(calls.labels != 0)
    with (folder / f"FastCallResults_{sample}.txt").open("w") as handle:
        handle.write("Chromosome\tStart\tEnd\tSegment\tCNF\tCN\tCall\tProbCall\n")
        for index in significant:
            start, end = starts[index], ends[index] - 1
            cnf = 2 * 2 ** original[index]
            handle.write(
                "\t".join(
                    [
                        rows[start, 0],
                        rows[start, 2],
                        rows[end, 3],
                        number(original[index]),
                        number(cnf),
                        number(np.rint(cnf)),
                        str(calls.labels[index]),
                        number(calls.probabilities[index]),
                    ]
                )
                + "\n"
            )
    for windows in [True, False]:
        records = []
        for index in significant:
            start, end = starts[index], ends[index] - 1
            cnf = np.round(2 * 2 ** original[index], 2)
            genotype = (
                f"1/1:{number(np.rint(cnf))}:{number(cnf)}:{calls.labels[index]}:"
                f"{number(np.round(calls.probabilities[index], 2))}"
            )
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
            wanted = {record[1] for record in group}
            # Legacy selects matching FRB rows in FRB order, not by reordering lookup.
            bases = reference[np.isin(reference[:, 0], list(wanted)), 1]
            if len(bases) != len(group):
                raise ValueError("legacy reference-base selection does not match call rows")
            for (_, start, end, genotype), base in zip(group, bases, strict=True):
                start, end = str(int(float(start))), str(int(float(end)))
                info = f"IMPRECISE;SVTYPE=CNV;END={end};SVLEN={int(end) - int(start) + 1};"
                fields = "GT:CN:CNF:FCL:FCP" + (":L2R" if windows else "")
                lines.append(
                    "\t".join(
                        [chromosome, start, ".", base, "<CNV>", ".", "PASS", info, fields, genotype]
                    )
                )
        kind = "Window" if windows else "Region"
        (folder / f"EXCAVATOR{kind}Call_{sample}.vcf").write_text("\n".join(lines) + "\n")
