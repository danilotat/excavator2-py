"""
Main CLI entry point for EXCAVATOR2.
"""

import click
from excavator2 import __version__


@click.group()
@click.version_option(version=__version__)
@click.option('--verbose', '-v', count=True, help='Increase verbosity')
@click.pass_context
def cli(ctx, verbose):
    """
    EXCAVATOR2: CNV detection from whole-exome sequencing data.

    This is a modern Python/C++ rewrite of the original EXCAVATOR2 tool.
    """
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose

    # Test that C++ module loads
    if verbose > 0:
        try:
            from excavator2._excavator_core import hello, has_openmp, __version__ as cpp_version
            click.echo(f"Python version: {__version__}")
            click.echo(f"C++ version: {cpp_version}")
            click.echo(f"OpenMP support: {has_openmp()}")
            click.echo(f"C++ module says: {hello()}")
        except ImportError as e:
            click.echo(f"Warning: C++ module not loaded: {e}", err=True)


@cli.command()
@click.option('--config', '-c', required=True, type=click.Path(exists=True),
              help='Path to configuration YAML file')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--force', '-f', is_flag=True, help='Overwrite existing output')
@click.pass_context
def target(ctx, config, output, force):
    """
    Initialize target regions (BED → HDF5 with MAP/GCC/FRB).

    Stage 1 of the EXCAVATOR2 pipeline.
    """
    click.echo("🚧 target command not yet implemented (Phase 6)")
    click.echo(f"  Config: {config}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Force: {force}")


@cli.command()
@click.option('--samples', '-s', required=True, type=click.Path(exists=True),
              help='Sample sheet YAML file')
@click.option('--target', '-t', required=True, type=click.Path(exists=True),
              help='Target HDF5 file (from target command)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--threads', '-@', default=1, type=int, help='Number of threads')
@click.option('--mapq', '-q', default=20, type=int, help='Minimum mapping quality')
@click.option('--reference', '-r', type=click.Path(exists=True),
              help='Reference FASTA (required for CRAM files)')
@click.option('--force', '-f', is_flag=True, help='Overwrite existing output')
@click.pass_context
def prepare(ctx, samples, target, output, threads, mapq, reference, force):
    """
    Read counting and normalization (BAM → NRC).

    Stage 2 of the EXCAVATOR2 pipeline. This command:

    1. Reads the sample sheet YAML mapping sample names to BAM/CRAM files
    2. Counts reads in each window defined by the target file
    3. Normalizes read counts for size, mappability, and GC biases
    4. Saves normalized counts to HDF5 files

    Example:

        excavator2 prepare -s samples.yaml -t target.h5 -o output/ -@ 4
    """
    import logging
    from pathlib import Path
    import yaml

    from excavator2.prepare import (
        ReadCountProcessor,
        ReadCountNormalizer,
        save_read_counts,
        save_normalized_counts,
    )

    # Setup logging
    verbose = ctx.obj.get('verbose', 0)
    log_level = logging.DEBUG if verbose > 1 else logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('excavator2.prepare')

    # Parse paths
    output_dir = Path(output)
    target_path = Path(target)

    # Check/create output directory
    if output_dir.exists() and not force:
        existing_files = list(output_dir.glob('*.h5'))
        if existing_files:
            raise click.ClickException(
                f"Output directory {output_dir} contains existing files. "
                "Use --force to overwrite."
            )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load sample sheet
    logger.info(f"Loading sample sheet: {samples}")
    with open(samples) as f:
        sample_sheet = yaml.safe_load(f)

    # Validate sample sheet format
    if not isinstance(sample_sheet, dict):
        raise click.ClickException(
            f"Invalid sample sheet format. Expected dictionary mapping "
            "sample names to BAM paths."
        )

    click.echo(f"Processing {len(sample_sheet)} samples")
    click.echo(f"Target file: {target_path}")
    click.echo(f"Output directory: {output_dir}")
    click.echo(f"MAPQ threshold: {mapq}")
    click.echo(f"Threads: {threads}")

    # Initialize processor
    try:
        processor = ReadCountProcessor(
            target_path=target_path,
            min_mapq=mapq,
            reference=reference
        )
    except FileNotFoundError as e:
        raise click.ClickException(str(e))
    except ValueError as e:
        raise click.ClickException(f"Invalid target file: {e}")

    click.echo(f"Loaded {len(processor.windows)} windows from target file")

    # Initialize normalizer
    normalizer = ReadCountNormalizer()

    # Process each sample
    for sample_name, bam_path in sample_sheet.items():
        click.echo(f"\nProcessing sample: {sample_name}")

        # Validate BAM path
        bam_path = Path(bam_path)
        if not bam_path.exists():
            click.echo(f"  WARNING: BAM file not found: {bam_path}", err=True)
            continue

        try:
            # Count reads
            click.echo(f"  Counting reads from: {bam_path}")
            sample_data = processor.process_sample(bam_path, sample_name)
            click.echo(f"  Total reads: {sample_data.total_reads:,}")
            click.echo(f"  Mean count per window: {sample_data.raw_counts.mean():.2f}")

            # Save raw counts
            raw_output = output_dir / f"{sample_name}.RC.h5"
            save_read_counts(sample_data, raw_output)
            click.echo(f"  Saved raw counts: {raw_output}")

            # Normalize
            click.echo("  Normalizing read counts...")
            norm_result = normalizer.normalize(sample_data)
            click.echo(f"  Normalized mean: {norm_result.normalized_counts.mean():.2f}")

            # Save normalized counts
            norm_output = output_dir / f"{sample_name}.NRC.h5"
            save_normalized_counts(norm_result, norm_output)
            click.echo(f"  Saved normalized counts: {norm_output}")

        except Exception as e:
            click.echo(f"  ERROR: Failed to process {sample_name}: {e}", err=True)
            if verbose > 0:
                import traceback
                traceback.print_exc()
            continue

    click.echo("\nData preparation complete.")


@cli.command()
@click.option('--samples', '-s', required=True, type=click.Path(exists=True),
              help='Sample file list YAML')
@click.option('--input', '-i', required=True, type=click.Path(exists=True),
              help='Input directory (from prepare command)')
@click.option('--target', '-t', required=True, type=click.Path(exists=True),
              help='Target directory')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--experiment', '-e', required=True,
              type=click.Choice(['paired', 'pooled']),
              help='Experimental design')
@click.option('--parameters', '-p', type=click.Path(exists=True),
              help='Parameters YAML (optional)')
@click.option('--threads', '-@', default=1, type=int)
@click.option('--force', '-f', is_flag=True)
@click.pass_context
def analyze(ctx, samples, input, target, output, experiment, parameters, threads, force):
    """
    Segmentation and CNV calling (HSLM + FastCall).

    Stage 3 of the EXCAVATOR2 pipeline.
    """
    click.echo("🚧 analyze command not yet implemented (Phase 7)")
    click.echo(f"  Samples: {samples}")
    click.echo(f"  Input: {input}")
    click.echo(f"  Target: {target}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Experiment: {experiment}")


if __name__ == '__main__':
    cli()
