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
    Initialize target regions (BED -> HDF5 with MAP/GCC/FRB).

    Stage 1 of the EXCAVATOR2 pipeline. This command:

    1. Reads the configuration YAML with Reference and Target sections
    2. Creates analysis windows from target BED regions
    3. Filters windows overlapping genomic gaps (centromeres, telomeres)
    4. Calculates GC content from reference FASTA
    5. Extracts mappability from BigWig file
    6. Saves all data to HDF5 format

    Example:

        excavator2 target -c config.yaml -o output/

    Config file format:

    \b
        Reference:
          Assembly: hg38
          FASTA: /path/to/ref.fasta
          BigWig: /path/to/mappability.bw
          Chromosomes: /path/to/chromosomes.txt
          Gaps: /path/to/gaps.txt

        Target:
          Name: SureSelectV7
          BED: /path/to/targets.bed
          Window: 30000
    """
    import logging
    from pathlib import Path
    import yaml

    from excavator2.target import (
        initialize_target,
        save_target_data,
    )

    # Setup logging
    verbose = ctx.obj.get('verbose', 0)
    log_level = logging.DEBUG if verbose > 1 else logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('excavator2.target')

    # Parse output path
    output_dir = Path(output)

    # Load configuration
    logger.info(f"Loading configuration: {config}")
    with open(config) as f:
        cfg = yaml.safe_load(f)

    # Validate configuration sections
    if 'Reference' not in cfg:
        raise click.ClickException("Configuration missing 'Reference' section")
    if 'Target' not in cfg:
        raise click.ClickException("Configuration missing 'Target' section")

    ref_cfg = cfg['Reference']
    target_cfg = cfg['Target']

    # Validate required Reference fields
    required_ref = ['Assembly', 'FASTA', 'BigWig', 'Chromosomes', 'Gaps']
    for field in required_ref:
        if field not in ref_cfg:
            raise click.ClickException(f"Reference section missing '{field}' field")

    # Validate required Target fields
    required_target = ['Name', 'BED', 'Window']
    for field in required_target:
        if field not in target_cfg:
            raise click.ClickException(f"Target section missing '{field}' field")

    # Extract configuration values
    assembly = ref_cfg['Assembly']
    fasta_path = Path(ref_cfg['FASTA'])
    bigwig_path = Path(ref_cfg['BigWig'])
    chromosome_path = Path(ref_cfg['Chromosomes'])
    gap_path = Path(ref_cfg['Gaps'])
    target_name = target_cfg['Name']
    bed_path = Path(target_cfg['BED'])
    window_size = int(target_cfg['Window'])

    # Validate window size
    if window_size < 10:
        raise click.ClickException(f"Window size must be at least 10bp, got {window_size}")

    # Validate input files exist
    for path, name in [(fasta_path, 'FASTA'), (bigwig_path, 'BigWig'),
                       (chromosome_path, 'Chromosomes'), (gap_path, 'Gaps'),
                       (bed_path, 'BED')]:
        if not path.exists():
            raise click.ClickException(f"{name} file not found: {path}")

    # Create output directory structure
    target_output_dir = output_dir / assembly / target_name / f"w_{window_size}"

    if target_output_dir.exists() and not force:
        existing_files = list(target_output_dir.glob('*.h5'))
        if existing_files:
            raise click.ClickException(
                f"Output directory {target_output_dir} contains existing files. "
                "Use --force to overwrite."
            )

    target_output_dir.mkdir(parents=True, exist_ok=True)

    # Display configuration
    click.echo(f"Assembly: {assembly}")
    click.echo(f"Target: {target_name}")
    click.echo(f"Window size: {window_size}")
    click.echo(f"Output: {target_output_dir}")
    click.echo(f"BED file: {bed_path}")
    click.echo(f"FASTA file: {fasta_path}")
    click.echo(f"BigWig file: {bigwig_path}")

    # Initialize target
    try:
        click.echo("\nInitializing target regions...")
        target_data = initialize_target(
            bed_path=bed_path,
            fasta_path=fasta_path,
            bigwig_path=bigwig_path,
            chromosome_path=chromosome_path,
            gap_path=gap_path,
            target_name=target_name,
            assembly=assembly,
            window_size=window_size
        )

        click.echo(f"\nTarget statistics:")
        click.echo(f"  Total windows: {target_data.n_windows}")
        click.echo(f"  IN-target windows: {target_data.n_in_target}")
        click.echo(f"  OUT-target windows: {target_data.n_out_target}")
        click.echo(f"  Chromosomes: {len(target_data.chromosomes)}")

        # Save to HDF5
        output_file = target_output_dir / f"{target_name}.h5"
        click.echo(f"\nSaving target data to: {output_file}")
        save_target_data(target_data, output_file)

        # Also save a settings file for reference
        settings_file = target_output_dir / "settings.yaml"
        settings = {
            'Reference': {
                'Assembly': assembly,
                'FASTA': str(fasta_path),
                'BigWig': str(bigwig_path),
                'Chromosomes': str(chromosome_path),
                'Gaps': str(gap_path)
            },
            'Target': {
                'Name': target_name,
                'BED': str(bed_path),
                'Window': window_size
            },
            'Output': {
                'n_windows': target_data.n_windows,
                'n_in_target': target_data.n_in_target,
                'n_out_target': target_data.n_out_target,
                'n_chromosomes': len(target_data.chromosomes)
            }
        }
        with open(settings_file, 'w') as f:
            yaml.dump(settings, f, default_flow_style=False)
        click.echo(f"Saved settings to: {settings_file}")

        click.echo("\nTarget initialization complete.")

    except Exception as e:
        logger.error(f"Target initialization failed: {e}")
        if verbose > 0:
            import traceback
            traceback.print_exc()
        raise click.ClickException(str(e))


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
