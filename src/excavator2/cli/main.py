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
              help='Target directory (from target command)')
@click.option('--output', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--threads', '-@', default=1, type=int, help='Number of threads')
@click.option('--mapq', '-q', default=20, type=int, help='Minimum mapping quality')
@click.option('--force', '-f', is_flag=True)
@click.pass_context
def prepare(ctx, samples, target, output, threads, mapq, force):
    """
    Read counting and normalization (BAM → NRC).

    Stage 2 of the EXCAVATOR2 pipeline.
    """
    click.echo("🚧 prepare command not yet implemented (Phase 5)")
    click.echo(f"  Samples: {samples}")
    click.echo(f"  Target: {target}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")


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
