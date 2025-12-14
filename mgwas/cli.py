import click
import logging
import os
from mgwas import core, utils, version

@click.command()
@click.option('-i', '--input', 'input_file', type=click.Path(exists=True), required=True,
              help="Tab-delimited file with shotgun data and metadata")
@click.option('-g', '--genome', required=True,
              help="NCBI RefSeq species name or custom reference genome FASTA path")
@click.option('-k', '--genbank', default=None,
              help="Path to gene annotations in GenBank format for custom genome")
@click.option('-o', '--output', default='./out/', show_default=True,
              help="Output folder path")
@click.option('-f', '--force', is_flag=True, default=False, show_default=True,
              help="Force overwrite outputs")
@click.option('-m', '--mode', type=click.Choice(['variant', 'gene', 'both']), default='both', show_default=True,
              help="GWAS testing mode")
@click.option('--filter', 'filter_level', type=click.Choice(["-1", "0", "1", "2", "3"]), default="-1", show_default=True,
              help="Impact level to filter variants on based on SnpEff annotations. 3: high impact; 2: moderate impact and above; 1: low impact and above; 0: modifier and above; -1: no filtering")
@click.option('--paired', is_flag=True, default=False, show_default=True,
              help="Flag for paired-end reads")
@click.option('--temp', default='./tmp/', show_default=True,
              help="Temporary folder path")
@click.option('--keep-temp', is_flag=True, default=False, help="Keep temporary folder after run.")
@click.option('--ploidy', type=int, default=2, show_default=True,
              help="Sample ploidy parameter for GATK HaplotypeCaller")
@click.option('--min-depth', type=int, default=1, show_default=True,
              help="Minimum depth cutoff")
@click.option('--min-maf', type=float, default=0.1, show_default=True,
              help="Minimum minor allele frequency cutoff ranging from 0 to 1")
@click.option('--min-sample-size', type=int, default=10, show_default=True,
              help="Minimum sample size cutoff")
@click.option('-q', '--qcutoff', type=float, default=0.05, show_default=True,
              help="Q-value cutoff for calling significant hits")
@click.option('-t', '--threads', type=int, default=4, show_default=True,
              help="Number of processors to use")
@click.version_option(version.__version__, '-v', '--version', prog_name='mgwas')
def main(**kwargs):
    # Basic setup
    os.makedirs(kwargs['output'], exist_ok=True)
    if not kwargs['force'] and any(os.scandir(kwargs['output'])):
        click.echo(f"Output folder \"{kwargs['output']}\" exists and is not empty. Use --force to overwrite.")
        return
    else:
        log_file = os.path.join(kwargs['output'], 'logfile')
        utils.setup_logging(log_file)
        click.echo(f"Log info will be written to {log_file}.")

    # Check required software
    missing_tools = utils.check_dependencies(['bowtie2', 'samtools', 'gatk'])
    if missing_tools:
        click.echo(f"Missing required tools: {', '.join(missing_tools)}. Please install them before running.")
        return
    
    # Run pipeline
    core.run_pipeline(**kwargs)
    
    final_message = "MGWAS finished successfully!"
    logging.info(final_message)
    click.secho(f"\n{final_message}\n", fg="green", bold=True)
