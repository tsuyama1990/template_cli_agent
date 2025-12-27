import click

from .orchestrator_cycle01 import run_cycle01_workflow

@click.command()
@click.option(
    '--config',
    'config_path',
    required=True,
    type=click.Path(exists=True),
    help="Path to the main YAML configuration file."
)
@click.option(
    '--structure',
    'structure_path',
    required=True,
    type=click.Path(exists=True),
    help="Path to the input atomic structure file (e.g., XYZ, CIF)."
)
def main(config_path: str, structure_path: str):
    """
    Command-line interface for running the MLIP-AutoPipe Cycle 01 workflow.
    """
    run_cycle01_workflow(config_path=config_path, structure_path=structure_path)

if __name__ == "__main__":
    main()
