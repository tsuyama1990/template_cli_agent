import click
from .orchestrator import run_cycle01_workflow

@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the Cycle 01 YAML configuration file.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the initial atomic structure file (e.g., CIF, POSCAR).",
)
def main(config_path: str, structure_path: str):
    """
    CLI entrypoint for the MLIP-AutoPipe Cycle 01 workflow.
    """
    run_cycle01_workflow(config_path, structure_path)

if __name__ == "__main__":
    main()
