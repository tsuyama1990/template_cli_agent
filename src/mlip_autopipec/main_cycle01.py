# Description: CLI entry point for the Cycle 01 workflow.
import click

from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the YAML configuration file.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the atomic structure file (e.g., CIF, XYZ).",
)
def main(config_path: str, structure_path: str):
    """
    Runs the full DFT labelling and MLIP training workflow for Cycle 01.
    """
    run_cycle01_workflow(config_path, structure_path)


if __name__ == "__main__":
    main()
