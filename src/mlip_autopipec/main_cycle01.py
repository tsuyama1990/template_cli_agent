import os

import click

from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the Cycle 01 YAML configuration file."
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the atomic structure file (e.g., XYZ, CIF)."
)
def main(config_path: str, structure_path: str):
    """
    CLI for running the MLIP-AutoPipe Cycle 01 workflow.

    This script takes a configuration file and an atomic structure file
    as input, runs the DFT labelling, and then trains a MACE model.
    """
    # Ensure paths are absolute to avoid issues with changing directories
    abs_config_path = os.path.abspath(config_path)
    abs_structure_path = os.path.abspath(structure_path)

    run_cycle01_workflow(
        config_path=abs_config_path,
        structure_path=abs_structure_path
    )

if __name__ == "__main__":
    main()
