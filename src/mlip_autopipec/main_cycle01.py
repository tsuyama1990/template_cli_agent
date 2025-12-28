import click

from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the YAML configuration file for Cycle 01.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the input atomic structure file (e.g., CIF, XYZ).",
)
def main(config_path: str, structure_path: str):
    """
    Command-line interface for running the MLIP-AutoPipe Cycle 01 workflow.

    This script takes an atomic structure and a configuration file,
    runs a DFT calculation to generate labels (energy, forces, stress),
    and then trains a basic Machine Learning Interatomic Potential (MLIP)
    using the MACE framework.
    """
    run_cycle01_workflow(config_path, structure_path)


if __name__ == "__main__":
    main()
