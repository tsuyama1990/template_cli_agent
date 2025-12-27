from pathlib import Path

import click

from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the YAML configuration file for Cycle 01.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the atomic structure file (e.g., CIF, XYZ).",
)
def main(config_path: Path, structure_path: Path):
    """
    CLI entry point for running the MLIP-AutoPipe Cycle 01 workflow.
    """
    run_cycle01_workflow(config_path, structure_path)


if __name__ == "__main__":
    main()
