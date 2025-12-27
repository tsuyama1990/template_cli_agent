"""
Main CLI entry point for the MLIP-AutoPipe application.
"""
import click
from pathlib import Path

from mlip_autopipec.orchestrator import run_workflow

@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the YAML configuration file.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the input atomic structure file (e.g., CIF, POSCAR).",
)
def main(config_path: Path, structure_path: Path):
    """
    Runs the MLIP-AutoPipe Cycle 01 workflow: DFT Labelling and MLIP Training.
    """
    run_workflow(config_path=config_path, structure_path=structure_path)

if __name__ == "__main__":
    main()
