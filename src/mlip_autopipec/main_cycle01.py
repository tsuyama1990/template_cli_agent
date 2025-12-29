# ruff: noqa: D101, D102, D103, D104, D105, D107
from pathlib import Path

import click

from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the Cycle 01 YAML configuration file.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the initial atomic structure file (e.g., CIF, POSCAR).",
)
@click.option(
    "--db",
    "db_path",
    default="mlip.db",
    type=str,
    show_default=True,
    help="Path to the ASE database file.",
)
def main(config_path: Path, structure_path: Path, db_path: str):
    """
    CLI for running the MLIP-AutoPipe Cycle 01 workflow.
    This workflow takes a single atomic structure, calculates its DFT properties
    using Quantum Espresso, and trains a basic MACE model on the resulting data.
    """
    run_cycle01_workflow(config_path, structure_path, db_path)


if __name__ == "__main__":
    main()
