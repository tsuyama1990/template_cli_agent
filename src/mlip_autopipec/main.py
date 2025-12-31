"""Main CLI entry point for the MLIP AutoPipe."""

import os

import click
import yaml
from ase.io import read as ase_read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.labelling_engine import LabellingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.utils.runner import SubprocessRunner


@click.group()
def cli():
    """MLIP-AutoPipe command-line interface."""
    pass


@cli.command("run-cycle-01")
@click.option("--config-file", default="config.yaml", help="Path to the configuration file.")
@click.option("--input-dir", default=".", help="Directory containing initial structure files.")
def run_cycle_01(config_file, input_dir):
    """
    Runs the full Cycle 01 workflow: data loading, labelling, and training.
    """
    click.echo("--- Starting MLIP-AutoPipe Cycle 01 ---")

    # Load configuration
    with open(config_file) as f:
        config = yaml.safe_load(f)

    db_path = config.get("db_path", "cycle01.db")

    # 1. Initialize Database and load initial data
    if os.path.exists(db_path):
        os.remove(db_path)  # Start with a clean DB for the cycle
    db = AseDB(db_path)

    structure_files = [f for f in os.listdir(input_dir) if f.endswith((".cif", ".poscar", ".xyz"))]
    initial_atoms = [ase_read(os.path.join(input_dir, f)) for f in structure_files]

    if not initial_atoms:
        click.echo("No initial structures found. Exiting.")
        return

    db.add_atoms(initial_atoms)
    click.echo(f"Loaded {len(initial_atoms)} structures into the database.")

    # 2. Run Labelling Engine
    click.echo("\n--- Running Labelling Engine ---")
    runner = SubprocessRunner()
    labelling_engine = LabellingEngine(db, runner, config)
    labelling_engine.run()
    click.echo("--- Labelling Engine Finished ---")

    # 3. Run Training Engine
    click.echo("\n--- Running Training Engine ---")
    training_engine = TrainingEngine(db, config)
    training_engine.run()
    click.echo("--- Training Engine Finished ---")

    click.echo("\n--- MLIP-AutoPipe Cycle 01 Finished Successfully ---")


if __name__ == "__main__":
    cli()
