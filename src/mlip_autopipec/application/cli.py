"""Main CLI entry point for the MLIP AutoPipe."""

import click
import os
import yaml
import traceback
from ase.io import read as ase_read
from typing import List
from ase import Atoms

from mlip_autopipec.config import Settings
from mlip_autopipec.application.services import LabellingService, TrainingService
from mlip_autopipec.infrastructure.database import AseDB
from mlip_autopipec.infrastructure.process_runner import SubprocessRunner

class ApplicationError(Exception):
    """Base exception for application-level errors."""
    pass


@click.group()
def cli():
    """MLIP-AutoPipe command-line interface."""
    pass


def load_structures_from_directory(input_dir: str) -> List[Atoms]:
    """Loads all supported atomic structure files from a directory."""
    supported_exts = ('.cif', '.poscar', '.xyz')
    files = [f for f in os.listdir(input_dir) if f.endswith(supported_exts)]
    if not files:
        raise ApplicationError(f"No structure files found in directory: {input_dir}")
    return [ase_read(os.path.join(input_dir, f)) for f in files]


@cli.command("run-cycle-01")
@click.option("--config-file", default="config.yaml", help="Path to the configuration file.")
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory containing initial structure files.")
def run_cycle_01(config_file: str, input_dir: str):
    """
    Runs the full Cycle 01 workflow: data loading, labelling, and training.
    """
    try:
        click.echo("--- Starting MLIP-AutoPipe Cycle 01 ---")

        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        settings = Settings(**config_data)

        db_adapter = AseDB(db_path=settings.db_path)
        runner_adapter = SubprocessRunner()

        labelling_service = LabellingService(db=db_adapter, runner=runner_adapter, settings=settings)
        training_service = TrainingService(db=db_adapter, settings=settings)

        if os.path.exists(settings.db_path):
            os.remove(settings.db_path)

        initial_atoms = load_structures_from_directory(input_dir)
        db_adapter.add_atoms(initial_atoms)
        click.echo(f"Loaded {len(initial_atoms)} structures into the database.")

        click.echo("\n--- Running Labelling Service ---")
        labelling_service.run()
        click.echo("--- Labelling Service Finished ---")

        click.echo("\n--- Running Training Service ---")
        training_service.run()
        click.echo("--- Training Service Finished ---")

        click.echo("\n--- MLIP-AutoPipe Cycle 01 Finished Successfully ---")

    except ApplicationError as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"An unexpected error occurred: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    cli()
