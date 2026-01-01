import click
from .workflow import WorkflowOrchestrator
from .database import AseDBWrapper
from ase.build import molecule
import os

@click.group()
def main():
    """MLIP-AutoPipe: Automated MLIP Generation Pipeline."""
    pass

@main.command()
@click.option('--id', 'uid', required=True, type=int, help="The ID of the structure to label.")
@click.option('--config', 'config_path', default='config.yaml', help="Path to the configuration file.")
def label(uid, config_path):
    """Run the DFT labeling for a single structure."""
    if not os.path.exists(config_path):
        click.echo(f"Error: Configuration file '{config_path}' not found.")
        return

    orchestrator = WorkflowOrchestrator(config_path)
    orchestrator.run_labeling_for_id(uid)

@main.command()
@click.option('--config', 'config_path', default='config.yaml', help="Path to the configuration file.")
def train(config_path):
    """Train the MLIP model from the labeled dataset."""
    if not os.path.exists(config_path):
        click.echo(f"Error: Configuration file '{config_path}' not found.")
        return

    orchestrator = WorkflowOrchestrator(config_path)
    orchestrator.run_training()

@main.command()
@click.option('--db', 'db_path', default='asedb.db', help="Path to the database file.")
def init_db(db_path):
    """
    Initializes the database with a sample water molecule.
    Deletes the existing database file if it exists.
    """
    if os.path.exists(db_path):
        os.remove(db_path)
        click.echo(f"Removed existing database at '{db_path}'.")

    db_wrapper = AseDBWrapper(db_path)
    h2o = molecule("H2O")
    h2o.set_cell([10, 10, 10]) # Set a cell for DFT calculation
    h2o.center()

    uid = db_wrapper.add_atoms(h2o)
    click.echo(f"Initialized database at '{db_path}'.")
    click.echo(f"Added H2O molecule with ID: {uid}")

if __name__ == "__main__":
    main()
