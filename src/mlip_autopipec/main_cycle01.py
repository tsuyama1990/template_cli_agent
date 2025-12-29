import click
import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


@click.command()
@click.option(
    '--config',
    'config_path',
    required=True,
    type=click.Path(exists=True),
    help='Path to the YAML configuration file.'
)
@click.option(
    '--structure',
    'structure_path',
    required=True,
    type=click.Path(exists=True),
    help='Path to the atomic structure file (e.g., XYZ, CIF).'
)
def run_cycle01_workflow(config_path, structure_path):
    """
    Runs the simple, linear workflow for Cycle 01:
    1. Reads a structure file.
    2. Labels it using Quantum Espresso via the LabellingEngine.
    3. Trains a basic MACE model on the resulting data using the TrainingEngine.
    """
    click.echo("--- Starting MLIP-AutoPipe Cycle 01 Workflow ---")

    # Load configuration
    with open(config_path) as f:
        config = yaml.safe_load(f)

    qe_command = config['labelling']['qe_command']
    db_path = config['database']['path']
    training_config_data = config['training']

    # Initialize components
    db = AseDB(db_path)
    training_config = TrainingConfig(**training_config_data)
    labeller = LabellingEngine(qe_command=qe_command, db=db)
    trainer = TrainingEngine(config=training_config, db=db)

    # 1. Read initial structure
    click.echo(f"Reading structure from: {structure_path}")
    initial_structure = read(structure_path)

    # 2. Label the structure
    click.echo("Executing Labelling Engine...")
    db_id = labeller.execute(
        atoms=initial_structure,
        parameters=config['labelling']['parameters'],
        pseudopotentials=config['labelling']['pseudopotentials'],
        kpoints=tuple(config['labelling']['kpoints'])
    )
    click.echo(f"Labelling complete. Data saved to database with ID: {db_id}")

    # 3. Train the model
    click.echo("Executing Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    click.echo(f"Training complete. Model saved to: {trained_model_path}")

    click.echo("--- Workflow complete. ---")

if __name__ == '__main__':
    run_cycle01_workflow()
