"""Main CLI entry point for the MLIP-AutoPipe Cycle 1 workflow."""
import click
import yaml
from ase.io import read
from .modules.c_labelling_engine import LabellingEngine
from .modules.d_training_engine import TrainingEngine
from .data.database import AseDB
from .data.models import TrainingConfig

@click.command()
@click.option('--config', 'config_path', required=True, type=click.Path(exists=True), help='Path to the YAML configuration file.')
@click.option('--structure', 'structure_path', required=True, type=click.Path(exists=True), help='Path to the atomic structure file (e.g., POSCAR, cif, xyz).')
def main(config_path, structure_path):
    """
    Executes the Cycle 01 workflow: DFT labelling of a single structure
    followed by training a rudimentary MLIP model.
    """
    click.echo("--- Starting MLIP-AutoPipe Cycle 01 Workflow ---")

    # 1. Load configuration
    click.echo(f"Loading configuration from: {config_path}")
    with open(config_path, 'r') as f:
        config_data = yaml.safe_load(f)

    # 2. Initialize components
    db_path = config_data.get("database_path", "mlip.db")
    db = AseDB(db_path)

    training_config = TrainingConfig(**config_data["training"])

    labeller = LabellingEngine(
        qe_command=config_data["qe_command"],
        parameters=config_data["dft_parameters"],
        pseudos=config_data["pseudopotentials"],
        db=db
    )

    trainer = TrainingEngine(
        config=training_config,
        db=db
    )

    # 3. Load initial structure
    click.echo(f"Loading initial structure from: {structure_path}")
    initial_structure = read(structure_path)

    # 4. Execute Labelling Engine
    click.echo("Executing Labelling Engine...")
    try:
        db_id = labeller.execute(initial_structure)
        click.echo(f"  -> Labelling complete. Data stored in DB with ID: {db_id}")
    except Exception as e:
        click.echo(f"  -> ERROR during labelling: {e}", err=True)
        return

    # 5. Execute Training Engine
    click.echo("Executing Training Engine...")
    try:
        model_path = trainer.execute(ids=[db_id])
        click.echo(f"  -> Training complete. Model saved to: {model_path}")
    except Exception as e:
        click.echo(f"  -> ERROR during training: {e}", err=True)
        return

    click.echo("--- Workflow Finished Successfully ---")

if __name__ == '__main__':
    main()
