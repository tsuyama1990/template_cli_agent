from pathlib import Path

import ase.io
import click

from mlip_autopipec.configs.config_loader import load_config
from mlip_autopipec.modules.labelling_engine import LabellingEngine
from mlip_autopipec.modules.training_engine import TrainingEngine
from mlip_autopipec.utils.ase_db import AseDB


@click.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the main YAML configuration file.",
)
@click.option(
    "--structure",
    "structure_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the atomic structure file (e.g., POSCAR, CIF).",
)
def main(config_path: Path, structure_path: Path):
    """
    Main orchestrator for the Cycle 01 workflow.
    """
    # Load main configuration
    config = load_config(config_path)

    # Read atomic structure
    atoms = ase.io.read(structure_path)

    # --- Labelling Step ---
    labelling_engine = LabellingEngine()
    dft_result, was_successful = labelling_engine.execute(atoms)

    # --- Database Step ---
    db_path = Path("mlip.db")
    ase_db = AseDB(db_path)
    uid = ase_db.write(
        atoms, dft_result, was_successful, key_value_pairs={"source": str(structure_path)}
    )
    click.echo(f"Wrote structure {uid} to database at {db_path}")

    if not was_successful:
        click.echo("Labelling failed. Skipping training.")
        return

    # --- Training Step ---
    training_engine = TrainingEngine()
    model_path = training_engine.execute(config.training, db_path)

    click.echo(f"Workflow complete. Model saved to: {model_path}")


if __name__ == "__main__":
    main()
