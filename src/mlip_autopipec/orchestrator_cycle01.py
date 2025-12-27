from pathlib import Path

import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def run_cycle01_workflow(config_path: Path, structure_path: Path):
    """
    Runs the simple, linear workflow for Cycle 01.

    Args:
        config_path: Path to the YAML configuration file.
        structure_path: Path to the atomic structure file (e.g., CIF, XYZ).
    """
    with open(config_path) as f:
        config_data = yaml.safe_load(f)

    db_path = Path(config_data["database"]["path"])
    db = AseDB(db_path)

    # Setup Labelling Engine
    labeller = LabellingEngine(
        qe_command=config_data["dft"]["qe_command"],
        db=db,
        parameters=config_data["dft"]["parameters"],
        pseudopotentials=config_data["dft"]["pseudopotentials"],
    )

    # Setup Training Engine
    training_config = TrainingConfig(**config_data["training"])
    trainer = TrainingEngine(config=training_config, db=db)

    # 1. Read initial structure
    print(f"Reading initial structure from: {structure_path}")
    initial_structure = read(structure_path)

    # 2. Run Labelling Engine
    print("Executing Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    print(f"Labelling complete. Data saved with ID: {db_id}")

    # Check if labelling was successful before proceeding
    _, kvp = db.get(db_id)
    if not kvp.get("was_successful", False):
        print(f"ERROR: Labelling failed. Reason: {kvp.get('error_message', 'Unknown')}")
        print("Workflow terminated.")
        return

    # 3. Run Training Engine
    print("Executing Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    print(f"Workflow complete. Model saved to: {trained_model_path}")
