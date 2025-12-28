# Description: Simple orchestrator for the Cycle 01 workflow.
from pathlib import Path

import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Executes the full workflow for Cycle 01.

    Args:
        config_path: Path to the YAML configuration file.
        structure_path: Path to the atomic structure file (e.g., CIF, POSCAR).
    """
    # Load configuration
    config = yaml.safe_load(Path(config_path).read_text())

    # Initialize database and engines
    db_path = "mlip.db"
    db = AseDB(db_path)
    labeller = LabellingEngine(qe_command=config["labelling_engine"]["qe_command"], db=db)

    training_config_data = config["training_engine"]
    training_config = TrainingConfig(**training_config_data)
    trainer = TrainingEngine(config=training_config, db=db)

    print(f"Processing structure from: {structure_path}")
    initial_structure = read(structure_path)

    print("Step 1: Running Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    print(f"Labelling complete. Data stored with ID: {db_id}")

    # Verify that the labelling was successful before training
    _, kvp = db.get(db_id)
    if not kvp.get("was_successful", False):
        print(
            "ERROR: Labelling failed. "
            f"Reason: {kvp.get('error_message', 'Unknown error')}. "
            "Skipping training."
        )
        return

    print("Step 2: Running Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    print(f"Workflow complete. Model saved to: {trained_model_path}")
