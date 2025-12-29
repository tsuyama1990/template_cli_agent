import os
from typing import Any

import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def load_config(config_path: str) -> dict[str, Any]:
    """Loads the YAML configuration file."""
    with open(config_path) as f:
        return yaml.safe_load(f)

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Orchestrates the Cycle 01 workflow: Label -> Train.

    Args:
        config_path: Path to the YAML configuration file.
        structure_path: Path to the initial atomic structure file (e.g., XYZ, CIF).
    """
    print("--- Starting Cycle 01 Workflow ---")

    # 1. Load configuration
    print(f"Loading configuration from {config_path}...")
    config = load_config(config_path)

    # Get the directory of the config file to resolve relative paths
    config_dir = os.path.dirname(config_path)

    # 2. Initialize components
    db_path = os.path.join(config_dir, config['db_path'])
    print(f"Initializing database at {db_path}...")
    db = AseDB(db_path)

    print("Initializing Labelling Engine...")
    labeller = LabellingEngine(
        qe_command=config['qe_command'],
        parameters=config['dft_parameters'],
        db=db
    )

    print("Initializing Training Engine...")
    training_config_model = TrainingConfig(**config['training'])
    trainer = TrainingEngine(config=training_config_model, db=db)

    # 3. Read initial structure
    print(f"Reading initial structure from {structure_path}...")
    initial_structure = read(structure_path)

    # 4. Execute Labelling
    print("Executing Labelling Engine...")
    try:
        db_id = labeller.execute(initial_structure)
        print(f"Labelling complete. Result saved to database with ID: {db_id}")
    except Exception as e:
        print(f"An error occurred during labelling: {e}")
        return

    # 5. Execute Training
    print("Executing Training Engine...")
    # FIX: Pass an absolute path for the output directory
    model_output_dir = os.path.join(config_dir, "models")
    try:
        trained_model_path = trainer.execute(ids=[db_id], output_dir=model_output_dir)
        print(f"Workflow complete. Model saved to: {trained_model_path}")
    except ValueError as e:
        print(f"Skipping training due to an issue with the data: {e}")
    except Exception as e:
        print(f"An error occurred during training: {e}")

    print("--- Cycle 01 Workflow Finished ---")
