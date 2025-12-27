import yaml
from ase.io import read
from typing import Dict, Any

from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Runs the simple, linear workflow for Cycle 01.
    """
    print("--- Starting Cycle 01 Workflow ---")

    # 1. Load configuration
    print(f"Loading configuration from {config_path}...")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # 2. Initialize components
    db_path = config.get('db_path', 'mlip.db')
    db = AseDB(db_path)

    training_config_data = config.get('training', {})
    training_config = TrainingConfig(**training_config_data)

    labeller = LabellingEngine(
        qe_command=config.get('qe_command', 'pw.x'),
        db=db,
        default_params=config.get('dft_params', {}),
        pseudos=config.get('pseudos', {})
    )
    trainer = TrainingEngine(config=training_config, db=db)

    # 3. Load initial structure
    print(f"Loading initial structure from {structure_path}...")
    initial_structure = read(structure_path)

    # 4. Run the Labelling Engine
    print("Running Labelling Engine...")
    try:
        db_id = labeller.execute(initial_structure)
        print(f"Labelling complete. Data saved to database with ID: {db_id}")
    except Exception as e:
        print(f"ERROR: An exception occurred during labelling: {e}")
        # In a more robust system, we would have more granular error handling.
        # For Cycle 1, we just print the error and exit.
        return

    # Check if labelling was successful before proceeding
    _, kvp = db.get(db_id)
    if not kvp.get('was_successful', False):
        print(f"ERROR: Labelling failed for {structure_path}. Reason: {kvp.get('error_message', 'Unknown')}")
        print("Aborting workflow.")
        return

    # 5. Run the Training Engine
    print("\nRunning Training Engine...")
    try:
        trained_model_path = trainer.execute(ids=[db_id])
        print(f"Workflow complete. Model saved to: {trained_model_path}")
    except Exception as e:
        print(f"ERROR: An exception occurred during training: {e}")

    print("--- Cycle 01 Workflow Finished ---")
