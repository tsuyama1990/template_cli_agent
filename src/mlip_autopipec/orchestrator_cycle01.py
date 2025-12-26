from pathlib import Path
import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.models import TrainingConfig

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Orchestrates the Cycle 01 workflow: Label -> Train.
    """
    print("--- Starting Cycle 01 Workflow ---")

    # 1. Load Configuration
    print(f"Loading configuration from {config_path}...")
    with open(config_path, 'r') as f:
        config_dict = yaml.safe_load(f)

    qe_command = config_dict["qe_command"]
    training_config = TrainingConfig(**config_dict["training"])
    db_path = config_dict.get("db_path", "mlip.db")

    # 2. Initialize Components
    db = AseDB(db_path)
    labeller = LabellingEngine(qe_command, db)
    trainer = TrainingEngine(training_config, db)

    # 3. Load Initial Structure
    print(f"Loading initial structure from {structure_path}...")
    initial_structure = read(structure_path)

    # 4. Run Labelling Engine
    print("Executing Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    label_result = db.get_row(db_id)
    if not label_result.get("was_successful", False):
        print(f"ERROR: Labelling failed. Reason: {label_result.get('error_message', 'Unknown')}")
        return

    print(f"Labelling complete. DFT results stored with ID: {db_id}")

    # 5. Run Training Engine
    print("Executing Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])

    print("\n--- Workflow complete ---")
    print(f"Trained model saved to: {trained_model_path}")
