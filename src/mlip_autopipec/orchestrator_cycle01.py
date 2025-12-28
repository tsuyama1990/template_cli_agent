import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Orchestrates the Cycle 01 workflow: label a structure and train a model.

    Args:
        config_path: Path to the YAML configuration file.
        structure_path: Path to the atomic structure file (e.g., XYZ, CIF).
    """
    print("--- Starting MLIP-AutoPipe Cycle 01 Workflow ---")

    # 1. Load Configuration
    print(f"Loading configuration from {config_path}...")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    qe_config = config.get("quantum_espresso", {})
    train_config_dict = config.get("training", {})

    # Validate training config with Pydantic
    training_config = TrainingConfig(**train_config_dict)

    # 2. Initialize Components
    print("Initializing components (Database, Labelling Engine, Training Engine)...")
    db = AseDB("mlip.db")
    labeller = LabellingEngine(qe_command=qe_config.get("command"), db=db)
    trainer = TrainingEngine(config=training_config, db=db)

    # 3. Read Initial Structure
    print(f"Reading initial structure from {structure_path}...")
    try:
        initial_structure = read(structure_path)
    except FileNotFoundError:
        print(f"Error: Structure file not found at {structure_path}")
        return

    # 4. Run Labelling Engine
    print("Executing Labelling Engine to perform DFT calculation...")
    try:
        db_id = labeller.execute(initial_structure)
        print(f"Labelling complete. Data stored in database with ID: {db_id}")
    except Exception as e:
        print(f"An error occurred during labelling: {e}")
        return

    # 5. Run Training Engine
    print("Executing Training Engine to train the MLIP model...")
    try:
        # Verify the calculation was successful before training
        _, kvp = db.get(db_id)
        if not kvp.get("was_successful", False):
            print(f"Skipping training because labelling failed. Reason: {kvp.get('error_message')}")
            print("--- Workflow Finished (with errors) ---")
            return

        trained_model_path = trainer.execute(ids=[db_id])
        print("Training complete.")
        print(f"Workflow complete. Model saved to: {trained_model_path}")
        print("--- Workflow Finished Successfully ---")

    except Exception as e:
        print(f"An error occurred during training: {e}")
        print("--- Workflow Finished (with errors) ---")
