import yaml
from ase.io import read
from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.models import TrainingConfig

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Orchestrates the Cycle 01 workflow:
    1. Reads configuration and initial structure.
    2. Runs the Labelling Engine to perform DFT calculation.
    3. Runs the Training Engine to train an MLIP.
    """
    # Load configuration
    with open(config_path, 'r') as f:
        config_data = yaml.safe_load(f)

    # Initialize database
    db = AseDB(config_data['database_path'])

    # Initialize Labelling Engine
    labeller = LabellingEngine(
        qe_command=config_data['qe_command'],
        parameters=config_data['dft_parameters'],
        db=db
    )

    # Read initial structure
    initial_structure = read(structure_path)

    # Run Labelling
    print("Running Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    print(f"Labelling complete. Data stored with ID: {db_id}")

    # Check if labelling was successful
    db_entry = db.read(db_id)
    if not db_entry.get('was_successful', False):
        print(f"Labelling failed: {db_entry.get('error_message')}")
        print("Skipping training.")
        return

    # Initialize Training Engine
    training_config = TrainingConfig(**config_data['training'])
    trainer = TrainingEngine(config=training_config, db=db)

    # Run Training
    print("Running Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    print(f"Workflow complete. Model saved to: {trained_model_path}")
