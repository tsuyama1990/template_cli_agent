import yaml
from ase.io import read
from .modules.c_labelling_engine import LabellingEngine
from .modules.d_training_engine import TrainingEngine
from .data.database import AseDB
from .data.models import TrainingConfig

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Executes the simple, linear workflow for Cycle 01.

    1. Reads the configuration file.
    2. Initializes the database and engine modules.
    3. Reads the initial atomic structure.
    4. Calls the LabellingEngine to perform a DFT calculation.
    5. Calls the TrainingEngine to train an MLIP on the new data.
    """
    # 1. Read configuration
    with open(config_path, 'r') as f:
        config_dict = yaml.safe_load(f)

    qe_command = config_dict["labelling_engine"]["qe_command"]
    pseudo_dir = config_dict["labelling_engine"]["pseudo_dir"]
    pseudopotentials = config_dict["labelling_engine"]["pseudopotentials"]
    training_config = TrainingConfig(**config_dict["training_engine"])

    # 2. Initialize modules
    db = AseDB("mlip.db")
    labeller = LabellingEngine(
        db=db,
        qe_command=qe_command,
        pseudo_dir=pseudo_dir,
        pseudopotentials=pseudopotentials,
    )
    trainer = TrainingEngine(config=training_config, db=db)

    # 3. Read initial structure
    initial_structure = read(structure_path)

    # 4. Run Labelling Engine
    print("Executing Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    print(f"Labelling complete. Data stored with ID: {db_id}")

    # 5. Run Training Engine
    print("Executing Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    print(f"Workflow complete. Model saved to: {trained_model_path}")
