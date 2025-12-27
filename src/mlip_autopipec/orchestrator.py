"""
Orchestrator for the Cycle 01 workflow.
"""
import yaml
from pathlib import Path
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.data.models import TrainingConfig

def run_workflow(config_path: Path, structure_path: Path):
    """
    Executes the full Cycle 01 workflow: Label -> Train.

    Args:
        config_path: Path to the YAML configuration file.
        structure_path: Path to the input atomic structure file (e.g., CIF, POSCAR).
    """
    # 1. Load configuration
    print(f"Loading configuration from {config_path}...")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # 2. Initialize components
    db_path = Path("mlip_pipe.db")
    db = AseDB(db_path)
    print(f"Database initialized at {db_path.resolve()}")

    labelling_engine = LabellingEngine(
        db=db,
        qe_command=config["qe_command"],
        pseudo_dir=config["pseudo_dir"],
        pseudos=config["pseudos"],
        kpts=tuple(config["kpts"]),
        ecutwfc=config["ecutwfc"],
    )

    training_config = TrainingConfig(**config["training"])
    training_engine = TrainingEngine(
        config=training_config,
        db=db,
    )

    # 3. Load initial structure
    print(f"Loading initial structure from {structure_path}...")
    initial_structure = read(structure_path)

    # 4. Execute Labelling Engine
    print("Executing Labelling Engine...")
    db_id = labelling_engine.execute(initial_structure)

    # Check if labelling was successful before proceeding
    _, result = db.get(db_id)
    if not result.was_successful:
        print(f"Labelling failed: {result.error_message}")
        print("Workflow terminated.")
        return

    print(f"Labelling complete. Result stored in DB with ID: {db_id}")

    # 5. Execute Training Engine
    print("Executing Training Engine...")
    model_path = training_engine.execute(ids=[db_id])

    print("\n--- Workflow Complete ---")
    print(f"Trained model saved to: {Path(model_path).resolve()}")
    print("-------------------------")
