# ruff: noqa: D101, D102, D103, D104, D105, D107
from pathlib import Path

import yaml
from ase.io import read

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine


def run_cycle01_workflow(config_path: Path, structure_path: Path, db_path: str = "mlip.db"):
    """
    Orchestrates the Cycle 01 workflow: label a structure and train a model.
    """
    # Load configuration
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Initialize database
    db = AseDB(db_path)

    # Initialize Labelling Engine
    labeller = LabellingEngine(
        qe_command=config["qe_command"],
        db=db,
        pseudo_potentials=config["pseudo_potentials"],
        k_points=tuple(config["dft_params"]["k_points"]),
        ecutwfc=config["dft_params"]["ecutwfc"],
    )

    # Initialize Training Engine
    training_config = TrainingConfig(**config["training"])
    trainer = TrainingEngine(config=training_config, db=db)

    # --- Execute Workflow ---
    print(f"Reading initial structure from: {structure_path}")
    initial_structure = read(structure_path)

    print("Executing Labelling Engine...")
    db_id = labeller.execute(initial_structure)

    # Check if labelling was successful before proceeding
    _, kvp = db.get(db_id)
    if not kvp.get("was_successful", False):
        error_msg = kvp.get("error_message", "Unknown error")
        print(f"ERROR: Labelling failed for {structure_path.name}. Reason: {error_msg}")
        return

    print(f"Labelling complete. Data stored in DB with ID: {db_id}")

    print("Executing Training Engine...")
    trained_model_path = trainer.execute(ids=[db_id])
    print(f"Workflow complete. Model saved to: {trained_model_path}")
