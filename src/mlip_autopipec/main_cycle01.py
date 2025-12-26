import argparse
import yaml
from ase.io import read
from pathlib import Path

from mlip_autopipec.data.database import AseDB
from mlip_autopipec.data.models import TrainingConfig
from mlip_autopipec.modules.c_labelling_engine import LabellingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine

def run_workflow(config: dict, structure_path: str):
    """
    Orchestrates the Cycle 01 workflow: Label -> Train.

    Args:
        config: A dictionary containing the workflow configuration.
        structure_path: Path to the initial atomic structure file (e.g., xyz, cif).
    """
    print("--- Starting MLIP-AutoPipe Cycle 01 Workflow ---")

    # 1. Initialize components from config
    print(f"Loading configuration and initializing components...")
    db_path = config["database"]["path"]
    db = AseDB(db_path)

    dft_config = config["dft"]
    labeller = LabellingEngine(
        qe_command=dft_config["qe_command"],
        db=db,
        pseudos=dft_config["pseudos"],
        kpts=tuple(dft_config["kpts"]),
        ecutwfc=dft_config["ecutwfc"],
    )

    training_cfg_model = TrainingConfig(**config["training"])
    trainer = TrainingEngine(config=training_cfg_model, db=db)

    # 2. Load initial structure
    print(f"Reading initial structure from: {structure_path}")
    try:
        initial_structure = read(structure_path)
    except FileNotFoundError:
        print(f"Error: Structure file not found at {structure_path}")
        return

    # 3. Run the Labelling Engine
    print("Executing Labelling Engine (DFT calculation)...")
    # Note: This will execute a real DFT calculation if pw.x is in the PATH.
    # For testing, this should be run in an environment where pw.x is a mock script.
    db_id = labeller.execute(initial_structure)

    # Check if labelling was successful
    result_row = db.get(db_id)
    if not result_row.key_value_pairs["was_successful"]:
        print("Labelling failed. See database for error details.")
        print("--- Workflow Terminated ---")
        return

    print(f"Labelling complete. Result saved to database with ID: {db_id}")

    # 4. Run the Training Engine
    print("Executing Training Engine (MLIP training)...")
    trained_model_path = trainer.execute(ids=[db_id])

    print(f"\nWorkflow complete. Final model saved to: {trained_model_path}")
    print("--- MLIP-AutoPipe Cycle 01 Workflow Finished ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the MLIP-AutoPipe Cycle 01 workflow."
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the YAML configuration file.",
    )
    parser.add_argument(
        "--structure",
        type=str,
        required=True,
        help="Path to the initial atomic structure file (e.g., .cif, .xyz).",
    )
    args = parser.parse_args()

    # Load config file
    config_path = Path(args.config)
    if not config_path.is_file():
        print(f"Error: Config file not found at {config_path}")
    else:
        with open(config_path, "r") as f:
            loaded_config = yaml.safe_load(f)
        run_workflow(loaded_config, args.structure)
