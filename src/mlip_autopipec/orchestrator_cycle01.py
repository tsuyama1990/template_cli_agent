import yaml
from pathlib import Path

from ase.io import read

from .modules.c_labelling_engine import LabellingEngine
from .modules.d_training_engine import TrainingEngine
from .data.database import AseDB
from .data.models import TrainingConfig

def run_cycle01_workflow(config_path: str, structure_path: str):
    """
    Orchestrates the simple, linear workflow for Cycle 01.
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    qe_command = config.get("qe_command", "pw.x")
    training_config = TrainingConfig(**config.get("training", {}))
    db_path = Path(config.get("db_path", "./mlip.db"))

    db = AseDB(db_path)

    labeller = LabellingEngine(qe_command=qe_command, db=db)
    trainer = TrainingEngine(config=training_config, db=db)

    print(f"Reading initial structure from: {structure_path}")
    initial_structure = read(structure_path)

    print("Executing Labelling Engine...")
    db_id = labeller.execute(initial_structure)
    print(f"Labelling complete. Result stored in database with ID: {db_id}")

    _, dft_result = db.get(db_id)

    if dft_result and dft_result.was_successful:
        print("Executing Training Engine...")
        trained_model_path = trainer.execute(ids=[db_id])
        if trained_model_path:
            print(f"Workflow complete. Model saved to: {trained_model_path}")
        else:
            print("Training Engine failed to produce a model.")
    else:
        print("Labelling was not successful. Skipping training.")
