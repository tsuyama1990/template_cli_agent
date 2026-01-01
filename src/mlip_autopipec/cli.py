import click

from mlip_autopipec.config import DFTInputConfig, MLIPTrainingConfig, ModelType
from mlip_autopipec.workflow import WorkflowOrchestrator

# Dummy configurations for CYCLE01
# In a future cycle, these will be loaded from a user-provided YAML file.
DUMMY_DFT_CONFIG = DFTInputConfig(
    pseudopotentials={"H": "H.UPF", "O": "O.UPF"},
    kpoints=(1, 1, 1),
    ecutwfc=60,
    control={"calculation": "scf"},
)

DUMMY_TRAINING_CONFIG = MLIPTrainingConfig(
    model_type=ModelType.ACE,
    r_cut=5.0,
    loss_weights={"energy": 1.0, "forces": 100.0},
)

# This would also come from a config file
QE_COMMAND = "pw.x"
DB_PATH = "asedb.db"


@click.group()
def app():
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass


@app.command()
@click.option("--id", required=True, type=int, help="The ID of the structure to label.")
def label(id: int):
    """Runs the DFT labeling for a single structure."""
    orchestrator = WorkflowOrchestrator(
        dft_config=DUMMY_DFT_CONFIG,
        training_config=DUMMY_TRAINING_CONFIG,
        db_path=DB_PATH,
        qe_command=QE_COMMAND,
    )
    orchestrator.run_labeling_for_id(id)


@app.command()
def train():
    """Trains the MLIP model on the existing labeled data."""
    orchestrator = WorkflowOrchestrator(
        dft_config=DUMMY_DFT_CONFIG,
        training_config=DUMMY_TRAINING_CONFIG,
        db_path=DB_PATH,
        qe_command=QE_COMMAND,
    )
    orchestrator.run_training()


if __name__ == "__main__":
    app()
