import click

from mlip_autopipec.factories import create_workflow_orchestrator


@click.group()
def app():
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass


@app.command()
@click.option(
    "--id",
    required=True,
    type=int,
    help="The ID of the structure to label.",
)
@click.option(
    "--db-path",
    default="asedb.db",
    help="Path to the ASE database file.",
)
def label(id: int, db_path: str):
    """Runs the DFT labeling for a single structure."""
    orchestrator = create_workflow_orchestrator(db_path)
    orchestrator.label_structure_by_id(id)


@app.command()
@click.option(
    "--db-path",
    default="asedb.db",
    help="Path to the ASE database file.",
)
def train(db_path: str):
    """Trains the MLIP model on the existing labeled data."""
    orchestrator = create_workflow_orchestrator(db_path)
    orchestrator.run_training()


if __name__ == "__main__":
    app()
