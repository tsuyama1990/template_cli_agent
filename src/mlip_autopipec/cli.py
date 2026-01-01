import click
from pydantic import ValidationError

from mlip_autopipec.config import Settings
from mlip_autopipec.factories import create_workflow_orchestrator


@click.group()
@click.pass_context
def app(ctx):
    """MLIP-AutoPipe: Automated MLIP Generation."""
    try:
        settings = Settings()
        ctx.obj = settings
    except ValidationError as e:
        raise click.ClickException(f"Configuration error: {e}") from e


@app.command()
@click.option(
    "--id",
    required=True,
    type=int,
    help="The ID of the structure to label.",
)
@click.option(
    "--db-path",
    help="Path to the ASE database file.",
)
@click.pass_context
def label(ctx, id: int, db_path: str | None):
    """Runs the DFT labeling for a single structure."""
    settings = ctx.obj
    db_path = db_path or settings.db_path
    orchestrator = create_workflow_orchestrator(db_path, settings)
    orchestrator.label_structure_by_id(id)


@app.command()
@click.option(
    "--db-path",
    help="Path to the ASE database file.",
)
@click.pass_context
def train(ctx, db_path: str | None):
    """Trains the MLIP model on the existing labeled data."""
    settings = ctx.obj
    db_path = db_path or settings.db_path
    orchestrator = create_workflow_orchestrator(db_path, settings)
    orchestrator.run_training()


if __name__ == "__main__":
    app()
