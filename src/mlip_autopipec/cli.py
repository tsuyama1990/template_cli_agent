import click

from mlip_autopipec.app import create_app


@click.group()
@click.pass_context
def app(ctx):
    """MLIP-AutoPipe: Automated MLIP Generation."""
    ctx.obj = create_app()


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
    application = ctx.obj
    application.label_structure(id, db_path)


@app.command()
@click.option(
    "--db-path",
    help="Path to the ASE database file.",
)
@click.pass_context
def train(ctx, db_path: str | None):
    """Trains the MLIP model on the existing labeled data."""
    application = ctx.obj
    application.train(db_path)


if __name__ == "__main__":
    app()
