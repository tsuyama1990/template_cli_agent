import click

from mlip_autopipec.workflow import WorkflowOrchestrator


@click.command()
@click.option(
    "--config-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the main YAML configuration file.",
)
@click.option(
    "--database-file",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to the ASE database file to use or create.",
)
@click.option(
    "--input-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the XYZ file containing initial atomic structures.",
)
def mlip_pipe(config_file: str, database_file: str, input_file: str):
    """
    Main entrypoint for the MLIP-AutoPipe workflow.
    """
    click.echo("Workflow initialized via CLI.")

    try:
        orchestrator = WorkflowOrchestrator(
            config_path=config_file,
            db_path=database_file,
            input_file_path=input_file,
        )
        orchestrator.run_workflow()
        click.secho("✅ Workflow completed successfully!", fg="green")
    except Exception as e:
        click.secho(f"❌ Workflow failed: {e}", fg="red")
        # In a real app, you might want more specific error handling
        # and possibly a traceback in debug mode.
        raise click.Abort()

if __name__ == "__main__":
    mlip_pipe()
