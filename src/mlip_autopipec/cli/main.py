from pathlib import Path

import typer
import yaml
from rich.console import Console

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator

app = typer.Typer()
console = Console()

@app.command()
def run(
    config_path: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Path to the YAML configuration file.",
    ),
    db_path: Path = typer.Option(
        "structures.db",
        "--db-path",
        "-db",
        help="Path to the output ASE database file.",
    ),
):
    """
    Runs the full MLIP-AutoPipe data generation pipeline.
    """
    console.log(f"Loading configuration from: [cyan]{config_path}[/cyan]")
    with open(config_path) as f:
        config_dict = yaml.safe_load(f)

    try:
        config = FullConfig(**config_dict)
    except Exception as e:
        console.log("[bold red]Configuration Error:[/bold red]")
        console.print(e)
        raise typer.Exit(code=1)

    console.log("Configuration loaded successfully.")

    orchestrator = PipelineOrchestrator(config, db_path)
    orchestrator.run_pipeline()


if __name__ == "__main__":
    app()
