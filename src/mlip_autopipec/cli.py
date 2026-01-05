"""Command-line interface for MLIP-AutoPipe."""
from pathlib import Path

import typer
import yaml
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer()
console = Console()


@app.command()
def run(config_path: Path = typer.Option(..., "--config", help="Path to the YAML configuration file.")):
    """Runs the MLIP-AutoPipe data generation pipeline."""
    if not config_path.exists():
        console.print(f"[bold red]Error:[/] The specified configuration file '{config_path}' was not found.")
        raise typer.Exit(code=1)

    try:
        with config_path.open("r") as f:
            config_data = yaml.safe_load(f)
        config = FullConfig(**config_data)
        runner = PipelineRunner(config)
        runner.run()
    except Exception as e:
        console.print(f"[bold red]An error occurred:[/] {e}")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
