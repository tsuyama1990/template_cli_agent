"""Typer-based Command-Line Interface application."""
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
    """Run the MLIP-AutoPipe pipeline."""
    if not config_path.exists():
        console.print(f"[bold red]Error:[/] The specified configuration file '{config_path}' was not found.")
        raise typer.Exit(code=1)

    with open(config_path) as f:
        config_dict = yaml.safe_load(f)

    try:
        config = FullConfig(**config_dict)
    except Exception as e:
        console.print(f"[bold red]Error:[/] Configuration validation failed. {e}")
        raise typer.Exit(code=1)

    runner = PipelineRunner(config)
    runner.run()


if __name__ == "__main__":
    app()
