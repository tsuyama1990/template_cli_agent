"""Command-Line Interface for MLIP-AutoPipe."""
from pathlib import Path
from typing import Annotated

import typer
import yaml
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer()
console = Console()


@app.command()
def run(
    config: Annotated[
        Path,
        typer.Option(
            "--config",
            "-c",
            help="Path to the configuration YAML file.",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ]
) -> None:
    """Run the MLIP-AutoPipe data generation pipeline."""
    console.print("[bold green]Starting MLIP-AutoPipe...[/bold green]")

    try:
        with open(config) as f:
            config_dict = yaml.safe_load(f)

        console.print("Configuration loaded successfully. Validating...")
        full_config = FullConfig(**config_dict)
        console.print("[green]Configuration validation successful.[/green]")

    except (yaml.YAMLError, ValueError) as e:
        console.print(f"[bold red]Error:[/bold red] Failed to load or validate configuration: {e}")
        raise typer.Exit(code=1)

    runner = PipelineRunner(full_config)
    runner.run()

    console.print("[bold green]Pipeline finished successfully.[/bold green]")
