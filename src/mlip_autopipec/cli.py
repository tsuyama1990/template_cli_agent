import logging
from pathlib import Path
from typing import Annotated

import typer
import yaml
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer()
console = Console()
logging.basicConfig(level=logging.INFO)


@app.command()
def run(
    config_path: Annotated[
        Path,
        typer.Option("--config", "-c", help="Path to the YAML configuration file."),
    ],
) -> None:
    """Run the MLIP-AutoPipe pipeline."""
    if not config_path.exists():
        console.print(f"[bold red]Error: Config file not found at '{config_path}'[/bold red]")
        raise typer.Exit(code=1)

    try:
        with config_path.open("r") as f:
            raw_config = yaml.safe_load(f)
        config = FullConfig.model_validate(raw_config)

        console.print(
            f"Successfully loaded configuration from [bold green]{config_path}[/bold green]"
        )
        runner = PipelineRunner(config)
        runner.run()
        console.print("[bold green]Pipeline finished successfully.[/bold green]")

    except Exception as e:
        console.print(f"[bold red]An error occurred: {e}[/bold red]")
        raise typer.Exit(code=1) from e
