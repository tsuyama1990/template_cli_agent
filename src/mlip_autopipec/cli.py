"""Command-line interface (CLI) for the MLIP-AutoPipe application.

This module provides the main entry point for users to interact with the
MLIP-AutoPipe tool. It is built using the Typer library, which provides a clean,
modern interface with automatic help generation and robust argument parsing.
The CLI is responsible for parsing the user's configuration file, validating it,
and initiating the main data generation pipeline.
"""

from pathlib import Path

import typer
import yaml
from pydantic import ValidationError
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer(
    name="mlip-autopipec",
    help="An automated pipeline for generating MLIP training data.",
    add_completion=False,
)

console = Console()


def load_config(config_path: Path) -> FullConfig | None:
    """Load, parse, and validate the YAML configuration file.

    This function serves as the primary input boundary for the application. It
    handles file reading, YAML parsing, and, most importantly, validation
    against the Pydantic `FullConfig` schema.

    Args:
        config_path: The path to the user-provided YAML configuration file.

    Returns:
        A validated `FullConfig` object if the file is valid, otherwise `None`.

    """
    if not config_path.exists():
        console.print(f"[bold red]Error:[/] Config file not found at: {config_path}")
        return None

    try:
        with config_path.open("r") as f:
            config_dict = yaml.safe_load(f)
        return FullConfig(**config_dict)
    except ValidationError as e:
        console.print("[bold red]Error:[/] Configuration validation failed.")
        console.print(e)
        return None
    except yaml.YAMLError as e:
        console.print(f"[bold red]Error:[/] Failed to parse YAML file: {e}")
        return None


@app.command()
def run(
    config: Path = typer.Option(
        ...,
        "--config",
        "-c",
        help="Path to the YAML configuration file.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
) -> None:
    """Run the full data generation pipeline.

    This command orchestrates the entire MLIP-AutoPipe workflow, from loading
    the configuration to generating structures and saving them to a database.
    """
    console.print(f"Loading configuration from: {config}")
    full_config = load_config(config)

    if full_config is None:
        raise typer.Exit(code=1)

    try:
        runner = PipelineRunner(config=full_config)
        runner.run()
    except Exception as e:
        console.print(
            f"[bold red]An unexpected error occurred during the pipeline execution:[/] {e}"
        )
        raise typer.Exit(code=1)

    console.print("[bold green]MLIP-AutoPipe finished successfully![/]")


if __name__ == "__main__":
    app()
