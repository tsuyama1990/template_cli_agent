# -*- coding: utf-8 -*-
"""
Typer-based Command-Line Interface application for MLIP-AutoPipe.

This module provides the main entry point for the user to interact with the
data generation pipeline.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
import yaml
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer()
console = Console()


def load_config(config_path: Path) -> FullConfig:
    """
    Load and validate the configuration file from the given path.

    Args:
        config_path: The path to the YAML configuration file.

    Returns:
        A validated FullConfig object.

    Raises:
        typer.Exit: If the configuration file is not found.
    """
    if not config_path.exists():
        console.print(f"[bold red]Error: Configuration file not found at {config_path}[/bold red]")
        raise typer.Exit(code=1)
    with config_path.open() as f:
        config_dict = yaml.safe_load(f)
    return FullConfig(**config_dict)


@app.command()
def run(
    config_path_str: Optional[str] = typer.Option(  # noqa: B008
        None, "--config", "-c", help="Path to the YAML configuration file."
    )
) -> None:
    """
    Run the MLIP-AutoPipe data generation pipeline.

    This command initializes and executes the full workflow, including structure
    generation, exploration, sampling, and storage, based on the provided
    YAML configuration file.
    """
    if config_path_str is None:
        console.print("[bold red]Error: --config option is required.[/bold red]")
        raise typer.Exit(code=1)

    config_path = Path(config_path_str)

    console.print("[bold green]Starting MLIP-AutoPipe...[/bold green]")
    try:
        config = load_config(config_path)
        runner = PipelineRunner(config)
        runner.run()
        console.print("[bold green]MLIP-AutoPipe finished successfully.[/bold green]")
    except Exception as e:
        console.print(f"[bold red]An error occurred: {e}[/bold red]")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
