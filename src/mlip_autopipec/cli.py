"""Command-line interface for MLIP-AutoPipe."""
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
    """
    Loads, parses, and validates the YAML configuration file.

    Args:
        config_path: Path to the configuration file.

    Returns:
        A validated FullConfig object, or None if validation fails.
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
    """
    Run the full data generation pipeline.
    """
    console.print(f"Loading configuration from: {config}")
    full_config = load_config(config)

    if full_config is None:
        raise typer.Exit(code=1)

    try:
        runner = PipelineRunner(config=full_config)
        runner.run()
    except Exception as e:
        console.print(f"[bold red]An unexpected error occurred during the pipeline execution:[/] {e}")
        raise typer.Exit(code=1)

    console.print("[bold green]MLIP-AutoPipe finished successfully![/]")


if __name__ == "__main__":
    app()
