# src/mlip_autopipec/cli/main.py
"""
The main Command-Line Interface (CLI) for the MLIP-AutoPipe application.
"""

from pathlib import Path

import typer
import yaml
from rich.console import Console

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator

app = typer.Typer(
    name="mlip-autopipec",
    help="A tool to automate the generation of training data for MLIPs.",
    add_completion=False,
)

console = Console()


def load_config(config_path: Path) -> FullConfig | None:
    """
    Loads and validates the YAML configuration file against the Pydantic model.
    """
    try:
        with config_path.open("r") as f:
            raw_config = yaml.safe_load(f)
        return FullConfig.model_validate(raw_config)
    except FileNotFoundError:
        console.print(
            f"[bold red]Error: Configuration file not found at '{config_path}'[/bold red]"
        )
        return None
    except Exception as e:
        console.print(f"[bold red]Error parsing or validating configuration file:\n{e}[/bold red]")
        return None


@app.command()
def run_pipeline(
    config_file: str = typer.Argument(
        ...,
        help="Path to the YAML configuration file.",
    ),
) -> None:
    """
    Runs the full data generation pipeline based on a configuration file.
    """
    config_path = Path(config_file)
    if not config_path.exists():
        console.print(
            f"[bold red]Error: Configuration file not found at '{config_path}'[/bold red]"
        )
        raise typer.Exit(code=1)

    console.print(f"Loading configuration from: {config_path}")
    config = load_config(config_path)

    if config is None:
        raise typer.Exit(code=1)

    orchestrator = PipelineOrchestrator(config)
    orchestrator.run()


if __name__ == "__main__":
    app()
