# src/mlip_autopipec/cli/main.py
"""The main Command-Line Interface for MLIP-AutoPipe."""

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


@app.command()
def run(
    config_path: str = typer.Argument(
        ...,
        help="Path to the YAML configuration file.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
) -> None:
    """
    Runs the full MLIP-AutoPipe data generation pipeline.
    """
    console.print("🚀 [bold green]Starting MLIP-AutoPipe pipeline...[/bold green]")
    console.print(f"   - Loading configuration from: {config_path}")

    try:
        with open(config_path) as f:
            config_data = yaml.safe_load(f)

        config = FullConfig(**config_data)

        orchestrator = PipelineOrchestrator(config)
        orchestrator.run()

    except Exception as e:
        console.print("🔥 [bold red]An error occurred:[/bold red]")
        console.print(e)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
