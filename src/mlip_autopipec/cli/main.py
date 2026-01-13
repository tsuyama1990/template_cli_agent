from pathlib import Path

import typer
import yaml
from rich.console import Console

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator

app = typer.Typer(
    name="mlip-autopipec",
    help="An automated pipeline for generating MLIP training data.",
    add_completion=False,
)

console = Console()


@app.command()
def run(
    config_path: Path = typer.Option(  # noqa: B008
        ...,
        "--config",
        "-c",
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to the YAML configuration file.",
    ),
) -> None:
    """
    Run the full data generation pipeline based on a configuration file.
    """
    if not config_path.is_file():
        console.print(f"[bold red]Error:[/bold red] Config file not found at: {config_path}")
        raise typer.Exit(code=1)

    console.print(f"🚀 Loading configuration from: [bold cyan]{config_path}[/bold cyan]")
    try:
        with config_path.open() as f:
            raw_config = yaml.safe_load(f)

        # Validate the configuration using Pydantic
        config = FullConfig.model_validate(raw_config)
        console.print("✅ Configuration loaded and validated successfully.")

        # Initialize and run the pipeline
        orchestrator = PipelineOrchestrator(config)
        orchestrator.run()

        console.print("\n🎉 [bold green]Pipeline finished successfully![/bold green]")

    except yaml.YAMLError as e:
        console.print(f"[bold red]Error parsing YAML file:[/bold red]\n{e}")
        raise typer.Exit(code=1) from e
    except Exception as e:
        console.print(f"[bold red]An unexpected error occurred:[/bold red]\n{e}")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
