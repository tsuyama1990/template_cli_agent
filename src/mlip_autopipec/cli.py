from pathlib import Path
import yaml

import typer
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

app = typer.Typer()
console = Console()


def load_config(config_path: Path) -> FullConfig:
    """Loads and validates the configuration file."""
    with open(config_path, "r") as f:
        config_dict = yaml.safe_load(f)
    return FullConfig(**config_dict)


@app.command()
def run(
    config_path: Path = typer.Option(
        ...,
        "--config",
        "-c",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
) -> None:
    """Runs the MLIP-AutoPipe data generation pipeline."""
    try:
        console.print(f"Loading configuration from: {config_path}", style="bold green")
        config = load_config(config_path)
        runner = PipelineRunner(config)
        runner.run()
    except Exception as e:
        console.print(f"An error occurred: {e}", style="bold red")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
