from pathlib import Path

import typer
import yaml

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator

app = typer.Typer()


@app.command()
def run(config_path: str = typer.Argument(..., help="Path to the YAML configuration file.")):
    """
    Runs the MLIP AutoPipe pipeline with the given configuration file.
    """
    try:
        config_file = Path(config_path)
        with config_file.open("r") as f:
            config_dict = yaml.safe_load(f)

        config = FullConfig(**config_dict)

        orchestrator = PipelineOrchestrator(config)
        orchestrator.run()

    except FileNotFoundError as e:
        typer.echo(f"Error: Configuration file not found at {config_path}", err=True)
        raise typer.Exit(code=1) from e
    except Exception as e:
        typer.echo(f"An error occurred: {e}", err=True)
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
