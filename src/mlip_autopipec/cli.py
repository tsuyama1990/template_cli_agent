"""Command-Line Interface for MLIP-AutoPipe."""
import typer
import yaml
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner
from pydantic import ValidationError

app = typer.Typer(pretty_exceptions_show_locals=False)
console = Console()


def load_config(config_path: str) -> FullConfig:
    """Loads, parses, and validates the YAML configuration file."""
    try:
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return FullConfig(**config_dict)
    except FileNotFoundError:
        console.print(f"[bold red]Error: Configuration file not found at '{config_path}'[/bold red]")
        raise typer.Exit(code=1)
    except ValidationError as e:
        console.print("[bold red]Error: Configuration validation failed.[/bold red]")
        # Print a simplified, user-friendly version of the Pydantic error.
        for error in e.errors():
            field = ".".join(map(str, error['loc']))
            message = error['msg']
            console.print(f"- [yellow]Field '{field}': {message}[/yellow]")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]An unexpected error occurred: {e}[/bold red]")
        raise typer.Exit(code=1)


@app.command()
def run(config: str = typer.Option(..., "--config", "-c", help="Path to the YAML configuration file.")):
    """
    Runs the MLIP-AutoPipe data generation pipeline.
    """
    console.print("[bold green]Starting MLIP-AutoPipe...[/bold green]")

    # Load and validate the configuration
    full_config = load_config(config)

    # Instantiate and run the pipeline
    runner = PipelineRunner(config=full_config)
    runner.run()

    console.print("[bold green]...MLIP-AutoPipe finished successfully.[/bold green]")


if __name__ == "__main__":
    app()
