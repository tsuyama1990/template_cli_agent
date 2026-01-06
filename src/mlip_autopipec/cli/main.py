"""
Main Command-Line Interface (CLI) for MLIP-AutoPipe.
"""

import typer
import yaml
import logging
from pathlib import Path
from rich.console import Console
from rich.traceback import install
from typing_extensions import Annotated

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator
from mlip_autopipec.storage.database_manager import DatabaseManager
from mlip_autopipec.explorers.md_engine import MDEngine
from mlip_autopipec.factories import create_generator, create_sampler
from mlip_autopipec.common.exceptions import PipelineError, ConfigurationError

from ase.calculators.emt import EMT

# Set up rich console for pretty printing and exception handling
console = Console()
install(show_locals=False) # Rich exception handler

app = typer.Typer(
    name="mlip-autopipec",
    help="A tool to automate the generation of training data for MLIPs.",
    add_completion=False,
)

# Configure top-level logger
logging.basicConfig(
    level="INFO",
    format="%(asctime)s - [%(levelname)s] - %(name)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def load_config(config_path: Path) -> FullConfig:
    """Loads and validates the YAML configuration file."""
    logger.info(f"Loading configuration from: {config_path}")
    if not config_path.exists():
        msg = f"Configuration file not found at: {config_path}"
        raise ConfigurationError(msg)
    try:
        with config_path.open("r") as f:
            config_dict = yaml.safe_load(f)
        return FullConfig(**config_dict)
    except yaml.YAMLError as e:
        msg = f"Error parsing YAML file: {e}"
        raise ConfigurationError(msg) from e
    except Exception as e: # Catches Pydantic's ValidationError
        msg = f"Configuration validation error: {e}"
        raise ConfigurationError(msg) from e


@app.command()
def run_pipeline(
    config_file_str: Annotated[
        str,
        typer.Argument(
            help="Path to the YAML configuration file.",
        ),
    ],
    db_path_str: Annotated[
        str,
        typer.Option(
            "--db-path", "-db",
            help="Path to the output ASE database file. [default: output.db]",
        ),
    ] = "output.db",
    debug: Annotated[
        bool,
        typer.Option("--debug", help="Enable debug logging."),
    ] = False,
) -> None:
    """
    Runs the full MLIP-AutoPipe data generation pipeline.
    """
    config_file = Path(config_file_str)
    if not config_file.is_file():
        console.print(f"[bold red]Error:[/bold red] Config file not found at: {config_file}")
        raise typer.Exit(code=1)

    db_path = Path(db_path_str)

    if debug:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.debug("Debug mode enabled.")

    console.rule("[bold green]MLIP-AutoPipe: Data Generation Pipeline[/bold green]")

    try:
        # 1. Load and validate configuration
        config = load_config(config_file)
        console.log("✅ Configuration loaded and validated.")

        # 2. Instantiate components
        calculator = EMT() # type: ignore
        console.log("🔧 Using placeholder EMT calculator for this run.")

        generator = create_generator(config)
        sampler = create_sampler(config)
        explorer = MDEngine(config.exploration, calculator)
        storage = DatabaseManager(db_path)
        console.log("🔧 Core components instantiated.")

        # 3. Instantiate and run the orchestrator
        orchestrator = PipelineOrchestrator(
            config=config,
            generator=generator,
            explorer=explorer,
            sampler=sampler,
            storage=storage,
        )

        orchestrator.run_pipeline()
        console.rule("[bold green]Pipeline Complete[/bold green]")
        console.log(f"✅ Successfully generated dataset at: {db_path}")

    except PipelineError as e:
        console.print(f"[bold red]Pipeline Error:[/bold red] {e}", style="red")
        logger.exception("Pipeline failed due to a controlled error:")
        raise typer.Exit(code=1) from e
    except Exception as e:
        console.print("[bold red]An unexpected critical error occurred.[/bold red]")
        # Rich will print the traceback automatically
        logger.critical("Pipeline failed due to an unexpected error.", exc_info=True)
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
