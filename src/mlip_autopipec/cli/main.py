from pathlib import Path

import typer
import yaml
from rich.console import Console

from mlip_autopipec import factories
from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator
from mlip_autopipec.storage.database_manager import DatabaseManager

app = typer.Typer()
console = Console()


@app.command()
def run(
    config_path: Path = typer.Option(
        ..., "--config", "-c", help="Path to the YAML configuration file."
    )
):
    """
    Runs the MLIP-AutoPipe data generation pipeline.
    """
    try:
        with open(config_path) as f:
            config_dict = yaml.safe_load(f)
        config = FullConfig(**config_dict)

        output_dir = Path.cwd()
        db_path = output_dir / "results.db"

        # Create components using factories
        generator = factories.create_generator(config)
        explorer = factories.create_explorer(config)
        sampler = factories.create_sampler(config)
        db_manager = DatabaseManager(db_path=str(db_path))

        # Instantiate and run the orchestrator
        orchestrator = PipelineOrchestrator(
            config=config,
            generator=generator,
            explorer=explorer,
            sampler=sampler,
            db_manager=db_manager,
            output_dir=output_dir,
        )
        orchestrator.run()

        console.print("[bold green]Pipeline completed successfully![/bold green]")

    except Exception as e:
        console.print(f"[bold red]An error occurred: {e}[/bold red]")
        raise typer.Exit(code=1)
