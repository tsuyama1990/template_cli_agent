from pathlib import Path

from rich.console import Console

from mlip_autopipec.common.errors import PhysicsViolationError
from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.explorers.md_engine import MDEngine
from mlip_autopipec.factories import create_generator, create_sampler
from mlip_autopipec.storage.database_manager import DatabaseManager

console = Console()

class PipelineOrchestrator:
    """
    Orchestrates the entire MLIP data generation pipeline from start to finish.
    """

    def __init__(self, config: FullConfig, db_path: Path):
        self.config = config
        self.db_path = db_path

    def run_pipeline(self):
        """
        Executes the full four-stage pipeline:
        1. Generation: Create initial seed structures.
        2. Exploration: Run MD simulations to generate trajectories.
        3. Sampling: Select a diverse subset of structures.
        4. Storage: Save the final dataset to a database.
        """
        # --- 1. Generation ---
        console.rule("[bold green]Stage 1: Generation[/bold green]")
        generator = create_generator(self.config)
        try:
            initial_structures = generator.generate()
            console.log(f"Generated {len(initial_structures)} initial structure(s).")
        except PhysicsViolationError as e:
            console.log(f"[bold red]Error during generation: {e}[/bold red]")
            return

        # --- 2. Exploration ---
        console.rule("[bold green]Stage 2: Exploration[/bold green]")
        explorer = MDEngine(self.config)
        # For now, we run exploration on the first generated structure.
        # A more advanced implementation would parallelize this.
        if not initial_structures:
            console.log("[yellow]No valid initial structures to explore. Aborting.[/yellow]")
            return

        trajectory = explorer.run(initial_structures[0])
        console.log(f"Exploration complete. Trajectory has {len(trajectory)} steps.")

        # --- 3. Sampling ---
        console.rule("[bold green]Stage 3: Sampling[/bold green]")
        sampler = create_sampler(self.config)
        sampled_structures = sampler.sample(trajectory)
        console.log(f"Sampled {len(sampled_structures)} structures from trajectory.")

        # --- 4. Storage ---
        console.rule("[bold green]Stage 4: Storage[/bold green]")
        with DatabaseManager(self.db_path) as db:
            db.write_structures(sampled_structures)
        console.log(f"Saved final dataset to [cyan]{self.db_path}[/cyan]")

        console.rule("[bold blue]Pipeline completed successfully![/bold blue]")
