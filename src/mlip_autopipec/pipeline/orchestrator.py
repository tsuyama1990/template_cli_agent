# src/mlip_autopipec/pipeline/orchestrator.py
"""
The main orchestrator that runs the entire MLIP data generation pipeline.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

# For Cycle 1, we will use a simple, widely available calculator.
from ase.calculators.emt import EMT
from rich.console import Console

from mlip_autopipec import factories
from mlip_autopipec.storage.database_manager import DatabaseManager

if TYPE_CHECKING:
    from mlip_autopipec.common.pydantic_models import FullConfig


class PipelineOrchestrator:
    """
    Manages the sequential execution of the data generation workflow.
    """

    def __init__(self, config: FullConfig) -> None:
        """
        Initializes the orchestrator with a validated configuration.

        Args:
            config: The Pydantic model containing all pipeline settings.
        """
        self.config = config
        self.console = Console()

    def run(self) -> None:
        """
        Executes the full four-stage pipeline:
        1. Generation: Create initial seed structures.
        2. Exploration: Run MD simulations to generate trajectories.
        3. Sampling: Select a diverse subset of structures.
        4. Storage: Save the final dataset to a database.
        """
        self.console.print("[bold green]Starting MLIP-AutoPipe Pipeline...[/bold green]")

        # --- 1. Generation ---
        self.console.print("Step 1: [cyan]Generating initial structures...[/cyan]")
        generator = factories.create_generator(self.config)
        initial_structures = generator.generate()
        if not initial_structures:
            self.console.print(
                "[bold red]Error: No valid initial structures were generated.[/bold red]"
            )
            return
        self.console.print(f"Generated {len(initial_structures)} initial structure(s).")

        # --- 2. Exploration ---
        # For Cycle 1, we run exploration only on the first generated structure.
        self.console.print("Step 2: [cyan]Running MD exploration...[/cyan]")
        # In a real scenario, you'd inject a calculator like MACE.
        # EMT is a placeholder for testing and basic functionality.
        calculator = EMT()  # type: ignore[no-untyped-call]
        # The MDEngine is not yet implemented, so we will mock its output for now
        # This will be replaced in the next step.
        from mlip_autopipec.explorers.md_engine import MDEngine

        md_engine = MDEngine(config=self.config, calculator=calculator)
        trajectory = md_engine.explore(initial_structures[0])
        self.console.print(f"MD exploration produced a trajectory of {len(trajectory)} frames.")

        # --- 3. Sampling ---
        self.console.print("Step 3: [cyan]Sampling structures from trajectory...[/cyan]")
        sampler = factories.create_sampler(self.config)
        sampled_structures = sampler.sample(trajectory)
        self.console.print(f"Sampled {len(sampled_structures)} structures.")

        # --- 4. Storage ---
        self.console.print(
            f"Step 4: [cyan]Storing final dataset in '{self.config.db_path}'...[/cyan]"
        )
        db_manager = DatabaseManager()
        with db_manager.connect(self.config.db_path) as connection:
            db_manager.write_structures(connection, sampled_structures)

        self.console.print("[bold green]Pipeline completed successfully.[/bold green]")
