# src/mlip_autopipec/pipeline/orchestrator.py
"""The main orchestrator for the MLIP-AutoPipe workflow."""

from rich.console import Console

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.explorers.md_engine import MDEngine
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.samplers.random import RandomSampler
from mlip_autopipec.storage.database_manager import DatabaseManager

console = Console()


class PipelineOrchestrator:
    """Coordinates the entire four-stage data generation pipeline."""

    def __init__(self, config: FullConfig) -> None:
        self.config = config
        # In a more advanced implementation, a factory would be used here
        self.generator = AlloyGenerator(self.config.system)
        self.explorer = MDEngine(self.config.exploration)
        self.sampler = RandomSampler(self.config.sampling)
        self.storage = DatabaseManager(self.config.db_path)

    def run(self) -> None:
        """Executes the full pipeline from generation to storage."""
        # Step 1: Generation
        console.print("Stage 1: Generating initial structures...")
        initial_structures = self.generator.generate()
        console.print(f"  > Generated {len(initial_structures)} seed structure(s).")

        # Step 2: Exploration
        console.print("Stage 2: Exploring structures via simulation...")
        trajectory = self.explorer.run(initial_structures)
        console.print(
            f"  > Simulation produced a trajectory of {len(trajectory)} structures.",
        )

        # Step 3: Sampling
        console.print("Stage 3: Sampling diverse structures...")
        sampled_structures = self.sampler.sample(trajectory)
        console.print(f"  > Sampled {len(sampled_structures)} final structures.")

        # Step 4: Storage
        console.print("Stage 4: Storing data in the database...")
        self.storage.write_structures(sampled_structures)
        console.print(f"  > Data saved to '{self.config.db_path}'.")
        console.print("\n[bold green]Pipeline finished successfully![/bold green]")
