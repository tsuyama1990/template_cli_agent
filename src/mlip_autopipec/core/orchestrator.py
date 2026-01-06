"""The main orchestrator for the data generation pipeline."""
from __future__ import annotations

from typing import TYPE_CHECKING

from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

if TYPE_CHECKING:
    from mlip_autopipec.config.models import FullConfig


class PipelineRunner:
    """Orchestrates the entire data generation pipeline."""

    def __init__(self, config: FullConfig) -> None:
        """
        Initialise the PipelineRunner with the full configuration.

        Args:
            config: The validated Pydantic model for the entire configuration.
        """
        self.config = config

    def run(self) -> None:
        """
        Execute the full Generation -> Exploration -> Sampling -> Storage pipeline.
        """
        print("Starting MLIP-AutoPipe pipeline...")

        # 1. Generation Stage
        print("[1/4] Starting structure generation...")
        generator = AlloyGenerator(config=self.config.system)
        initial_structures = generator.generate()
        print(f"Generated {len(initial_structures)} initial structures.")

        # 2. Exploration Stage
        print("[2/4] Executing exploration stage...")
        explorer = MDExplorer(config=self.config.exploration)
        explored_structures = explorer.run_md(structures=initial_structures)
        print(f"Exploration complete. {len(explored_structures)} structures available.")

        # 3. Sampling Stage
        print("[3/4] Sampling structures...")
        sampler = RandomSampler(config=self.config.sampling)
        sampled_structures = sampler.sample(structures=explored_structures)
        print(f"Sampled {len(sampled_structures)} structures.")

        # 4. Storage Stage
        print("[4/4] Writing to database...")
        db_path = f"{self.config.project_name}.db"
        with AseDBWrapper(db_path=db_path) as db:
            db.write_structures(structures=sampled_structures)
        print(f"Successfully wrote structures to '{db_path}'.")

        print("Pipeline complete.")
