"""The main orchestrator for the MLIP-AutoPipe data generation pipeline.

This module contains the `PipelineRunner`, the central class that manages the
entire workflow of the application. It is responsible for initializing all the
necessary components based on the user's configuration and executing each stage
of the pipeline in the correct sequence. This centralized orchestration simplifies
the application's logic and ensures a clear, predictable data flow.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

if TYPE_CHECKING:
    from mlip_autopipec.config.models import FullConfig


class PipelineRunner:
    """Orchestrates the entire data generation pipeline.

    This class acts as the main controller for the application. It takes the
    full, validated configuration and uses it to drive the four main stages:
    1. Generation: Create initial seed structures.
    2. Exploration: Run simulations to generate diverse configurations.
    3. Sampling: Select an information-rich subset of structures.
    4. Storage: Write the final dataset to a database.

    Attributes:
        config (FullConfig): The validated Pydantic model for the entire
            pipeline configuration.

    """

    def __init__(self, config: FullConfig) -> None:
        """Initialise the PipelineRunner with the full configuration.

        Args:
            config: The validated Pydantic model for the entire configuration.

        """
        self.config = config

    def run(self) -> None:
        """Execute the full Generation -> Exploration -> Sampling -> Storage pipeline."""
        logging.info("Starting MLIP-AutoPipe pipeline...")

        # 1. Generation Stage
        logging.info("[1/4] Starting structure generation...")
        generator = AlloyGenerator(config=self.config.system)
        initial_structures = generator.generate()
        logging.info(f"Generated {len(initial_structures)} initial structures.")

        # 2. Exploration Stage
        logging.info("[2/4] Executing exploration stage...")
        explorer = MDExplorer(config=self.config.exploration)
        explored_structures = explorer.run_md(structures=initial_structures)
        logging.info(f"Exploration complete. {len(explored_structures)} structures available.")

        # 3. Sampling Stage
        logging.info("[3/4] Sampling structures...")
        sampler = RandomSampler(config=self.config.sampling)
        sampled_structures = sampler.sample(structures=explored_structures)
        logging.info(f"Sampled {len(sampled_structures)} structures.")

        # 4. Storage Stage
        logging.info("[4/4] Writing to database...")
        db_path = f"{self.config.project_name}.db"
        with AseDBWrapper(db_path=db_path) as db:
            db.write_structures(structures=sampled_structures)
        logging.info(f"Successfully wrote structures to '{db_path}'.")

        logging.info("Pipeline complete.")
