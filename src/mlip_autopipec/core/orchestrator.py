# -*- coding: utf-8 -*-
"""
Contains the main PipelineRunner class, which orchestrates the workflow.

This module is the heart of the application, connecting the individual
components (generation, exploration, sampling, storage) into a cohesive
pipeline.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

if TYPE_CHECKING:
    from mlip_autopipec.config.models import FullConfig

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class PipelineRunner:
    """
    The main orchestrator for the MLIP-AutoPipe workflow.

    This class is responsible for executing the entire data generation pipeline
    from start to finish. It initializes the necessary components based on the
    user's configuration and manages the flow of data between each stage.
    """

    def __init__(self, config: FullConfig) -> None:
        """
        Initialize the PipelineRunner.

        Args:
            config: The full, validated Pydantic configuration model for the run.
        """
        self.config = config

    def run(self) -> None:
        """
        Execute the full Generation -> Exploration -> Sampling -> Storage sequence.

        This method orchestrates the pipeline by calling each component in the
        correct order. It also includes error handling to gracefully manage
        unexpected issues during execution.

        Raises:
            Exception: Propagates any exception that occurs during the pipeline stages.
        """
        logging.info("Starting MLIP-AutoPipe pipeline...")

        try:
            # 1. Generation
            logging.info("[1/4] Starting structure generation...")
            generator = AlloyGenerator(config=self.config.system)
            structures = generator.generate()
            logging.info(f"Generated {len(structures)} initial structures.")

            # 2. Exploration (Placeholder)
            logging.info("[2/4] Executing exploration stage...")
            explorer = MDExplorer()
            structures = explorer.run_md(structures)
            logging.info("Exploration stage complete.")

            # 3. Sampling
            logging.info("[3/4] Sampling structures...")
            sampler = RandomSampler(config=self.config.sampling)
            sampled_structures = sampler.sample(structures)
            logging.info(f"Sampled {len(sampled_structures)} structures.")

            # 4. Storage
            logging.info("[4/4] Writing to database...")
            db_path = Path(f"{self.config.project_name}.db")
            with AseDBWrapper(db_path) as db:
                db.write_structures(sampled_structures)
            logging.info(f"Successfully wrote structures to {db_path}.")

            logging.info("Pipeline complete.")

        except Exception as e:
            logging.error(f"An unexpected error occurred during the pipeline execution: {e}")
            raise
