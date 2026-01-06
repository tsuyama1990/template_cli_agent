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

logger = logging.getLogger(__name__)


class PipelineRunner:
    """The main orchestrator for the MLIP-AutoPipe pipeline."""

    def __init__(self, config: FullConfig):
        """
        Initialise the PipelineRunner.

        Args:
            config: The full, validated configuration object.
        """
        self.config = config

    def run(self) -> None:
        """Execute the full Generation -> Exploration -> Sampling -> Storage pipeline."""
        logger.info("Starting MLIP-AutoPipe pipeline.")

        # 1. Generation
        generator = AlloyGenerator(self.config.system)
        structures = generator.generate()

        # 2. Exploration
        explorer = MDExplorer(self.config.exploration)
        structures = explorer.explore(structures)

        # 3. Sampling
        sampler = RandomSampler(self.config.sampling)
        structures = sampler.sample(structures)

        # 4. Storage
        db_path = Path(f"{self.config.project_name}.db")
        with AseDBWrapper(db_path) as db:
            db.write_structures(structures)

        logger.info("MLIP-AutoPipe pipeline finished successfully.")
