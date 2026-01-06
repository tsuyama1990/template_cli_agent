"""
The central coordinating class for the MLIP-AutoPipe workflow.
"""

import logging
from pathlib import Path
from tempfile import TemporaryDirectory

from mlip_autopipec.common.exceptions import (
    ExplorationError,
    GenerationError,
    SamplingError,
    StorageError,
)
from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.explorers.md_engine import MDEngine
from mlip_autopipec.interfaces import IGenerator, ISampler
from mlip_autopipec.storage.database_manager import DatabaseManager

# Configure logger for this module
logger = logging.getLogger(__name__)

def _raise_generation_error() -> None:
    msg = "Generation stage produced no initial structures."
    raise GenerationError(msg)

def _raise_exploration_error() -> None:
    msg = "Exploration stage produced no trajectory frames."
    raise ExplorationError(msg)

class PipelineOrchestrator:
    """
    Manages the sequential execution of the four-stage pipeline.
    """

    def __init__(
        self,
        config: FullConfig,
        generator: IGenerator,
        explorer: MDEngine,
        sampler: ISampler,
        storage: DatabaseManager,
    ) -> None:
        self.config = config
        self.generator = generator
        self.explorer = explorer
        self.sampler = sampler
        self.storage = storage
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    def run_pipeline(self) -> None:
        """
        Executes the full four-stage pipeline with detailed error handling.
        """
        logger.info("🚀 Starting MLIP-AutoPipe pipeline...")

        # Stage 1: Generation
        try:
            logger.info("--- Stage 1: Generation ---")
            initial_structures = self.generator.generate()
            if not initial_structures:
                _raise_generation_error()
            logger.info(f"Generated {len(initial_structures)} initial structure(s).")
        except Exception as e:
            logger.error(f"Error in Generation stage: {e}", exc_info=True)
            msg = f"Generation failed: {e}"
            raise GenerationError(msg) from e

        with TemporaryDirectory() as temp_dir:
            trajectory_path = Path(temp_dir) / "trajectory.traj"

            # Stage 2: Exploration
            try:
                logger.info("--- Stage 2: Exploration ---")
                full_trajectory = self.explorer.explore(initial_structures, trajectory_path)
                if not full_trajectory:
                    _raise_exploration_error()
                logger.info(f"Exploration produced a trajectory with {len(full_trajectory)} frames.")
            except Exception as e:
                logger.error(f"Error in Exploration stage: {e}", exc_info=True)
                msg = f"Exploration failed: {e}"
                raise ExplorationError(msg) from e

            # Stage 3: Sampling
            try:
                logger.info("--- Stage 3: Sampling ---")
                sampled_structures = self.sampler.sample(full_trajectory)
                if not sampled_structures:
                    logger.warning("Sampling stage selected no structures. This may be expected.")
                else:
                    logger.info(f"Sampled {len(sampled_structures)} final structures.")
            except Exception as e:
                logger.error(f"Error in Sampling stage: {e}", exc_info=True)
                msg = f"Sampling failed: {e}"
                raise SamplingError(msg) from e

            # Stage 4: Storage
            try:
                logger.info("--- Stage 4: Storage ---")
                self.storage.write_structures(sampled_structures)
            except Exception as e:
                logger.error(f"Error in Storage stage: {e}", exc_info=True)
                msg = f"Storage failed: {e}"
                raise StorageError(msg) from e

        logger.info("✅ Pipeline completed successfully.")
