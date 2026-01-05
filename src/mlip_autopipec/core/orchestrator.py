"""The main orchestrator for the data generation pipeline."""
import logging
from pathlib import Path

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

# Configure basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class PipelineRunner:
    """Orchestrates the Generation -> Exploration -> Sampling -> Storage workflow."""

    def __init__(self, config: FullConfig):
        self.config = config
        self.db_path = Path(f"{config.project_name}.db")

    def run(self):
        """Executes the full data generation pipeline."""

        # 1. Generation Stage
        logging.info("Starting structure generation...")
        generator = AlloyGenerator(self.config.system)
        initial_structures = generator.generate()
        logging.info(f"Generated {len(initial_structures)} initial structures.")

        # 2. Exploration Stage
        logging.info("Starting exploration stage...")
        explorer = MDExplorer()
        explored_structures = explorer.run_md(initial_structures)
        logging.info(f"Exploration produced {len(explored_structures)} structures.")

        # 3. Sampling Stage
        logging.info("Starting sampling stage...")
        sampler = RandomSampler(self.config.sampling)
        sampled_structures = sampler.sample(explored_structures)
        logging.info(f"Sampled {len(sampled_structures)} structures.")

        # 4. Storage Stage
        logging.info(f"Writing final structures to database: {self.db_path}...")
        with AseDBWrapper(self.db_path) as db:
            db.write_structures(sampled_structures)
        logging.info("Pipeline complete.")
