import logging
from pathlib import Path

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.explorers import MDEngine
from mlip_autopipec.factories import create_generator, create_sampler
from mlip_autopipec.storage import DatabaseManager

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class PipelineOrchestrator:
    """
    Orchestrates the entire MLIP data generation pipeline.
    """

    def __init__(self, config: FullConfig) -> None:
        self.config = config
        self.output_dir = Path("output")
        self.output_dir.mkdir(exist_ok=True)

    def run(self) -> None:
        """Executes the four stages of the pipeline: Generate, Explore, Sample, Store."""
        logging.info("Starting MLIP-AutoPipe workflow...")

        # 1. Generation Stage
        logging.info("Stage 1: Generating initial structures...")
        generator = create_generator(self.config)
        initial_structures = generator.generate()
        logging.info(f"Generated {len(initial_structures)} initial structure(s).")

        # 2. Exploration Stage
        logging.info("Stage 2: Exploring structures via MD...")
        explorer = MDEngine(self.config.exploration)
        full_trajectory = []
        for i, structure in enumerate(initial_structures):
            logging.info(f"  - Exploring structure {i + 1}/{len(initial_structures)}")
            trajectory = explorer.explore(structure)
            full_trajectory.extend(trajectory)
        logging.info(f"Exploration resulted in a trajectory of {len(full_trajectory)} frames.")

        # 3. Sampling Stage
        logging.info("Stage 3: Sampling diverse structures from trajectory...")
        sampler = create_sampler(self.config)
        sampled_structures = sampler.sample(full_trajectory)
        logging.info(f"Sampled {len(sampled_structures)} structures.")

        # 4. Storage Stage
        logging.info("Stage 4: Storing final structures in the database...")
        db_path = self.output_dir / "final_structures.db"
        with DatabaseManager(db_path=str(db_path)) as db:
            db.write_structures(sampled_structures)
        logging.info(f"Successfully saved dataset to {db_path}")

        logging.info("MLIP-AutoPipe workflow finished successfully!")
