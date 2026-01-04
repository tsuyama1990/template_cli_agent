import logging
from pathlib import Path

from ase.io import write

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.interfaces import IExplorer, ISampler, IStructureGenerator
from mlip_autopipec.storage.database_manager import DatabaseManager

logging.basicConfig(level=logging.INFO)


class PipelineOrchestrator:
    def __init__(
        self,
        config: FullConfig,
        generator: IStructureGenerator,
        explorer: IExplorer,
        sampler: ISampler,
        db_manager: DatabaseManager,
        output_dir: Path,
    ):
        self.config = config
        self.generator = generator
        self.explorer = explorer
        self.sampler = sampler
        self.db_manager = db_manager
        self.output_dir = output_dir
        self.output_dir.mkdir(exist_ok=True, parents=True)

    def run(self):
        """Executes the full data generation pipeline."""
        logging.info("Starting pipeline orchestration...")

        # 1. Generation
        logging.info("Stage 1: Generating initial structures...")
        initial_structures = self.generator.generate()
        initial_structures_path = self.output_dir / "initial_structures.xyz"
        write(initial_structures_path, initial_structures, format="extxyz")
        logging.info(f"Generated {len(initial_structures)} structures.")

        # 2. Exploration
        logging.info("Stage 2: Running MD simulations...")
        trajectory_frames = self.explorer.run(initial_structures)
        trajectory_path = self.output_dir / "trajectory.xyz"
        write(trajectory_path, trajectory_frames, format="extxyz")
        logging.info(f"MD simulation produced a trajectory of {len(trajectory_frames)} frames.")

        # 3. Sampling
        logging.info("Stage 3: Sampling structures from trajectory...")
        sampled_structures = self.sampler.sample(trajectory_frames)
        logging.info(f"Sampled {len(sampled_structures)} structures.")

        # 4. Storage
        logging.info("Stage 4: Writing sampled structures to database...")
        self.db_manager.write_atoms(sampled_structures)
        logging.info(f"Successfully wrote {len(sampled_structures)} structures to the database.")

        logging.info("Pipeline finished successfully.")
