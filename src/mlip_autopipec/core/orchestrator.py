from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

if TYPE_CHECKING:
    from mlip_autopipec.config.models import FullConfig


class PipelineRunner:
    """Orchestrates the MLIP-AutoPipe workflow."""

    def __init__(self, config: FullConfig):
        self.config = config

    def run(self) -> None:
        """Executes the full data generation pipeline."""
        # 1. Generation
        generator = AlloyGenerator(self.config.system)
        structures = generator.generate()

        # 2. Exploration
        explorer = MDExplorer()
        structures = explorer.run_md(structures)

        # 3. Sampling
        sampler = RandomSampler(self.config.sampling)
        sampled_structures = sampler.sample(structures)

        # 4. Storage
        db_path = Path(f"{self.config.project_name}.db")
        storage = AseDBWrapper(db_path)
        storage.write_structures(sampled_structures)

        print(f"Pipeline complete. Database saved to {db_path}")
