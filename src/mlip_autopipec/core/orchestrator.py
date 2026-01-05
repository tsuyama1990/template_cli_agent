"""The main orchestrator for the MLIP-AutoPipe pipeline."""

from pathlib import Path

from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import write_structures

console = Console()


class PipelineRunner:
    """Orchestrates the data generation pipeline."""

    def __init__(self, config: FullConfig) -> None:
        self.config = config

    def run(self):
        """Executes the full data generation pipeline."""
        console.print("[1/4] Starting structure generation...")
        generator = AlloyGenerator(self.config.system)
        structures = generator.generate()

        console.print("[2/4] Executing exploration stage...")
        explorer = MDExplorer()
        explored_structures = explorer.run_md(structures)

        console.print("[3/4] Sampling structures...")
        sampler = RandomSampler(self.config.sampling)
        sampled_structures = sampler.sample(explored_structures)

        console.print("[4/4] Writing to database...")
        db_path = Path(f"{self.config.project_name}.db")
        write_structures(db_path, sampled_structures)

        console.print("Pipeline complete.")
