"""Contains the main PipelineRunner class."""
from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper


class PipelineRunner:
    """The main orchestrator for the MLIP-AutoPipe pipeline."""

    def __init__(self, config: FullConfig):
        """Initialize the PipelineRunner."""
        self.config = config

    def run(self):
        """Run the full pipeline."""
        print("[1/4] Starting structure generation...")
        generator = AlloyGenerator(self.config.system)
        structures = generator.generate()

        print("[2/4] Executing exploration stage...")
        explorer = MDExplorer()
        structures = explorer.run_md(structures)

        print("[3/4] Sampling structures...")
        sampler = RandomSampler(self.config.sampling)
        structures = sampler.sample(structures)

        print("[4/4] Writing to database...")
        db_wrapper = AseDBWrapper(f"{self.config.project_name}.db")
        db_wrapper.write_structures(structures)

        print("Pipeline complete.")
