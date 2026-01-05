"""Main pipeline orchestrator."""
from rich.console import Console

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.exploration.engine import MDExplorer
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.sampling.samplers import RandomSampler
from mlip_autopipec.storage.database import AseDBWrapper

console = Console()


class PipelineRunner:
    """Orchestrates the data generation pipeline."""

    def __init__(self, config: FullConfig):
        """
        Initialize the runner with the full configuration.

        Args:
            config: The FullConfig Pydantic model.
        """
        self.config = config

    def run(self) -> None:
        """Execute the full Generation -> Exploration -> Sampling -> Storage pipeline."""
        console.print("[bold cyan][1/4][/bold cyan] Starting structure generation...")
        generator = AlloyGenerator(self.config.system)
        structures = generator.generate()

        console.print("[bold cyan][2/4][/bold cyan] Executing exploration stage...")
        explorer = MDExplorer()
        structures = explorer.run_md(structures)

        console.print("[bold cyan][3/4][/bold cyan] Sampling structures...")
        sampler = RandomSampler(self.config.sampling)
        sampled_structures = sampler.sample(structures)

        console.print("[bold cyan][4/4][/bold cyan] Writing to database...")
        db_path = f"{self.config.project_name}.db"
        with AseDBWrapper(db_path) as db:
            db.write_structures(sampled_structures)

        console.print("[bold green]Pipeline complete.[/bold green]")
