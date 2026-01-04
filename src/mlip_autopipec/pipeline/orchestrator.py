from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.storage.database_manager import DatabaseManager


class PipelineOrchestrator:
    """
    Manages the entire MLIP data generation pipeline.
    """

    def __init__(self, config: FullConfig) -> None:
        self.config = config
        self.generator = AlloyGenerator(config.system)
        self.db_manager = DatabaseManager("mlip_data.db")

    def run(self) -> None:
        """
        Executes the full pipeline.
        """
        # Step 1: Generation
        initial_structures = self.generator.generate(num_structures=10)

        # Step 2: Exploration (Placeholder)
        # In a real scenario, this would involve running MD simulations.
        # For now, we'll just use the initial structures.
        explored_structures = initial_structures

        # Step 3: Sampling (Placeholder)
        # Here, we would apply a sampling method like FPS.
        # For now, we'll take a subset of the explored structures.
        num_samples = self.config.sampling.number_of_samples
        sampled_structures = explored_structures[:num_samples]

        # Step 4: Storage
        self.db_manager.write_structures(sampled_structures)
