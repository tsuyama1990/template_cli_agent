from mlip_autopipec.data.database import AseDBWrapper
from mlip_autopipec.data.models import FullConfig
from mlip_autopipec.modules.a_structure_generator import StructureGenerator
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine
from mlip_autopipec.modules.d_training_engine import TrainingEngine
from mlip_autopipec.orchestrator import Orchestrator


def initialize_and_run_workflow(config: FullConfig):
    """
    Initializes all components and runs the full MLIP-AutoPipe workflow.

    This function acts as the "composition root" of the application, wiring
    together all the necessary objects based on the provided configuration.

    Args:
        config: A fully validated FullConfig object.
    """
    # 1. Initialize the Database Wrapper
    # For now, we'll hardcode the database name. In a future cycle, this
    # could come from the config.
    db_wrapper = AseDBWrapper(db_path="mlip_autopipec.db")

    # 2. Initialize all the Engine Modules
    structure_generator = StructureGenerator(
        config=config.generation, db_wrapper=db_wrapper
    )
    labeling_engine = LabelingEngine(config=config.dft_compute)
    training_engine = TrainingEngine(config=config.mlip_training)

    # 3. Initialize the Orchestrator with all the engines
    orchestrator = Orchestrator(
        db_wrapper=db_wrapper,
        structure_generator=structure_generator,
        labeling_engine=labeling_engine,
        training_engine=training_engine,
    )

    # 4. Run the main workflow
    orchestrator.run_full_pipeline()
