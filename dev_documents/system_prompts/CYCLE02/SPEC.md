# CYCLE 02 SPEC.md: Advanced Features - Simulation, Sampling, and Storage

## 1. Summary

This document provides the detailed technical specification for the second and final development cycle of the MLIP-AutoPipe project. This cycle builds directly upon the foundational CLI, configuration, and generation services that were successfully established in Cycle 01. The primary, overarching goal of Cycle 02 is to implement the sophisticated scientific and data-handling capabilities that represent the core value proposition of the application. This will transform the tool from a simple structure generator into a fully functional, end-to-end pipeline capable of producing a curated, diverse, and analysis-ready database of atomic structures suitable for training high-quality Machine Learning Interatomic Potentials. The development in this cycle is strategically focused on three major new components: first, a powerful **Exploration service** that leverages Molecular Dynamics (MD) and Monte Carlo (MC) simulations to generate a vast ensemble of atomic configurations; second, an intelligent **Sampling service** that employs advanced algorithms to select the most scientifically valuable and diverse structures from the raw simulation trajectories; and third, a robust **Storage service** that saves the final, curated dataset into a structured and queryable ASE SQLite database.

Upon the successful completion of this cycle, the MLIP-AutoPipe will be feature-complete. It will deliver on its full promise, providing researchers with a powerful, automated tool to significantly accelerate the development of next-generation MLIPs. The implementation will involve tackling complex challenges, including the management of parallel processes for computationally intensive simulations, the efficient handling of large data volumes, and the integration of external scientific libraries for tasks like structural analysis and descriptor calculation. The architectural patterns of modularity and service-based design established in Cycle 01 will be critical to managing this complexity. Each new service will be a self-contained, testable unit, and the `WorkflowOrchestrator` will be extended to manage the full, four-stage data flow. This cycle represents the culmination of the project's design phase, turning the architectural blueprint into a tangible, powerful, and complete scientific software application.

## 2. System Architecture

The work in this cycle involves a significant expansion of the existing project structure. While the foundational architecture remains the same, we will be adding the new service modules for exploration, sampling, and storage. The existing `config.py` module will be modified to incorporate the configuration schemas for these new stages, and the `WorkflowOrchestrator` will be enhanced to manage the complete, four-stage pipeline. A new `domain/models.py` file will also be introduced to house the Pydantic models used for internal data transfer, ensuring type safety and clear contracts between the different service layers. The test suite will be correspondingly expanded with new files to ensure that each new service is thoroughly and independently verified.

**File Structure (Cycle 02 Additions/Modifications in bold):**

```
src/
└── mlip_autopipec/
    ├── __init__.py
    ├── cli.py
    ├── **config.py**              # **Modified** to add Exploration, Sampling, and Storage sections
    ├── **domain/**
    │   ├── **__init__.py**
    │   └── **models.py**          # New file for core Pydantic data models (e.g., DFTResult)
    ├── services/
    │   ├── __init__.py
    │   ├── generation.py
    │   ├── **exploration.py**     # New service for the MD/MC simulation engine
    │   ├── **sampling.py**        # New service for FPS and Random sampling logic
    │   └── **storage.py**         # New service for database interaction (AseDBWrapper)
    └── orchestrators/
        ├── __init__.py
        └── **workflow.py**        # **Modified** to orchestrate the full four-stage pipeline
tests/
└── unit/
    ├── __init__.py
    ├── test_cli.py
    ├── test_config.py
    ├── test_generation.py
    ├── **test_exploration.py**  # New test suite for the Exploration service
    ├── **test_sampling.py**     # New test suite for the Sampling service
    └── **test_storage.py**      # New test suite for the Storage service
```

This extended architecture fully realizes the service-oriented design envisioned in the project's planning phase. Each major scientific function of the application—Generation, Exploration, Sampling, and Storage—is now cleanly encapsulated in its own dedicated module within the `services/` directory. This high degree of modularity is crucial for long-term maintainability and extensibility. The `WorkflowOrchestrator` remains the high-level conductor, but its logic will be expanded to handle the more complex, multi-stage data flow, ensuring that the output of each service is correctly passed as the input to the next.

## 3. Design Architecture

The design for Cycle 02 remains firmly rooted in the schema-first approach established in Cycle 01. The `config.py` module will be significantly extended to define the schemas for all the new, advanced stages of the pipeline. The core of this cycle's design, however, lies in the detailed architecture of the complex new service modules, which will encapsulate the majority of the project's scientific and computational logic.

**`config.py` - Pydantic Schema Extensions Blueprint:**

The `FullConfig` Pydantic model will be extended to become the single, comprehensive configuration object for the entire application. This provides a unified, validated, and self-documenting interface for the user.

```python
# (Imports and GenerationConfig from Cycle 01 remain)

class ExplorationConfig(BaseModel):
    """Configuration for the thermodynamic exploration phase."""
    engine: str = Field("ASE_MD", description="The simulation engine to use.")
    num_steps: int = Field(1000, gt=0, description="The number of simulation steps to perform for each structure.")
    temperature_k: float = Field(300.0, gt=0, description="The target simulation temperature in Kelvin.")
    mlip_model: str = Field("MACE", description="The MLIP model to use for calculating forces and energy.")
    # ... other MD parameters like timestep, pressure, etc.

class SamplingConfig(BaseModel):
    """Configuration for the structure sampling phase."""
    method: str = Field("FPS", description="The sampling method to use ('Random' or 'FPS').")
    num_samples: int = Field(100, gt=0, description="The target number of structures to select for the final dataset.")
    # ... other sampling-specific parameters

class StorageConfig(BaseModel):
    """Configuration for the final data storage phase."""
    db_path: str = Field("output_dataset.db", description="The file path for the output ASE database.")

class FullConfig(BaseModel):
    """The root configuration model for the entire MLIP-AutoPipe pipeline."""
    generation: GenerationConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    storage: StorageConfig
```

**`services/exploration.py` - ExplorationEngine Detailed Design:**
This service is the most computationally intensive and complex part of the application.
-   **`ExplorationEngine` class:**
    -   `__init__(self, config: ExplorationConfig)`: The constructor will take the relevant, validated configuration object.
    -   `run(self, initial_structures_path: str) -> List[str]`: The main public method will accept the file path to the initial structures generated in Cycle 01. It will return a list of file paths to the generated trajectory files.
-   **Parallel Execution Management:** The `run` method will be built around Python's `concurrent.futures.ProcessPoolExecutor`. This will allow it to launch multiple, independent MD simulations in parallel, with each worker process handling one of the initial seed structures. This is essential for achieving high throughput on multi-core systems.
-   **"Late Binding" of MLIP Calculator:** A private helper method, `_initialize_calculator(model_name: str) -> ase.calculators.Calculator`, will be responsible for instantiating the actual MLIP calculator (e.g., loading a pre-trained MACE model). Crucially, this method will be called *inside* the worker function that is executed by each process, not in the main orchestrator process. This "late binding" is a critical design decision to avoid the often-fatal errors and significant performance overhead associated with attempting to pickle and transfer large PyTorch models between processes.
-   **Automatic Thermodynamic Ensemble Switching:** The worker function will contain logic to inspect the periodic boundary conditions (`pbc`) of the incoming `ase.Atoms` object. If all three dimensions are periodic (`[True, True, True]`), it will intelligently select an NPT (isothermal-isobaric) ensemble, which allows the simulation cell to relax. Otherwise (e.g., for a surface slab with `[True, True, False]`), it will use an NVT (canonical) ensemble. This prevents unphysical simulation artifacts, such as the collapse of vacuum layers in surface simulations.

**`services/sampling.py` - Sampler Detailed Design:**
-   **`BaseSampler` (ABC):** An abstract base class will define the sampling contract with a single abstract method: `sample(self, trajectory_paths: List[str]) -> List[ase.Atoms]`.
-   **`FPSSampler` class:**
    -   This class will implement the `sample` method. Its implementation will first read all the simulation trajectory files into a single, large list of `ase.Atoms` objects.
    -   It will then use a library like `dscribe` to compute SOAP (Smooth Overlap of Atomic Positions) descriptors for every single frame in the aggregated trajectory. The SOAP descriptor is a mathematical fingerprint that represents the local atomic environment, allowing for a quantitative comparison of structural similarity.
    -   Finally, it will implement the iterative Farthest Point Sampling algorithm. This algorithm will operate on the high-dimensional space of the SOAP descriptors to select a subset of `num_samples` structures that are maximally diverse, ensuring the final dataset is not redundant and covers the explored conformational space as efficiently as possible.

**`services/storage.py` - Storage Service Detailed Design:**
-   **`AseDBWrapper` class:** This class will provide a clean, high-level API for all database operations, abstracting away the low-level details of the `ase.db` library.
    -   `__init__(self, db_path: str)`: The constructor takes the path to the database file.
    -   `write_structures(self, structures: List[ase.Atoms])`: The main public method will iterate through the final list of curated structures. It will use a `with ase.db.connect(self.db_path) as db:` context manager to ensure the database connection is handled safely and efficiently. Inside a loop, it will call `db.write(atoms)` for each structure. It will also be responsible for attaching any relevant metadata, such as the calculated potential energy and atomic forces, to each database entry.

## 4. Implementation Approach

The implementation of Cycle 02 will be carefully staged to manage its inherent complexity. We will build and test each new service independently before integrating them into the main `WorkflowOrchestrator`.

1.  **Extend `config.py`:** The first step is to expand the existing configuration schema. The new Pydantic models—`ExplorationConfig`, `SamplingConfig`, and `StorageConfig`—will be added. The root `FullConfig` model will be updated to include these new models, creating the comprehensive schema for the entire end-to-end pipeline.
2.  **Implement `services/storage.py` (Build from the end):** We will begin by implementing the final stage of the pipeline.
    -   The `AseDBWrapper` class will be created.
    -   The `write_structures` method will be implemented. The core logic will be encapsulated within a `with ase.db.connect(self.db_path) as db:` context manager to ensure robust handling of the database connection. The implementation will loop through a list of `ase.Atoms` objects and call `db.write(atoms)` for each one.
3.  **Implement `services/exploration.py` (The core challenge):**
    -   The `ExplorationEngine` class will be created.
    -   The `run` method will be implemented, setting up the `ProcessPoolExecutor` to manage the parallel execution of the simulations.
    -   A dedicated, private worker function (`_run_single_md`) will be written. This function will be the target of the process pool and will contain the main MD simulation loop (e.g., using `ase.md.langevin.Langevin`). It will be responsible for initializing the calculator, running the dynamics, and saving the resulting trajectory to a unique file.
    -   The `_initialize_calculator` helper method will be implemented. For initial development and testing, this will be hardcoded to use a fast, simple potential like ASE's built-in EMT potential. This allows for rapid iteration before introducing the complexity of loading real MLIP models.
4.  **Implement `services/sampling.py` (The refinement step):**
    -   The `BaseSampler` abstract base class will be created to define the common interface.
    -   The `FPSSampler` will be implemented. This will be a multi-step process: first, reading all trajectory files; second, computing SOAP descriptors for all frames; and third, implementing the iterative FPS algorithm to select the most diverse structures.
5.  **Update `orchestrators/workflow.py` (The integration):**
    -   The `WorkflowOrchestrator`'s `run` method will be extended to execute the full, four-stage pipeline.
    -   It will now call each service in a strict sequence, managing the intermediate data (file paths) and passing the output of one stage as the input to the next:
        1.  `GenerationService` -> `initial_structures.xyz` (path)
        2.  `ExplorationEngine` -> `List[trajectory_paths]`
        3.  `SamplingService` -> `List[ase.Atoms]` (in-memory)
        4.  `AseDBWrapper` -> `output_dataset.db` (final output)
    -   The orchestrator will be responsible for creating and cleaning up temporary directories for the intermediate trajectory files.

## 5. Test Strategy

Testing in Cycle 02 is of paramount importance due to the introduction of complex, computationally intensive services and their intricate interactions. The strategy will rely heavily on the principle of isolation, using extensive mocking to create fast, deterministic, and reliable unit tests.

**Unit Testing Approach:**
-   **`test_storage.py`:** This suite will rigorously test the `AseDBWrapper`. To avoid filesystem dependencies, it will use an in-memory SQLite database by specifying the database path as `":memory:"`. The tests will call the `write_structures` method with a predefined list of `ase.Atoms` objects. After the write operation, it will use the ASE DB API directly to query the in-memory database and assert that the correct number of structures were written and that their stored properties (e.g., energy, atomic numbers, cell parameters) exactly match the input data.
-   **`test_exploration.py`:** Testing the `ExplorationEngine` presents the greatest challenge and will rely heavily on mocking. The external dependency—the MLIP calculator—will be completely mocked using `unittest.mock.patch`. The mock calculator will be configured with a simple, deterministic `get_forces` method that returns predefined force vectors. This allows us to test the engine's complex internal logic without the overhead and non-determinism of a real simulation. We will write specific tests to verify: (1) that the `ProcessPoolExecutor` is correctly instantiated and called, (2) that the automatic logic for selecting the NPT vs. NVT ensemble is correct based on the input structure's periodic boundary conditions, and (3) that the output trajectory file is created and contains the expected number of frames. We will not be testing the scientific correctness of the trajectory, only that the service correctly executes its orchestration workflow.
-   **`test_sampling.py`:** The `FPSSampler` will be tested in a deterministic manner. A small, pre-generated trajectory file containing a few simple, distinct structures will be used as a fixed input. We will pre-calculate the expected SOAP descriptors and the exact selection order that the FPS algorithm should produce for this input. The test will then run the sampler on this file and assert that it returns the exact list of `ase.Atoms` objects in the precise, expected order of diversity. This provides a rigorous validation of the sampling algorithm's implementation.

**Integration Testing Approach:**
The capstone of the Cycle 02 test strategy will be a full end-to-end integration test of the entire pipeline. This test will verify that all four services, now fully implemented, work together harmoniously.
1.  **Setup:** The test will create a temporary directory and a complete, valid `config.yaml` file. This configuration will be meticulously crafted for a fast test run: it will specify generating only one or two initial structures, running a very short MD simulation (e.g., 10-20 steps), and sampling a small number of final frames (e.g., 5). Critically, the `exploration` section of the config will specify the use of a very fast, non-MLIP calculator, such as ASE's built-in EMT potential, to ensure the entire test completes in seconds, not hours.
2.  **Execution:** The test will use the `typer.testing.CliRunner` to execute the main `run` command from the CLI, pointing it to the temporary configuration file. This single command will trigger the full, sequential execution of the pipeline: the `WorkflowOrchestrator` will call the real Generation, Exploration, Sampling, and Storage services.
3.  **Verification:** After the command completes successfully, the test will perform a comprehensive verification of the final output. It will first check for the existence of the `output_dataset.db` file, as specified in the configuration. It will then connect to this SQLite database file using `ase.db.connect`. A series of queries will be performed on the database to validate its contents. Detailed assertions will check: (1) that the number of rows in the database exactly matches the `num_samples` specified in the configuration, (2) that the stored structures have the correct chemical composition, and (3) that essential metadata like potential energy and atomic forces are present and have the correct format. This comprehensive test provides the highest level of confidence that the entire application dataflow is correct and that all components are properly integrated.
