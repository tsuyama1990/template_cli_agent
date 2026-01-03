# Specification: CYCLE01 - Core Pipeline and CLI Foundation

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary and overriding objective of CYCLE01 is to construct the foundational architecture and implement a minimal, yet fully functional, end-to-end data generation pipeline. This initial cycle is arguably the most critical, as it will establish the core patterns, interfaces, and workflows upon which all future enhancements will be built. The principal deliverable of this cycle will be a robust command-line interface (CLI) tool. Using this CLI, a user will be able to define a material system via a configuration file and execute the complete four-stage workflow—Generation, Exploration, Sampling, and Storage—to produce a final, MLIP-ready dataset without any manual intervention during the run. This cycle is about building a solid, reliable, and extensible backbone for the entire project.

To manage complexity and ensure a high-quality initial delivery, the scope for CYCLE01 is deliberately and strategically focused. All functionality will be centered around a single, well-understood, and important class of materials: multi-component alloys. By concentrating on alloys, we can perfect the core data flow and orchestration logic of the pipeline without being distracted by the diverse and often complex physical constraints of other material types (like charge neutrality in ionic crystals). The exploration stage, the computational heart of the pipeline, will be implemented with a standard Molecular Dynamics (MD) simulation operating under the NVT (isothermal-isochoric, i.e., constant temperature and volume) ensemble. This provides a solid and widely understood baseline for exploring the potential energy surface. For the sampling stage, a straightforward random sampling approach will be implemented. While this method is less sophisticated than the advanced techniques planned for later cycles, it is sufficient to produce a valid and usable dataset, allowing us to validate the complete pipeline.

The key deliverables for this cycle form a complete, self-contained system. This includes the working CLI (`main.py`), the central `PipelineOrchestrator` which acts as the system's brain, a durable data persistence layer using a dedicated `AseDBWrapper`, and a highly structured, type-safe configuration system based on the powerful combination of Hydra and Pydantic. Alongside these are the concrete service implementations required for the alloy use case: the `AlloyGenerator` for creating initial structures, the `MDEngine` for running simulations, and the `RandomSampler` for selecting the final data points. Crucially, a comprehensive suite of unit and integration tests will be developed in parallel with the application code. This test-first mentality is non-negotiable and guarantees the correctness, reliability, and maintainability of this foundational system. By the successful completion of CYCLE01, a user will have in their hands a powerful command-line tool. They will be able to define an arbitrary alloy system (e.g., CuAu, NiCrAl) in a simple YAML configuration file, run a single command, and receive, as output, a professionally structured ASE database containing a valid training dataset. This represents the successful delivery of the project's core promise: end-to-end automation.

## 2. System Architecture

The architecture for CYCLE01 is designed to establish the fundamental, modular structure of the application, which is a blueprint for all future development. All work in this cycle will be contained within the `src/mlip_autopipec/` directory, creating a well-defined and installable Python package. The file structure is a direct reflection of the logical, layered architecture of the system, promoting a clear separation of concerns.

**File Structure for CYCLE01:**

The following ASCII tree provides a detailed blueprint of the files to be created or modified during this cycle. Files marked in **bold** represent the primary focus of implementation for CYCLE01. This structure is not merely a suggestion but a mandatory layout to ensure consistency and maintainability.

```
src/mlip_autopipec/
├── cli/
│   ├── __init__.py
│   └── **main.py**              # CLI entrypoint using Click. Handles user interaction.
├── core/
│   ├── __init__.py
│   ├── **factories.py**         # Contains Factory classes for component dispatch.
│   ├── **interfaces.py**        # Defines the Abstract Base Class contracts for all services.
│   └── **pipeline_orchestrator.py** # The central class that orchestrates the pipeline stages.
├── domain/
│   ├── __init__.py
│   └── **configuration.py**     # Contains all Pydantic models for Hydra configuration schema.
└── services/
    ├── __init__.py
    ├── generation/
    │   ├── __init__.py
    │   ├── **alloy.py**         # Concrete implementation for AlloyGenerator.
    │   └── **base.py**          # Defines the BaseGenerator ABC.
    ├── exploration/
    │   ├── __init__.py
    │   └── **md_engine.py**     # Concrete implementation for a basic MDEngine.
    ├── sampling/
    │   ├── __init__.py
    │   └── **random.py**        # Concrete implementation for RandomSampler.
    └── storage/
        ├── __init__.py
        └── **ase_db_wrapper.py** # Concrete implementation for the AseDBWrapper service.
└── utils/
    ├── __init__.py
    └── **physics.py**           # Utility functions for physical validation (e.g., overlap_check).
```

**Component Blueprint:**

This section provides an in-depth blueprint for the core components to be built in this cycle.

*   **`cli/main.py`:** This is the user's entry point to the application. It will be implemented using the `click` library to create a clean and well-documented command-line interface. A main command, `run-pipeline`, will be defined. This command will accept the path to the Hydra configuration directory and a set of overrides. Its role is strictly limited to parsing user inputs, instantiating the `PipelineOrchestrator`, triggering the workflow, and handling any high-level exceptions to provide user-friendly feedback. It will not contain any business logic.

*   **`core/interfaces.py`:** This file is the cornerstone of the modular architecture. It will define the abstract contracts (interfaces) that all service components must adhere to. Using Python's `abc` module, we will define:
    *   `IGenerator`: Specifies a `generate()` method that must return a list of ASE `Atoms` objects.
    *   `IExplorer`: Specifies an `explore()` method that takes a list of file paths to input structures and must return a list of file paths to the output trajectory files.
    *   `ISampler`: Specifies a `sample()` method that takes a list of trajectory file paths and must return a list of the final, sampled ASE `Atoms` objects.
    *   `IStorage`: Specifies a `store()` method that takes a list of final `Atoms` objects and writes them to the persistent database.
    This interface-driven design is what allows the orchestrator to be decoupled from the concrete implementations.

*   **`core/pipeline_orchestrator.py`:** The `PipelineOrchestrator` class is the heart of the application. It contains the high-level logic for the entire workflow. Its main `run()` method will execute the four stages in the correct sequence. It will manage the flow of data between stages, which, by design, will be file paths. This loose coupling via the file system is a deliberate choice to enhance robustness and allow for easy resumption. The orchestrator will implement the checkpointing logic by checking for the existence of key output files (e.g., `initial_structures.xyz`) before starting a stage. If the file exists, the stage is skipped. This class depends only on the interfaces from `core/interfaces.py`, not on any concrete service.

*   **`services/generation/alloy.py`:** The `AlloyGenerator` will be the first concrete implementation of the `IGenerator` interface. It will be responsible for creating random alloy structures. Its logic will take the `AlloyGeneratorConfig` Pydantic model as input and use libraries like ASE to construct a crystal lattice, create a supercell, and then randomly assign elemental identities to the atomic sites based on the requested composition. A crucial final step will be to pass every generated structure through the validation functions in `utils/physics.py` to ensure its physical plausibility.

*   **`services/exploration/md_engine.py`:** The `MDEngine` will implement the `IExplorer` interface. For CYCLE01, its functionality will be to perform basic NVT MD simulations. It will use a `ProcessPoolExecutor` to farm out the individual simulations to a pool of worker processes, achieving parallelism. The most critical technical detail is the "late binding" of the MLIP calculator. The `MDEngine` will not create the calculator; instead, the worker function that runs in the separate process will be responsible for instantiating it. This avoids the severe performance and reliability problems associated with pickling large PyTorch models. The engine will write the output trajectories to disk for the next stage.

*   **`services/sampling/random.py`:** The `RandomSampler` will be a simple and straightforward implementation of the `ISampler` interface. Its job is to read the trajectory files produced by the exploration stage and select a pre-defined number of frames (atomic configurations) at random from across all trajectories.

*   **`services/storage/ase_db_wrapper.py`:** The `AseDBWrapper` will implement the `IStorage` interface. It will provide a clean, high-level API for all database interactions, completely encapsulating the underlying `ase.db` connection handling. It will include methods to connect to the database (creating it if it doesn't exist) and to write ASE `Atoms` objects along with all their calculated properties (energy, forces) and relevant metadata.

## 3. Design Architecture

The design of the system is centered around the rigorous and proactive use of Pydantic-based schemas. This approach, often termed "schema-first design," ensures that the data structures and configuration, which are the lifeblood of the application, are robust, clear, self-documenting, and validated at the system's boundaries. This prevents a vast class of potential runtime errors and provides immediate, helpful feedback to the user on configuration mistakes.

**Pydantic Schema Design (`domain/configuration.py`):**

The entire configuration structure for the application will be meticulously defined by a set of nested Pydantic models. Hydra will be configured to use these models as the schema for its YAML files. This means that as soon as the configuration is loaded, it is parsed, validated, and converted into a fully type-hinted Python object, which is a significant improvement over accessing configuration via untyped dictionaries.

```python
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Literal

class AlloyGeneratorConfig(BaseModel):
    elements: List[str] = Field(..., min_length=1, description="List of element symbols.")
    composition: Dict[str, float] = Field(..., description="Dictionary mapping element to its fraction.")
    lattice_type: Literal["fcc", "bcc", "hcp"] = Field("fcc", description="The crystal lattice type.")
    supercell_size: int = Field(5, gt=0, description="Size of the supercell matrix (e.g., 5 means 5x5x5).")
    num_structures: int = Field(10, gt=0, description="Number of initial seed structures to generate.")

    @validator('composition')
    def composition_must_sum_to_one(cls, v):
        if not abs(sum(v.values()) - 1.0) < 1e-6:
            raise ValueError("Composition fractions must sum to 1.0")
        return v

class GeneratorConfig(BaseModel):
    name: Literal["alloy"] = Field("alloy", description="The name of the generator to use.")
    params: AlloyGeneratorConfig

class MDEngineConfig(BaseModel):
    calculator: str = Field("mace_mp_0", description="Name of the MLIP model calculator to use.")
    temperature_K: float = Field(300.0, gt=0, description="Simulation temperature in Kelvin.")
    steps: int = Field(1000, gt=0, description="Number of MD steps to run.")
    timestep_fs: float = Field(1.0, gt=0, description="Timestep for MD integration in femtoseconds.")

class ExplorerConfig(BaseModel):
    name: Literal["md"] = Field("md", description="The name of the explorer engine to use.")
    params: MDEngineConfig
    max_workers: int = Field(4, gt=0, description="Number of parallel processes for exploration.")

class RandomSamplerConfig(BaseModel):
    num_samples: int = Field(100, gt=0, description="Total number of frames to sample.")

class SamplerConfig(BaseModel):
    name: Literal["random"] = Field("random", description="The name of the sampler to use.")
    params: RandomSamplerConfig

class StorageConfig(BaseModel):
    db_path: str = Field("results.db", description="Path to the output ASE SQLite database.")

class FullConfig(BaseModel):
    generator: GeneratorConfig
    explorer: ExplorerConfig
    sampler: SamplerConfig
    storage: StorageConfig

    class Config:
        extra = "forbid" # Disallow extra fields in the config
```

**Key Invariants, Constraints, and Validation Rules:**
The Pydantic models are not just data containers; they are active validators that enforce the business rules of the application at the configuration stage.
*   **Composition Sum:** A custom validator on `AlloyGeneratorConfig.composition` will ensure that the fractional compositions provided by the user sum to 1.0 (within a small tolerance), preventing a common and logical error.
*   **Element Consistency:** Further validators will ensure that the `elements` list and the keys of the `composition` dictionary are consistent with each other.
*   **Positive Values:** The `Field(gt=0)` directive on numeric fields like `temperature_K` and `steps` ensures that these physical quantities are strictly positive, which is a logical necessity.
*   **Strictness:** All models will be configured with `extra="forbid"`. This is a powerful feature that prevents users from passing unknown or misspelled configuration parameters. If a user types `temperture_K` instead of `temperature_K`, Hydra/Pydantic will immediately raise an error, pointing out the mistake.

**Data Consumers and Producers:**
*   **Producers:** The primary producer of the configuration is the end-user, who authors the YAML files.
*   **Consumers:** The `PipelineOrchestrator` is the primary consumer. After Hydra loads and parses the YAML files, the orchestrator receives a single, fully-validated `FullConfig` object. It then delegates the nested configuration objects (e.g., `config.generator.params`) to the appropriate factories and service components. Each service, therefore, receives a configuration object whose validity has already been guaranteed.

**Versioning, Extensibility, and Backward-Compatibility:**
This schema-first design provides a clear path for future evolution. For instance, when we add an `IonicGenerator` in CYCLE02, we will create a new `IonicGeneratorConfig` Pydantic model. The `GeneratorConfig.params` field will then be updated to be a `Union[AlloyGeneratorConfig, IonicGeneratorConfig]`, and the `name` field will be used as the discriminator. This approach ensures that any new additions are strongly typed and that the system can clearly distinguish between different component configurations. It makes extending the system's capabilities a matter of adding new, validated schemas, which is a highly robust and maintainable approach that ensures backward compatibility.

## 4. Implementation Approach

The implementation strategy for CYCLE01 is a logical, bottom-up progression. This approach ensures that each layer of the application is built upon a solid, well-tested foundation. We will start with the most independent components and gradually assemble them into the final, integrated application.

1.  **Project Scaffolding and Dependencies:** The first step is to set up the `pyproject.toml` file. This will define the project metadata, establish the `src/mlip_autopipec` package structure, and declare the initial set of core dependencies: `click` for the CLI, `hydra-core` for configuration, `pydantic` for schema validation, `ase` and `numpy` for the scientific computations.
2.  **Schema-First: Implement Pydantic Models:** Before any logic is written, we will create all the configuration models in `domain/configuration.py` as detailed in the section above. This acts as a formal contract and a single source of truth for the shape of the configuration data throughout the entire application. Writing these first forces clarity on the required parameters for each component.
3.  **Core Utilities:** Next, we will implement the stateless utility functions in `utils/physics.py`. The first and most important of these is a robust `check_overlap` function, which takes an ASE `Atoms` object and a tolerance, and returns `True` if any two atoms are closer than that tolerance. This utility can be developed and tested in complete isolation.
4.  **Define the Contracts: Core Interfaces:** With the data schemas defined, the next step is to formalize the contracts for our services. We will define the `IGenerator`, `IExplorer`, `ISampler`, and `IStorage` abstract base classes in `core/interfaces.py`. This is a critical step in our dependency inversion strategy.
5.  **Build the Services (Bottom-up):** Now we implement the concrete service classes, starting with the one with the fewest dependencies.
    *   **`AseDBWrapper`:** Implement the storage service in `services/storage/ase_db_wrapper.py`. This component can be fully developed and tested with only `ase` as a dependency.
    *   **`AlloyGenerator`:** Implement the alloy generation service in `services/generation/alloy.py`. It will depend on the Pydantic config and the `check_overlap` utility.
    *   **`RandomSampler`:** Implement the simple random sampling service in `services/sampling/random.py`. This is another straightforward component that primarily interacts with the file system.
    *   **`MDEngine`:** This is the most complex service in this cycle. Its development in `services/exploration/md_engine.py` will focus on two main challenges: correctly setting up the ASE `Langevin` dynamics for an NVT simulation, and robustly managing the parallel execution using `ProcessPoolExecutor`, including the critical "late binding" of the MLIP calculator inside the worker function to ensure both performance and reliability.
6.  **Implement the Factories:** With the services now created, we will implement the `GeneratorFactory` (and others as needed) in `core/factories.py`. The factory's job is to map the configuration (e.g., `name: "alloy"`) to the correct concrete implementation (`AlloyGenerator`).
7.  **Assemble the Orchestrator:** Now that all the individual "gears" of the machine are built, we will assemble them in the `PipelineOrchestrator` in `core/pipeline_orchestrator.py`. This class will instantiate the necessary services via the factories and will contain the main `run()` method that calls the service methods in the correct sequence, implementing the checkpointing logic.
8.  **Create the User Entrypoint (CLI):** The final step is to implement the user-facing CLI in `cli/main.py`. This file will be very lean; it will primarily be responsible for initializing Hydra, which in turn loads the configuration and instantiates the `PipelineOrchestrator`, and then simply calls its `run()` method.

## 5. Test Strategy

Testing is not a separate phase but an integral part of the development process for CYCLE01. We will adopt a test-driven mindset, writing tests in parallel with the implementation to ensure correctness, prevent regressions, and provide a safety net for future refactoring. We will use `pytest` as the testing framework, `pytest-mock` for isolation, and `pytest-cov` for monitoring coverage.

**Unit Testing Approach (Min 300 words):**
The primary goal of unit testing in this cycle is to rigorously verify the correctness of each individual component in complete isolation from the rest of the system. This ensures that the building blocks of our application are solid. We will use `pytest-mock` extensively to replace any external dependencies (like other services, the file system, or even random number generation) with controlled fakes, known as "mocks."

*   **Pydantic Models (`domain/`):** Our Pydantic configuration models are a crucial part of our defence against user error, and their validation logic must be tested. We will write specific tests to verify that our custom validators work as expected. For instance, we will create a test that attempts to initialize `AlloyGeneratorConfig` with a `composition` dictionary whose values sum to 0.9 and assert that this raises a `pydantic.ValidationError` with a helpful message. Similarly, we will test boundary conditions, like providing a `temperature_K` of -100, and confirm that the validation fails as expected.

*   **Services (`services/`):** Each service will have its own dedicated test file to verify its logic.
    *   The `AseDBWrapper` will be tested against a temporary, in-memory SQLite database. Tests will confirm that it can correctly create a new database file, write a list of `Atoms` objects to it, and then read them back, asserting that the data read is identical to the data written.
    *   The `AlloyGenerator` will be tested to ensure its output is correct. We will provide a fixed configuration and assert that the generated `Atoms` objects have the correct total number of atoms and the correct elemental composition (e.g., exactly 50% Cu and 50% Au atoms). We will also mock the `check_overlap` utility to ensure it is being called for every generated structure.
    *   The `MDEngine` presents the most interesting unit testing challenge. We will mock the expensive MLIP calculator with a simple, fast dummy calculator (like ASE's built-in EMT potential). The focus of the test will not be on the correctness of the physics but on the orchestration of the simulation. We will mock `ProcessPoolExecutor` to verify that the engine attempts to start the correct number of worker processes. We will test the worker function itself, asserting that it correctly initializes the calculator and runs the dynamics object for the specified number of steps, and that it successfully writes an output trajectory file.
    *   The `RandomSampler` will be tested by providing it with a dummy, multi-frame trajectory file on a mock filesystem and asserting that it selects the correct number of samples and that the returned `Atoms` objects are valid.

*   **Orchestrator (`core/`):** The `PipelineOrchestrator` will be tested by mocking all the service interfaces (`IGenerator`, `IExplorer`, etc.). The tests will verify that the orchestrator calls the `generate`, `explore`, `sample`, and `store` methods in the correct sequence. We will also test its resume/checkpointing logic thoroughly by using a mock filesystem. For example, we will pre-create the `initial_structures.xyz` file and then run the orchestrator, asserting that the `generate` method of our mock generator is *not* called, while the `explore` method *is* called.

**Integration Testing Approach (Min 300 words):**
While unit tests are essential for verifying individual components, they cannot guarantee that those components will work together correctly. That is the role of integration tests. These tests are designed to run the entire pipeline, or significant parts of it, to ensure that the "seams" between the components are correct and that data flows through the system as expected. These tests are inherently slower and more comprehensive than unit tests.

*   **End-to-End CLI Test:** The cornerstone of our integration testing strategy will be a single, powerful test that invokes the application from the very top layer—the CLI—just as a user would. We will use the `click.testing.CliRunner` to programmatically call our `run-pipeline` command from within a `pytest` test function. This test will execute the entire four-stage pipeline on a small, well-defined, and computationally cheap test case.

*   **Test Environment:** To ensure the test is hermetic and does not interfere with the developer's system, it will be run in an isolated temporary directory, which is created and torn down automatically by `pytest`'s `tmp_path` fixture. Within this directory, the test will programmatically create a minimal set of Hydra configuration files needed for the test case.

*   **Performance Considerations:** A key requirement for an effective integration test is that it must be fast enough to run frequently, ideally as part of a CI/CD pipeline. Running a full MD simulation with a real MLIP would be far too slow. Therefore, the test configuration will be carefully crafted for speed. It will specify a very small system (e.g., a 2x2x2 supercell of a binary alloy), and it will instruct the `MDEngine` to use a very fast, non-ML calculator, such as ASE's built-in EMT (Effective Medium Theory) potential. The MD simulation itself will be configured to run for a minimal number of steps (e.g., 10-20), just enough to produce a multi-frame trajectory file.

*   **Comprehensive Assertions:** After the `CliRunner.invoke` call completes, the test will perform a comprehensive series of assertions on the state of the temporary directory. It will check:
    1.  That the expected intermediate files (`initial_structures.xyz`, `traj_0.xyz`, etc.) were created.
    2.  That the final ASE database (e.g., `results.db`) was created.
    3.  It will then use the `ase.db` library to connect to the output database and rigorously query its contents.
    4.  It will assert that the database contains the exact number of structures specified in the `num_samples` configuration parameter.
    5.  It will perform sanity checks on the data within the database, for example by fetching a row and asserting that it contains valid, non-zero energy and force data.

This single end-to-end test provides an exceptionally high degree of confidence that the entire system is functioning correctly, from the initial configuration parsing at the CLI layer all the way down to the final data persistence in the database.
