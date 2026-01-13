# Specification: MLIP-AutoPipe Cycle 1

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to deliver a robust, end-to-end, command-line-driven pipeline for the automated generation of training data for Machine Learning Interatomic Potentials (MLIPs). This foundational cycle focuses on implementing the core workflow, establishing a solid architectural base, and delivering immediate value to the user by providing a fully functional tool for creating high-quality datasets. The scope of this cycle encompasses the four key stages of the pipeline: initial structure **Generation**, basic thermodynamic **Exploration** using Molecular Dynamics (MD), simple random **Sampling**, and final **Storage** into a structured database.

The user will interact with the system exclusively through a Command-Line Interface (CLI), which will be built using the Typer or Click library. All aspects of the pipeline's behaviour will be controlled by a set of YAML configuration files, parsed and managed by the Hydra framework. This approach provides a powerful and flexible way for users to define the material system and the desired simulation parameters. The core philosophy is to build a reliable and deterministic system. Given the same configuration and random seed, the pipeline must produce the exact same output, ensuring reproducibility, which is a cornerstone of scientific research.

The architecture will be strictly modular, with a clear separation of concerns between the different components. A central `PipelineRunner` class will orchestrate the workflow, but the specific logic for each stage will be encapsulated within its own dedicated module. For instance, the `AlloyGenerator` will handle the creation of initial alloy structures, while the `MDEngine` will manage the execution of MD simulations. All data persistence will be centralised in a `WorkflowOrchestrator` service, which will be the sole component responsible for interacting with the file system and the final ASE database. This decoupling is crucial for creating a system that is easy to test, maintain, and extend in future cycles. By the end of this cycle, a user will be able to define a simple binary alloy system in a YAML file, run a single command, and receive a curated ASE database of atomic structures ready for MLIP training.

## 2. System Architecture

The architecture for Cycle 1 is focused on building the foundational command-line pipeline. The file structure is designed to be modular and scalable, laying the groundwork for the more advanced features that will be added in Cycle 2. All development will take place within the `src/mlip_autopipec/` directory.

The file structure below details the key files to be created or modified in this cycle. Files marked in **bold** are the primary targets for creation and implementation during Cycle 1.

```text
src/mlip_autopipec/
├── __init__.py
├── **cli.py**                  # Main CLI entry point using Click/Typer
│
├── core/
│   ├── __init__.py
│   ├── **models.py**           # Core Pydantic data models (e.g., DFTResult, TrainingConfig)
│   ├── **pipeline.py**         # The main PipelineRunner and WorkflowOrchestrator classes
│   └── **interfaces.py**       # Abstract base classes for engines/generators
│
├── components/
│   ├── __init__.py
│   ├── generators/
│   │   ├── __init__.py
│   │   ├── **base.py**         # BaseGenerator abstract class
│   │   └── **alloy.py**        # Concrete implementation for alloys
│   ├── exploration/
│   │   ├── __init__.py
│   │   └── **md_engine.py**    # Core MD simulation logic
│   ├── sampling/
│   │   ├── __init__.py
│   │   └── **random.py**       # Random sampling implementation
│   └── storage/
│       ├── __init__.py
│       └── **database.py**     # ASE DB wrapper and data persistence logic
│
└── utils/
    ├── __init__.py
    └── **config_loader.py**    # Hydra configuration management utilities
```

**Architectural Blueprint:**

1.  **`cli.py`**: This file will serve as the main entry point for the application. It will use Typer to define a simple CLI command, for example, `mlip-autopipec run --config-path <path>`. This script will be responsible for initialising the Hydra configuration, setting up logging, instantiating the `PipelineRunner`, and executing the main pipeline workflow.

2.  **Pydantic Models (`core/models.py`)**: This is the cornerstone of our schema-first design. We will define a series of nested Pydantic models that precisely mirror the structure of our YAML configuration files. This provides a single source of truth for all settings and ensures that any user-provided configuration is rigorously validated before the pipeline starts. For example, we will have a `SystemConfig` model to define the chemical elements and composition, and an `ExplorationConfig` model to define the MD simulation parameters like temperature and pressure. We will also define data-carrying models like `DFTResult` to standardize how simulation results are passed between components.

3.  **Pipeline Orchestration (`core/pipeline.py`)**: This file will contain two key classes: `PipelineRunner` and `WorkflowOrchestrator`.
    *   `PipelineRunner`: This class orchestrates the high-level workflow. It will be initialised with the main configuration object. It will have methods like `run()`, which will sequentially call the generation, exploration, sampling, and storage stages. It acts as the "conductor" of the pipeline.
    *   `WorkflowOrchestrator`: This class will handle all I/O operations. It will be responsible for creating directories, writing intermediate files (like the initial `.xyz` structures), and saving the final results to the ASE database. By centralizing I/O here, we keep our computational components (the "engines") pure and testable.

4.  **Interfaces (`core/interfaces.py`)**: To ensure a modular and extensible design, we will define abstract base classes (ABCs) for our key components. For example, we will define a `BaseStructureGenerator` interface with an abstract `generate()` method. This ensures that any new generator we add in the future will conform to the same contract, making it instantly compatible with the `PipelineRunner`.

5.  **Component Implementations (`components/`)**:
    *   `alloy.py`: This will be our first concrete implementation of a structure generator. The `AlloyGenerator` class will inherit from `BaseStructureGenerator` and implement the logic for creating random alloy structures, including performing physical validity checks.
    *   `md_engine.py`: This will contain the `MDEngine` class. It will be responsible for taking a set of input structures, setting up an ASE `Atoms` object with the correct calculator (e.g., MACE), and running an MD simulation using ASE's dynamics modules. It will be designed to run multiple simulations in parallel using Python's `multiprocessing` library.
    *   `random.py`: This module will implement a simple `RandomSampler` class, which will read a trajectory file and select a specified number of frames at random.
    *   `database.py`: This module will contain the `DatabaseManager` class, providing a clean wrapper around the ASE database functionality. It will have methods like `connect()` and `write_structures()`.

This architecture ensures a clean separation of concerns, making the system robust, testable, and ready for future expansion.

## 3. Design Architecture

The design of MLIP-AutoPipe Cycle 1 is firmly rooted in a schema-first philosophy, with Pydantic serving as the bedrock for data validation, configuration management, and defining the data contracts between different modules. This approach ensures that the system is robust, self-documenting, and less prone to runtime errors.

**Pydantic-Based Schema Design:**

1.  **Configuration Schema (`core/models.py`)**: The entire configuration hierarchy will be defined as a set of nested Pydantic `BaseModel` classes. This is the single source of truth for the system's settings.
    *   `FullConfig`: The top-level model that aggregates all other configuration models.
    *   `SystemConfig`: Defines the physical system. Its fields will include `elements: list[str]`, `composition: dict[str, float]`, and `crystal_structure: str`. It will have validators to ensure, for example, that the elements are valid chemical symbols and that the composition fractions sum to 1.0.
    *   `GenerationConfig`: Defines parameters for the initial structure generation. It will include fields like `num_initial_structures: int` with validation `gt=0` (greater than zero).
    *   `ExplorationConfig`: Contains all MD simulation parameters, such as `temperature: float`, `pressure: float`, `ensemble: Literal["NVT", "NPT"]`, and `mlip_model_path: FilePath`. The use of `Literal` and `FilePath` from Pydantic provides powerful, automatic validation.
    *   `SamplingConfig`: Defines the sampling strategy and its parameters. In Cycle 1, this will be simple, e.g., `method: Literal["random"]` and `num_samples: int`.

2.  **Data Transfer Objects (DTOs)**: We will also use Pydantic models as DTOs to pass structured data between components.
    *   `DFTResult`: A key model to represent the output of a single-point calculation. It will contain fields for `energy: float`, `forces: npt.NDArray[np.float64]`, and `stress: npt.NDArray[np.float64]`. We will implement custom JSON encoders/decoders to handle the serialization of NumPy arrays, ensuring that this structured data can be easily passed between processes or potentially stored.

**Key Invariants and Constraints:**
*   **Immutability**: Configuration objects, once loaded and validated, should be treated as immutable throughout the pipeline's execution to prevent unexpected side effects.
*   **Physical Validity**: The `AlloyGenerator` must enforce the constraint that no two atoms in a generated structure can be closer than a specified threshold (e.g., 0.7 times the sum of their covalent radii). This is a critical invariant for ensuring the physical realism of the initial structures.
*   **Data Integrity**: The `WorkflowOrchestrator` is the sole gatekeeper for data persistence. No other component is allowed to directly write to the final database. This ensures that all stored data has passed through the full pipeline and is consistent.

**Consumers and Producers:**
*   The **User** is the primary producer of the initial configuration (as YAML files).
*   The **`cli.py`** module consumes the user's configuration, validates it using the Pydantic models, and produces a `FullConfig` object.
*   The **`PipelineRunner`** consumes the `FullConfig` object.
*   The **Generators** produce a list of initial `ase.Atoms` objects.
*   The **`MDEngine`** consumes the `Atoms` objects and produces trajectory files containing many more `Atoms` objects.
*   The **Sampler** consumes the trajectory files and produces a smaller, selected list of `Atoms` objects.
*   The **`WorkflowOrchestrator`** consumes this final list and produces the populated ASE database.

**Versioning and Extensibility:**
The use of abstract base classes in `core/interfaces.py` is our primary strategy for ensuring extensibility. For example, in Cycle 2, we can create a new `FPS_Sampler` class that inherits from a `BaseSampler` interface. As long as it implements the required `sample()` method, the `PipelineRunner` can use it without any modification, thanks to polymorphism. Configuration versioning will be managed implicitly by the structure of the Pydantic models. Any backward-incompatible changes to the configuration schema would require a major version bump of the software itself.

## 4. Implementation Approach

The implementation of Cycle 1 will proceed in a logical, step-by-step manner, starting with the foundational components and progressively building up the full pipeline. The development will be guided by our schema-first and test-driven principles.

**Step 1: Project Scaffolding and Core Schema (`core/models.py`)**
The first task is to set up the directory structure as outlined in the System Architecture section. We will create the empty Python files and `__init__.py` files to define the package structure. Immediately after, we will implement the full Pydantic configuration schema in `core/models.py`. This is the most critical first step, as it defines the "language" of the application. We will create all the `...Config` models, complete with type annotations and validators. This allows us to define the expected inputs and outputs for every component before writing their implementation.

**Step 2: CLI and Configuration Loading (`cli.py`, `utils/config_loader.py`)**
With the schema defined, we will build the user's entry point. In `cli.py`, we will use Typer to create a simple `run` command. This command will take the path to the configuration directory as an argument. We will then implement the configuration loading logic in `utils/config_loader.py`, which will use Hydra to parse the YAML files and instantiate our `FullConfig` Pydantic model. At this stage, we can already test that the CLI can correctly load and validate configuration files, providing clear error messages for invalid inputs.

**Step 3: Database and Orchestrator (`storage/database.py`, `core/pipeline.py`)**
Next, we will focus on the end of the pipeline: storage. We will implement the `DatabaseManager` in `storage/database.py`, which will be a thin wrapper around the `ase.db` functionality. Then, we will create the `WorkflowOrchestrator` in `core/pipeline.py`. This class will be responsible for all file and database I/O. We will write its methods for creating directories and for connecting to and writing data into the database. By implementing this persistence layer early, we can use it in the testing of subsequent components.

**Step 4: Structure Generation (`components/generators/`)**
Now we will implement the first active component of the pipeline. We will define the `BaseStructureGenerator` interface in `core/interfaces.py`. Then, we will create the `AlloyGenerator` in `components/generators/alloy.py`. This class will contain the logic to generate random alloy structures based on the `SystemConfig`. The implementation will heavily rely on the ASE library for creating lattices and manipulating `Atoms` objects. A key part of this step is to implement the physical validation checks, especially the `overlap_check`.

**Step 5: MD Exploration (`components/exploration/md_engine.py`)**
This is the most computationally intensive part of the pipeline. We will implement the `MDEngine`. Its main responsibility is to take a structure, attach an MLIP calculator (like MACE, loaded from a file specified in the config), and run an MD simulation using one of ASE's dynamics algorithms (e.g., `Langevin`). We will implement this to be parallelizable from the start, using Python's `multiprocessing.Pool` to run simulations for different initial structures concurrently. Each process will write its output trajectory to a separate file to avoid race conditions.

**Step 6: Sampling and Tying it all together (`components/sampling/random.py`, `core/pipeline.py`)**
Finally, we will implement the `RandomSampler` in `components/sampling/random.py`. This will be a straightforward module that reads the trajectory files produced by the MD engine and randomly selects a configured number of frames. With all the individual components built, the last step is to implement the main `run()` method in the `PipelineRunner` class in `core/pipeline.py`. This method will instantiate and call each component in the correct order: Generator -> MDEngine -> Sampler -> WorkflowOrchestrator. This final step connects all the pieces into a functioning end-to-end pipeline.

## 5. Test Strategy

The test strategy for Cycle 1 is focused on ensuring the correctness and reliability of the core pipeline and its individual components. We will use `pytest` as our testing framework and aim for high test coverage for all new code. Our strategy is divided into two main pillars: unit testing and integration testing.

**Unit Testing Approach (Min 300 words):**
The goal of unit testing is to verify each piece of the application in isolation. This allows us to pinpoint failures accurately and makes the tests fast to run. We will make extensive use of mocking to isolate the unit under test from its dependencies.

1.  **Pydantic Models (`core/models.py`):** We will write a dedicated test suite for our configuration models. These tests will not check business logic but will verify the validation rules. For example, we will test that `SystemConfig` raises a `pydantic.ValidationError` if the composition percentages do not sum to 1.0, or if an invalid element symbol is provided. This ensures the robustness of our configuration system.

2.  **Generators (`components/generators/alloy.py`):** The `AlloyGenerator` will be tested by calling its `generate()` method and inspecting the returned `ase.Atoms` objects. We will assert that the number of structures is correct, that each structure has the right number and type of atoms, and that the cell dimensions are as expected. A critical test will be to verify the physical validation: we will programmatically check that no two atoms in the output structures are closer than the allowed threshold. We will mock any direct calls to `ase.build` to make the tests independent of the underlying library's implementation details.

3.  **Orchestrator (`core/pipeline.py`):** The `WorkflowOrchestrator` will be tested using a temporary, in-memory SQLite database. We will write tests that call its methods to write a list of mock `ase.Atoms` objects and then immediately read them back, asserting that the data (positions, energy, forces) is preserved perfectly. This ensures our database logic is correct without touching the actual file system.

4.  **MD Engine (`components/exploration/md_engine.py`):** Testing the MD engine requires careful mocking. We will not run actual MD simulations. Instead, we will mock the ASE calculator and the dynamics objects. The tests will focus on verifying the *setup* of the simulation. For example, we will test that the `MDEngine` correctly creates an ASE `Langevin` dynamics object with the temperature specified in the configuration. We will provide a mock `Atoms` object and assert that the correct calculator is attached to it.

**Integration Testing Approach (Min 300 words):**
Integration tests are designed to verify that the different, individually tested components work together as a cohesive whole. For Cycle 1, we will have a small number of high-value integration tests that cover the entire CLI workflow.

1.  **End-to-End CLI Test:** The primary integration test will simulate a user running the pipeline from the command line. We will use the `click.testing.CliRunner` (or a similar tool for Typer) to invoke the `cli.py` script from within a pytest test.
    *   **Setup:** The test will begin by creating a temporary directory using `pytest`'s `tmp_path` fixture. Inside this directory, it will create a minimal set of valid YAML configuration files for a very small system (e.g., a 4-atom Si-Ge alloy). It will also create a dummy MLIP model file for the calculator to load.
    *   **Execution:** The test will then call `runner.invoke()` to run the `mlip-autopipec run ...` command, pointing it to the temporary configuration directory. We will run with a very short MD simulation (e.g., 5 steps) to ensure the test completes quickly.
    *   **Verification:** After the command finishes, the test will inspect the contents of the temporary directory. It will assert that:
        a. An `initial_structures.xyz` file was created and is not empty.
        b. A trajectory file (e.g., `md_run_0.traj`) was created.
        c. A final ASE database file (e.g., `final_structures.db`) was created.
    *   **Database Inspection:** The test will then connect to the generated `.db` file and query its contents. It will assert that the database contains the correct number of structures, as specified in the sampling configuration. It will also retrieve one of the structures and verify that its metadata (e.g., energy, forces) is present.

This end-to-end test provides a high degree of confidence that the entire pipeline is wired together correctly, from configuration parsing to final database storage. It validates the data flow between all the major components of the system.
