# SPECIFICATION: Cycle 1 - Core CLI Pipeline

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to deliver a robust, fully functional command-line interface (CLI) that orchestrates the core data generation pipeline. This foundational version will serve as the Minimal Viable Product (MVP), establishing the architectural backbone of the entire application and delivering immediate value to users who are comfortable with command-line workflows. The core focus is on creating a complete, end-to-end process that is both scientifically sound and programmatically robust. This cycle will enable users to generate sophisticated training datasets for MLIPs through a completely automated, non-interactive workflow, laying the essential groundwork for the more advanced features, such as the web UI, planned in Cycle 2. The scope of this cycle is comprehensive, covering every stage from initial configuration to final storage. It begins with a user-defined YAML configuration file that specifies the chemical system and simulation parameters. The pipeline then proceeds through the sequential stages of initial structure **Generation**, where a diverse set of seed structures are created; MD-based **Exploration**, where these seeds are evolved under simulated thermodynamic conditions; simple **Sampling**, where a subset of structures is selected; and finally, persistent **Storage**, where the curated dataset is saved to a queryable database.

The key deliverable for Cycle 1 is a stable, reliable, and well-tested Python package that can be easily installed and executed from the command line. The software architecture will be a primary focus, establishing a clean, modular design with a strict separation of concerns. This layered architecture will distinguish between the high-level orchestration logic, the individual scientific services (generation, exploration, etc.), the underlying data models (defined with Pydantic for robustness), and the infrastructure components that interact with external systems like the filesystem and databases. This modularity is not merely an aesthetic choice; it is critical for the long-term health of the project, ensuring that the system is maintainable, extensible, and easy to debug. It allows for individual components to be tested in isolation and for new functionality to be added in the future without requiring a complete rewrite. For example, the foundation laid in this cycle will make it straightforward to add new structure generators or sampling algorithms in the future. Key components to be developed include a flexible factory for instantiating different structure generators, a parallelized MD exploration engine that can leverage multi-core CPUs, a baseline `RandomSampler`, and a dedicated wrapper for all interactions with the ASE database. Upon completion of this cycle, a computational scientist will be able to define a material system in a simple text file, execute a single command in their terminal, and receive a high-quality, diverse dataset, ready to be used for training a state-of-the-art Machine Learning Interatomic Potential.

## 2. System Architecture

The architecture for Cycle 1 is intentionally focused on the backend and the command-line execution flow, creating a solid foundation for future extensions. The user's primary interaction with the system is through the CLI, providing a YAML configuration file that dictates the behavior of the entire pipeline. The `PipelineOrchestrator` acts as the central nervous system, reading the configuration and directing the flow of data through the various service components. This design ensures a clear and linear data processing path, which is easy to understand, debug, and test. The architecture is explicitly designed to be stateless from the orchestrator's perspective; all state is managed through files on disk (intermediate structures, trajectories, and the final database), which makes the process robust against interruptions and allows for easy inspection of intermediate results. The use of factories to create components decouples the orchestration logic from the specific implementations of the services, allowing for a flexible system where new components can be added without modifying the core workflow. This adherence to SOLID principles is a cornerstone of the Cycle 1 architectural design.

**File Structure (Cycle 1 Focus):**

The following files represent the core deliverables to be created or significantly modified during this cycle. Bold entries indicate the files that will contain the bulk of the new code. The structure is designed to enforce the logical separation of concerns between presentation, orchestration, domain modeling, and service implementation.

```
src/mlip_autopipec/
├── __init__.py
├── **cli.py**                  # **Presentation Layer**: The main entry point for the CLI, built using the Click library. It handles argument parsing and initiates the pipeline.
│
├── core/                   # **Application Core / Orchestration Layer**: Manages the overall workflow and component lifecycle.
│   ├── __init__.py
│   ├── **pipeline_orchestrator.py** # Contains the central PipelineOrchestrator class that drives the A-Z process.
│   └── **factories.py**             # Implements the factory pattern to create service components (generators, samplers) based on the configuration.
│
├── domain/                   # **Domain Layer**: Contains the core data structures and business rules of the application.
│   ├── __init__.py
│   ├── **models.py**                # **Crucial File**: Defines all Pydantic models for configuration, ensuring type safety and validation.
│   └── **interfaces.py**            # Defines the abstract interfaces (ABCs) for services, establishing the contracts for modular components.
│
├── infrastructure/           # **Infrastructure Layer**: Provides wrappers for external systems and low-level operations.
│   ├── __init__.py
│   ├── **ase_db_wrapper.py**        # Encapsulates all interactions with the ASE SQLite database, providing a clean API.
│   └── **process_runner.py**        # A generic wrapper for executing external subprocesses, ensuring consistent error handling.
│
└── services/                 # **Service Layer**: Contains the concrete implementations of the scientific and business logic.
    ├── __init__.py
    ├── generation/              # Services related to initial structure generation.
    │   ├── __init__.py
    │   ├── **alloy_generator.py**   # Implementation for generating alloy structures.
    │   └── **ionic_generator.py**   # Implementation for generating charge-neutral ionic crystals.
    ├── exploration/             # Services related to MD-based exploration.
    │   ├── __init__.py
    │   └── **md_mc_explorer.py**    # The core MD engine, including parallel execution logic.
    └── sampling/                # Services related to dataset sampling.
        ├── __init__.py
        ├── **random_sampler.py**    # The baseline implementation of a random sampler.
```

**Component Blueprint:**

This blueprint details the specific responsibilities and implementation plan for each key file.

*   **`cli.py`**: This file will implement the user-facing command-line interface using the `click` library. It will define a main command, `mlip-autopipec run`, which will accept options like `--config` and `--output-db`. Its primary responsibilities are to validate the existence of the input config file, parse it using a helper function that leverages the Pydantic models, instantiate the main `PipelineOrchestrator`, and trigger the execution of the pipeline. It will also handle top-level exception catching and provide user-friendly error messages.

*   **`core/pipeline_orchestrator.py`**: This is the heart of the application's logic. The `PipelineOrchestrator` class will be initialized with a validated `FullConfig` object. Its main public method, `run()`, will orchestrate the entire workflow by calling the individual services in the correct sequence. It will manage the data flow by controlling the paths to intermediate files. For example, it will call the generator service, receive a list of `ase.Atoms` objects, save them to an `initial_structures.xyz` file, and then pass this file path to the explorer service. This explicit management of state on the filesystem makes the process transparent and robust.

*   **`services/generation/`**: This module will contain the logic for creating the initial seed structures. The `base_generator.py` file will define the `IStructureGenerator` abstract base class, which will enforce a `generate()` method signature for all implementations. The `alloy_generator.py` file will implement the logic for creating random alloy structures, handling various lattice types and compositions, and performing rigorous checks for atomic overlaps. Similarly, `ionic_generator.py` will implement a generator that ensures the created ionic compounds are charge-neutral.

*   **`services/exploration/md_mc_explorer.py`**: The `MDMCExplorer` will be the computational workhorse of the pipeline. Its core responsibility is to take a set of seed structures and run MD simulations on them. A key feature is its use of Python's `multiprocessing.ProcessPoolExecutor` to distribute the independent simulations across multiple CPU cores, dramatically speeding up the workflow. The implementation will carefully manage the "late binding" of the MLIP calculator, ensuring that the potentially large PyTorch-based model is loaded only within the worker processes to avoid serialization issues that can plague multiprocessing applications. This service will output trajectory files for each successful simulation.

*   **`services/sampling/random_sampler.py`**: For Cycle 1, a baseline `RandomSampler` will be implemented. This service will read the trajectory files produced by the explorer, iterate through all the frames, and randomly select a user-specified number of structures to be included in the final dataset. While simple, this provides the necessary functionality to complete the end-to-end pipeline and serves as a benchmark for the more advanced samplers to be introduced in Cycle 2.

## 3. Design Architecture

The design architecture for Cycle 1 is fundamentally centered around a **schema-first, Pydantic-driven** approach. This principle dictates that the data structures and configuration schemas are the most critical part of the design, and they must be defined first. By doing so, we establish a clear, validated, and self-documenting contract that the rest of the application can rely on. This approach significantly reduces the chances of runtime data-related errors and makes the entire system more robust and easier to reason about. All configuration files and data transfer objects within the application will be instances of these Pydantic models, guaranteeing that any data passed between components is valid and conforms to the expected schema. This is a modern software engineering practice that brings the benefits of static typing and validation into the often-dynamic world of Python applications.

*   **`domain/models.py`**: This file is the canonical source of truth for all data structures used in the application. It will contain a set of Pydantic models that meticulously define the entire configuration schema.
    *   **`BaseModelConfig`**: A base Pydantic model will be defined with `ConfigDict(extra="forbid")`. All other configuration models will inherit from this, which enforces a strict schema and causes the application to fail fast if the user provides any unknown fields in their YAML configuration file. This prevents common errors like typos in configuration keys.
    *   **`SystemConfig`**: This model will define the physical system to be simulated. It will have strongly typed fields: `elements: list[str]`, `composition: dict[str, float]`, `lattice_type: str`, and `num_initial_structures: int`. A custom validator will be implemented to ensure that the values in the `composition` dictionary sum to exactly 1.0, a critical physical invariant.
    *   **`ExplorationConfig`**: This model will specify all parameters for the MD simulation stage. Fields will include `model_path: str` (for the MLIP file), `temperature_k: float`, `pressure_bar: float | None` (an optional field for NPT simulations), `md_steps: int`, and `ensemble: Literal['nvt', 'npt']` to restrict the user to valid choices.
    *   **`SamplingConfig`**: For Cycle 1, this will be a simple model. It will include `method: Literal['random']` to enforce that only the random sampler can be used in this cycle, and `num_samples: int` to specify the size of the final dataset. A validator will ensure `num_samples` is a positive integer.
    *   **`FullConfig`**: This top-level model will compose the other models, creating a single, nested structure that represents a complete configuration for a pipeline run. It will have fields like `system: SystemConfig`, `exploration: ExplorationConfig`, and `sampling: SamplingConfig`.

*   **Data Flow and Consumers/Producers**:
    *   The **`cli.py`** module acts as the initial **producer** of the `FullConfig` object. It reads the raw YAML file from the user, and then uses Pydantic's parsing capabilities to deserialize it into a validated `FullConfig` instance.
    *   The **`PipelineOrchestrator`** is the primary **consumer** of the `FullConfig` object, using it to configure the entire workflow.
    *   The individual services are consumers of specific sub-sections of the configuration. The **Generator** services consume the `SystemConfig`, the **Explorer** service consumes the `ExplorationConfig`, and the **Sampler** service consumes the `SamplingConfig`. This clear flow of data ensures that each component only has access to the configuration it needs.
    *   The **`AseDBWrapper`** is the final consumer in the chain, taking the list of `ase.Atoms` objects produced by the sampler and persisting them to the database.

*   **Invariants and Validation**:
    *   Beyond the schema validation provided by Pydantic, the application will enforce other invariants at the service level. The most critical of these is the physical invariant enforced by the `BaseGenerator`: no two atoms in a generated structure can be closer than a specified minimum distance. The service will perform this check and raise a custom `PhysicsViolationError` if this invariant is broken, preventing unphysical structures from ever entering the pipeline. This combination of schema-level and service-level validation ensures a high degree of correctness and robustness for the entire system.

## 4. Implementation Approach

The implementation of Cycle 1 will follow a structured, bottom-up approach. This method ensures that the foundational components are built and tested before the higher-level orchestration logic is written. This strategy minimizes dependencies between developers working on different parts of the system and makes the development process more predictable and manageable.

1.  **Step 1: Define the Data Contracts (Pydantic Models)**: The very first task is to implement all the Pydantic models as specified in the Design Architecture section within the `domain/models.py` file. This includes creating the nested `FullConfig` model and adding all necessary field type hints and validators. This step is critical as it establishes the data contracts that the rest of the application will be built upon. A dedicated set of unit tests will be written for these models to ensure the validators work correctly (e.g., verifying that a composition that doesn't sum to 1.0 raises a `ValidationError`).

2.  **Step 2: Build the Infrastructure Layer**: With the data models defined, the next step is to build the low-level wrappers that interact with external systems.
    *   First, the `AseDBWrapper` will be implemented in `infrastructure/ase_db_wrapper.py`. This class will encapsulate all `ase.db` calls. It will have a clear, simple API: `connect(path)`, `write_atoms(list_of_atoms)`, and `disconnect()`. Unit tests for this component will use a temporary file to verify that data can be written and then read back perfectly.
    *   Next, a basic `ProcessRunner` will be implemented in `infrastructure/process_runner.py` to provide a standardized way of calling any potential external command-line tools, abstracting away the `subprocess` module.

3.  **Step 3: Implement Core Scientific Services (with Unit Tests)**: This is the most substantial part of the implementation. Each service will be developed in isolation and will have a comprehensive suite of unit tests.
    *   The `AlloyGenerator` and `IonicGenerator` will be implemented. Their unit tests will be crucial for verifying that the generated structures are always physically valid (e.g., no overlapping atoms, correct charge neutrality).
    *   The `RandomSampler` will be implemented. Its unit test will involve creating a dummy trajectory file and asserting that the sampler returns the correct number of `ase.Atoms` objects.
    *   The `MDMCExplorer` will be implemented. The focus here will be on the correct management of the `ProcessPoolExecutor`, including process startup, shutdown, and error handling. The unit tests will use a fast mock calculator to verify the parallel execution logic without the overhead of a real MLIP.

4.  **Step 4: Develop the Orchestration and CLI Layers**: Once all the individual services are built and tested, they will be assembled into a cohesive pipeline.
    *   The `PipelineOrchestrator` will be implemented in `core/pipeline_orchestrator.py`. This class will instantiate the required services (using the yet-to-be-built factory) and its `run()` method will call them in the correct sequence, managing the intermediate file I/O.
    *   The `factories.py` module will be implemented. It will contain functions that take the `FullConfig` object and return instances of the correct generator and sampler classes.
    *   Finally, the thin CLI layer will be implemented in `cli.py` using `click`. This will be the final piece that connects the user's command to the application's entry point.

5.  **Step 5: Write End-to-End Integration Tests**: After the entire pipeline is assembled, a final layer of integration tests will be written. These tests will use `click.testing.CliRunner` to invoke the CLI and run the entire pipeline with a minimal, fast-running configuration. The tests will assert that the final database is created and contains the expected data. This will provide the ultimate confidence that all components are correctly integrated and working together as designed.

## 5. Test Strategy

A comprehensive and multi-layered testing strategy is indispensable for ensuring the reliability, correctness, and scientific validity of the pipeline. The strategy for Cycle 1 is designed to build a foundation of trust in the software, combining rigorous unit tests for isolated components with holistic integration tests for the complete workflow.

**Unit Testing Approach (Min 300 words):**

The unit testing strategy is focused on verifying the correctness of each individual module or class in complete isolation from its dependencies. We will use the `pytest` framework for its powerful features and clear syntax, and `pytest-mock` for creating mock objects to replace external dependencies. This ensures that our unit tests are fast, deterministic, and precisely test the logic of the component in question.

*   **Pydantic Models**: Even our data models will have unit tests. We will write tests to specifically target the custom validators, ensuring they raise `pydantic.ValidationError` for invalid inputs. For example, we will test the `SystemConfig` model by providing a composition dictionary whose values do not sum to 1.0 and assert that the correct validation error is raised.

*   **Generators (`AlloyGenerator`, `IonicGenerator`)**: Testing the generators is critical for ensuring physical realism. The unit tests will have several key assertions. We will check that the returned `ase.Atoms` object has the correct number of atoms and the correct chemical composition as requested. Most importantly, we will have a utility function that checks for atomic overlap, and every test will assert that no two atoms are closer than a defined physical limit. We will use `pytest-mock` to mock Python's `random` module to ensure that the "random" structures generated during a test run are actually deterministic, making the tests reproducible. For the `IonicGenerator`, an additional assertion will confirm that the sum of oxidation states results in a neutral cell.

*   **Samplers (`RandomSampler`)**: Testing the `RandomSampler` is straightforward. The test setup will involve creating a dummy `ase.io.Trajectory` file containing a known number of atomic frames. The test will then instantiate the `RandomSampler` with a target sample count and call its `sample()` method. The primary assertion will be that the length of the returned list of `ase.Atoms` objects is exactly equal to the requested number of samples.

*   **Infrastructure (`AseDBWrapper`)**: The database wrapper will be tested against a temporary database file created in a temporary directory managed by `pytest`. The tests will cover the full CRUD (Create, Read, Update, Delete) lifecycle, although in our case, it's primarily Create and Read. A test will write a known list of `ase.Atoms` objects (with specific energies and forces) to the temporary database. It will then create a new instance of the wrapper, connect to the same database, read the structures back, and perform a detailed assertion to ensure that all the data (atomic numbers, positions, cell, energy, forces) is identical to the original data, guarding against data corruption or incorrect serialization.

**Integration Testing Approach (Min 300 words):**

While unit tests are essential for verifying individual components, they cannot guarantee that these components will work together correctly. The integration testing strategy is designed to fill this gap by validating the entire end-to-end workflow, from the command line to the final database output. These tests are inherently more complex and slower than unit tests, but they provide the highest level of confidence in the system's overall functionality.

*   **End-to-End CLI Test with `CliRunner`**: The centerpiece of our integration testing will be a test that simulates a real user running the application from the command line. We will use the `click.testing.CliRunner` utility, which allows us to invoke our CLI commands programmatically within a test. The test will be structured as follows:
    1.  **Setup**: The test will use a `pytest` fixture to create a temporary, isolated filesystem. It will then write a minimal but valid YAML configuration file into this directory. This configuration will be designed to run as quickly as possible: a very small system (e.g., 2-atom Si), a very short MD simulation (e.g., 10 steps), and the use of a fast classical potential instead of a slow MLIP.
    2.  **Execution**: The `CliRunner` will be used to invoke the `mlip-autopipec run` command, passing the path to the newly created configuration file. The test will capture the output and the exit code of the command.
    3.  **Assertions**: The test will perform a series of assertions to verify the successful execution of the pipeline. First, it will assert that the CLI command exited with a success code of 0 and that no exceptions were printed to the console. Next, it will inspect the temporary filesystem to verify that the expected intermediate files (e.g., `initial_structures.xyz`, `trajectory.traj`) and the final output database (`output.db`) were created. Finally, and most importantly, it will use our `AseDBWrapper` to connect to the output database, read the structures, and assert that the database contains the correct number of structures as specified in the configuration file's `num_samples` parameter. It will also perform a sanity check on one of the stored structures to ensure it contains valid atomic data.

This single, comprehensive test case validates the entire data flow through the orchestrated pipeline, confirming that the configuration is parsed correctly, the services are instantiated and called in the right order, and the final data is stored as expected. It is the ultimate guarantee that the core functionality of Cycle 1 is delivered.
