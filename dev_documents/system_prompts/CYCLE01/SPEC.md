# Specification: CYCLE01 - Core Pipeline and Foundational Components

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to establish a solid, functional foundation for the entire application. This initial phase is strategically focused on building the essential "plumbing" and architectural skeleton upon which the more complex scientific capabilities will be built in the subsequent cycle. Rather than implementing the computationally intensive algorithms at the outset, this cycle prioritizes the creation of a minimal viable product (MVP). This MVP will be a fully operational, end-to-end data generation pipeline, albeit with simplified, placeholder logic for the computationally demanding exploration and sampling stages. The core philosophy is to validate the overall architecture, data flow, and user interaction model early, mitigating the risk of discovering fundamental design flaws late in the development process. By focusing on the non-computational infrastructure first, we can ensure that the system is robust, user-friendly, and built on a maintainable codebase.

By the conclusion of this cycle, a developer or user will be able to perform a complete workflow: they can define a simple material system (specifically, a binary alloy) in a well-structured YAML configuration file, invoke the application through a clean and responsive command-line interface, and in return, receive a valid, queryable ASE (Atomic Simulation Environment) database as output. Successfully achieving this milestone will provide tangible proof of the viability of the core architecture and the chosen technology stack. The key deliverables for this cycle are comprehensive and foundational. They include the creation of the project's complete directory structure, a robust set of Pydantic-based data models that enforce a strict configuration schema, a functional and user-friendly CLI, a reliable and abstracted database interaction layer, the central pipeline orchestrator that manages the workflow, and a basic but functional structure generator for alloys. The more advanced components, such as the Molecular Dynamics exploration engine and the intelligent Farthest Point Sampling algorithm, will be implemented as simple "pass-through" placeholders. These placeholders will have the same interface as the final components but will perform trivial operations, allowing the full pipeline to be connected and tested without introducing premature complexity. This strategic, incremental approach is crucial for managing the complexity of the project and ensuring that we can validate our core design choices early, providing a stable and reliable base for subsequent work. A rigorous and comprehensive suite of unit tests for each new component is a critical part of this cycle's definition of done, ensuring that the foundation we build is not only functional but also exceptionally reliable and maintainable from day one.

## 2. System Architecture

The architecture for Cycle 1 is meticulously planned to create the core modules and establish their interactions, laying the groundwork for all future development. The file structure is designed to be modular and scalable, establishing clear boundaries between the different concerns of the application, such as configuration management, data generation, data processing, and storage. This separation is crucial for long-term maintainability and for enabling parallel development of different components in the future. The blueprint below details the specific files to be created and their roles within the system. Each component is designed to have high cohesion—doing one thing well—and low coupling, interacting with other components only through well-defined interfaces. This adherence to established software engineering principles is vital for building a complex scientific application that is both powerful and easy to understand. The architecture explicitly includes not just the application code but also the corresponding test files, reinforcing the project's commitment to a test-driven development methodology where every piece of logic is verifiable and validated from its inception.

**File Structure for Cycle 1:**

The following ASCII tree shows the complete file and directory structure for the project. Files to be created or modified in this cycle are marked in **bold** to provide a clear checklist of the required components. This structure is designed to be self-documenting, with package and module names clearly indicating their purpose.

```
.
├── pyproject.toml
├── src
│   └── mlip_autopipec
│       ├── **__init__.py**
│       ├── **cli.py**
│       ├── config
│       │   ├── **__init__.py**
│       │   └── **models.py**
│       ├── core
│       │   ├── **__init__.py**
│       │   └── **orchestrator.py**
│       ├── exploration
│       │   ├── **__init__.py**
│       │   └── **engine.py**
│       ├── generators
│       │   ├── **__init__.py**
│       │   ├── **base.py**
│       │   └── **alloy.py**
│       ├── sampling
│       │   ├── **__init__.py**
│       │   └── **samplers.py**
│       └── storage
│           ├── **__init__.py**
│           └── **database.py**
└── tests
    ├── **test_cli.py**
    ├── **test_config.py**
    ├── **test_orchestrator.py**
    ├── **test_generator.py**
    └── **test_storage.py**

```

**Component Blueprint:**

*   **`src/mlip_autopipec/`**: This is the root directory for the installable Python package. By organizing the code as a proper package, we gain the benefits of standard Python tooling for dependency management, distribution, and testing.
*   **`config/models.py`**: This is arguably the most critical file in this cycle. It will contain all the Pydantic `BaseModel` classes that define the strict, hierarchical structure of the user-facing YAML configuration files. This schema-first foundation of the project ensures that all configuration data is validated and type-safe before being used by the application, which prevents a large class of common runtime errors. It will define nested models for `SystemConfig`, `ExplorationConfig`, `SamplingConfig`, and the top-level `FullConfig`.
*   **`cli.py`**: This will be the main entry point for the user. It will be built using the `Typer` library for its simplicity and powerful features. It will contain a `run` command that takes the path to a configuration file as an argument. Its role is to parse the configuration, instantiate the main `PipelineRunner`, execute the pipeline, and provide clear, user-friendly feedback to the console using the `Rich` library for formatted output.
*   **`storage/database.py`**: This module will contain the `AseDBWrapper` class. This class will abstract all interactions with the ASE database (which is a SQLite file). It will handle connection management, transaction control, and provide a simple, high-level `write_structures` method. This isolates all database-specific code from the main application logic, making it easy to change the storage backend in the future if needed.
*   **`generators/base.py`**: This file will define the `BaseStructureGenerator`, an abstract base class using Python's `abc` module. This class establishes the contract that all generator classes must adhere to, which is a single `generate` method that returns a list of `ase.Atoms` objects. This is a key element of the extensible, plugin-style architecture.
*   **`generators/alloy.py`**: This will be the first concrete implementation of a generator. The `AlloyGenerator` will be responsible for creating initial atomic structures for alloy systems based on user specifications (e.g., lattice type, element composition, supercell size). It will include basic physical validation to ensure atoms are not placed too closely together, ensuring the quality of the initial data.
*   **`core/orchestrator.py`**: This module houses the `PipelineRunner` class, the central nervous system of the application. It will be initialised with a `FullConfig` object. Its `run` method will orchestrate the entire workflow by instantiating and calling the Generator, Explorer, Sampler, and Storage components in the correct sequence. It is responsible for managing the flow of data between these stages.
*   **`exploration/engine.py`**: In this cycle, this file will contain a placeholder `MDExplorer` class. Its `run_md` method will have the correct signature but will simply log a message and return the input structures unmodified. This acts as a crucial placeholder in the pipeline, allowing the full end-to-end workflow to be tested without the complexity of actual MD simulations.
*   **`sampling/samplers.py`**: Similarly, this file will contain a basic `RandomSampler` class. It will implement a simple algorithm to select a random subset of the structures passed to it, fulfilling the pipeline's structural requirements for this stage and allowing the final database to be populated correctly.

## 3. Design Architecture

The design of MLIP-AutoPipe is strictly schema-first, with Pydantic models serving as the unambiguous, machine-readable definition of all data structures. This approach is a cornerstone of the project's commitment to robustness and maintainability. It ensures data integrity throughout the pipeline and provides clear, self-documenting contracts between the different parts of the system, which is invaluable in a complex scientific application. By defining the data structures first, the logic of the application can be built with a clear understanding of the inputs and outputs of each component. This also facilitates a test-driven development approach, as tests can be written against the well-defined data schemas even before the business logic is fully implemented.

**Pydantic Schema Design:**

*   **`config.models.FullConfig`**: This is the root of the configuration hierarchy. It acts as a container for all other configuration objects, providing a single, coherent view of the entire pipeline configuration. Its fields will be:
    *   `system: SystemConfig`: Defines the physical system to be generated. This is a nested Pydantic model.
    *   `exploration: ExplorationConfig`: Defines the parameters for the exploration stage. This is also a nested model.
    *   `sampling: SamplingConfig`: Defines the parameters for the sampling stage.
    *   `project_name: str`: A simple string for naming the project or run, which will be used to name the output database file.

*   **`config.models.SystemConfig`**: This model will describe the material itself in detail. It is the most critical part of the user-facing configuration.
    *   `elements: list[str]`: A list of chemical symbols (e.g., `['Fe', 'Pt']`). A validator will ensure this list is not empty.
    *   `composition: dict[str, float]`: A mapping from element to its fractional composition (e.g., `{'Fe': 0.5, 'Pt': 0.5}`).
    *   `lattice: str`: The crystal lattice, restricted to a specific set of common values (e.g., `'fcc'`, `'bcc'`, `'hcp'`) using `typing.Literal` to provide immediate feedback on invalid inputs.
    *   `num_structures: int`: The number of initial seed structures to generate. This will be validated to be a positive integer.
    *   `validator`: A Pydantic `@model_validator` will be included to perform cross-field validation. This validator will ensure that the sum of the compositions is sufficiently close to 1.0 (within a small tolerance) and, critically, that the set of elements provided in the `composition` dictionary exactly matches the set of elements in the `elements` list. This prevents inconsistent and ambiguous user input.

*   **`config.models.ExplorationConfig`**: This model will contain settings for the (placeholder) exploration stage. Even though the logic is simple in this cycle, the full data model is defined to ensure forward compatibility.
    *   `temperature: float`: The simulation temperature in Kelvin. A `Field(gt=0)` validator from Pydantic will be used to ensure this is always a positive value.

*   **`config.models.SamplingConfig`**: This model will contain settings for the sampling stage.
    *   `method: typing.Literal['random']`: For Cycle 1, the only allowed sampling method is 'random'. Using `Literal` enforces this constraint at the validation level.
    *   `fraction: float`: The fraction of structures to select from the exploration stage. A `Field(gt=0, le=1)` validator will ensure this is a valid fraction between 0 and 1.

**Data Consumers and Producers:**

The flow of data through the system is well-defined, with each component acting as either a producer or a consumer of data.

*   **Producer (`cli.py`, Hydra):** The user's YAML file is the initial raw data. It is consumed by the CLI, which uses Hydra and Pydantic to parse and validate it, producing a rich, typed `FullConfig` object.
*   **Consumer (`PipelineRunner`):** The `PipelineRunner` is the primary consumer of the `FullConfig` object. It uses this object to configure all of its child components.
*   **Producer (`AlloyGenerator`):** The `AlloyGenerator` consumes the `SystemConfig` part of the configuration and produces a `list[ase.Atoms]` object, which represents the initial, physically-validated atomic structures.
*   **Producer (`RandomSampler`):** The `RandomSampler` consumes the `list[ase.Atoms]` object from the previous stage and the `SamplingConfig`, and produces a new, smaller, randomly-selected `list[ase.Atoms]` object.
*   **Consumer (`AseDBWrapper`):** The database wrapper is the final consumer in the pipeline. It consumes the final `list[ase.Atoms]` object and writes its contents to persistent storage, producing the final output file of the application.

**Extensibility:**

The use of an abstract base class (`BaseStructureGenerator`) is a key design choice that explicitly enables extensibility. This adheres to the open/closed principle of software design. In the future, new generation methods (e.g., for ionic crystals, surfaces, or interfaces) can be added to the system simply by creating a new class that inherits from `BaseStructureGenerator` and implements the `generate` method. The `PipelineRunner`, through the use of a factory pattern, will be able to instantiate and use these new generators without any changes to its own code. This makes the system easy to extend with new scientific capabilities in the future.

## 4. Implementation Approach

The implementation of Cycle 1 will proceed in a logical, step-by-step manner, building the system from the most foundational components (the data structures) up to the user-facing interface. This incremental approach ensures that each layer is built on a solid, tested foundation. The process is designed to be systematic and to produce a high-quality, maintainable codebase from the outset. Each step will be accompanied by the development of corresponding unit tests, following a Test-Driven Development (TDD) philosophy as closely as possible.

1.  **Project Structure Setup:** The very first step is to create the complete directory and file structure as outlined in the System Architecture section. All necessary directories and empty `__init__.py` files will be created to ensure Python treats the directories as packages. The `pyproject.toml` file will be updated to include the `src/mlip_autopipec` path in the `[tool.hatch.build.targets.wheel].packages` section. This crucial step makes the project an installable package, which is essential for `pytest` to correctly discover and run the tests and for standard dependency management tools to work as expected.
2.  **Pydantic Models:** The implementation will begin with the most fundamental part of the application: the data contracts. The `config/models.py` file will be fully implemented. All Pydantic models (`FullConfig`, `SystemConfig`, `ExplorationConfig`, `SamplingConfig`) will be defined with their fields, strict type hints, and all the validation logic (e.g., using `Field` and `@model_validator`). This schema-first approach ensures that the data contracts are clear and enforced before any business logic is written. Corresponding unit tests will be written in `tests/test_config.py` to validate this logic.
3.  **Database Wrapper:** Next, the `storage/database.py` module will be implemented. The `AseDBWrapper` class will be created. It will be designed to be a context manager, using `__enter__` and `__exit__` methods to handle the database connection and cursor automatically. This robust pattern ensures that the database connection is always closed properly, even if errors occur. It will have a single public method, `write_structures`, that takes a list of `ase.Atoms` objects. Tests in `tests/test_storage.py` will verify its functionality against a temporary database.
4.  **Generator Implementation:** The `generators/base.py` file will be created with the `BaseStructureGenerator` abstract class. Then, the `generators/alloy.py` file will be implemented. The `AlloyGenerator` will use functions from the `ase.build` module to create primitive cells and will then contain the logic to apply random element substitutions to achieve the target composition. It will also include a crucial call to a physical validation utility to ensure inter-atomic distances are not unrealistically small. Tests for this will be written in `tests/test_generator.py`.
5.  **Placeholder Components:** The placeholder components for the deferred exploration and sampling stages will be created. The `exploration/engine.py` file will contain the `MDExplorer` class with its pass-through `run_md` method. Similarly, the `sampling/samplers.py` file will contain the `RandomSampler` class with its simple random selection logic. These components, though simple, are vital for ensuring the pipeline can be connected end-to-end.
6.  **Orchestrator Logic:** With all the components (both real and placeholder) and data structures defined, the `core/orchestrator.py` module will be implemented. The `PipelineRunner`'s `run` method will be written to instantiate and call each component in the correct, hardcoded sequence: it will first call the generator, then pass the result to the explorer, then to the sampler, and finally, it will use the `AseDBWrapper` context manager to write the final list of structures to the database. Unit tests with mocks in `tests/test_orchestrator.py` will verify this orchestration logic.
7.  **CLI Entry Point:** Finally, with all the backend logic in place, the `cli.py` module will be created. A Typer application will be defined, with a single `run` command. This command will handle loading the YAML file using a helper function that integrates Hydra and Pydantic, creating the `FullConfig` object, instantiating the `PipelineRunner` with it, and then calling its `run` method. The `Rich` library will be used to print user-friendly status messages to the console, providing a good user experience. The CLI's functionality will be tested in `tests/test_cli.py`.

## 5. Test Strategy

The test strategy for Cycle 1 is heavily focused on comprehensive **unit testing** to ensure that each foundational component is correct, robust, and reliable before they are integrated. This approach is essential for building a stable system, as it allows for the precise isolation and verification of individual pieces of logic. By ensuring each component works perfectly on its own, we significantly simplify the process of debugging issues that may arise during integration. All tests will be written using the `pytest` framework and will be placed in the `tests/` directory, mirroring the structure of the application code.

**Unit Testing Approach:**

*   **`test_config.py`**:
    *   **Goal:** To rigorously verify the validation logic within the Pydantic models, as this is the first line of defense against invalid user input.
    *   **Test Cases:**
        1.  **Valid Configuration:** A test will be created to load a known-valid YAML configuration file and assert that the resulting `FullConfig` object is created successfully with all the correct attribute values and nested objects.
        2.  **Invalid Composition Sum:** A test will use `pytest.raises(ValidationError)` to ensure that a configuration where the chemical compositions do not sum to 1.0 (e.g., `{'Fe': 0.6, 'Pt': 0.5}`) correctly raises a Pydantic validation error. The error message will also be inspected to ensure it is informative.
        3.  **Mismatched Elements:** A test will be written to check that a validation error is raised if the elements in the `composition` dictionary do not match the `elements` list.
        4.  **Invalid Numeric Values:** Another test with `pytest.raises` will check that a negative temperature in the `ExplorationConfig` or a sampling fraction greater than 1.0 are correctly rejected by the `Field` validators.
        5.  **Incorrect Data Types:** A test will verify that providing incorrect data types in the configuration file (e.g., a string for `num_structures`) raises a clear validation error.

*   **`test_storage.py`**:
    *   **Goal:** To verify the `AseDBWrapper`'s database interactions in an isolated environment.
    *   **Test Cases:**
        1.  **Database Creation and Closing:** A test will use the `AseDBWrapper` as a context manager with a temporary file path (provided by the `tmp_path` fixture in pytest). It will assert that the database file is created upon entering the context and that the connection is closed upon exiting.
        2.  **Write and Read Integrity:** A test will create a simple, known `ase.Atoms` object. It will then write this object to a temporary database using the wrapper. Finally, it will use ASE's standard `ase.db.connect` function to read the object back and assert that the retrieved object is identical to the one that was written, verifying the integrity of the serialization process.

*   **`test_generator.py`**:
    *   **Goal:** To verify that the `AlloyGenerator` produces a correct and physically plausible set of atomic structures.
    *   **Test Cases:**
        1.  **Correct Count and Composition:** A test will run the generator for a simple binary alloy (e.g., FePt). It will then assert that the number of returned `ase.Atoms` objects exactly matches the `num_structures` parameter from the configuration. It will also iterate through each generated structure and assert that the number of atoms of each element correctly matches the specified composition.
        2.  **Physical Validity Check:** A test will assert that for every structure generated by the `AlloyGenerator`, the minimum distance between any two atoms is greater than a reasonable physical cutoff (e.g., 1.0 Angstrom), ensuring the prevention of atomic overlaps.

*   **`test_cli.py`**:
    *   **Goal:** To verify that the command-line interface behaves correctly from a user's perspective.
    *   **Test Cases:**
        1.  **Successful Invocation:** The test will use `click.testing.CliRunner` to invoke the `run` command with a valid configuration file. The `PipelineRunner` will be mocked using `unittest.mock.patch` to prevent it from actually running the full pipeline. The test will assert that the CLI returns an exit code of 0 (indicating success) and that the `PipelineRunner.run` method was called exactly once.
        2.  **File Not Found Error:** The test will invoke the CLI with a path to a non-existent configuration file. It will assert that the application exits with a non-zero status code and that a user-friendly "File not found" error message is printed to the standard output.

**Integration Testing Approach:**

No formal, automated integration tests will be created during Cycle 1. The complexity of setting up and running even a simplified version of the full pipeline is non-trivial, and the primary focus of this cycle is on ensuring the correctness of the individual units. The primary integration test will be the "manual" act of running the CLI with a sample configuration file at the very end of the development cycle. This will serve as the final smoke test to confirm that all the individually unit-tested components can be successfully instantiated and called in sequence without raising unexpected errors. True integration tests, which require the careful management of a computational environment and external potentials, are deferred to Cycle 2 when the computational components are fully implemented.
