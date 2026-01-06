# Specification: CYCLE 01 - Core Pipeline and Foundation

## 1. Summary

This document outlines the detailed technical specifications for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to establish a robust, reliable, and functional foundation upon which the more complex features of the application will be built in subsequent cycles. This foundational work will involve creating the core components necessary for a minimal, yet fully operational, end-to-end data generation workflow. The scope of this cycle is intentionally and strategically limited to the essential features required to prove the viability of the core architecture before introducing the computationally intensive scientific algorithms. By focusing on a "thin slice" of functionality, we can validate our design choices, establish a solid testing framework, and deliver a tangible asset early in the development process.

Specifically, this cycle will deliver a command-line interface (CLI) application that serves as the primary user entry point. This CLI will be capable of parsing a user-defined configuration file in YAML format, validating its contents against a strict schema, generating a set of initial atomic structures based on that configuration, and storing these structures persistently in a local database. The key deliverables include a well-defined project structure that promotes modularity and scalability, a comprehensive set of Pydantic models that act as the single source of truth for all configuration data, a dedicated database wrapper to abstract away database-specific logic and ensure safe and reliable data persistence, and a basic but functional implementation of the structure generation and storage engines.

The focus of this cycle is squarely on engineering quality, clean architecture, and the establishment of a thorough testing culture. To this end, we will implement a basic `AlloyGenerator` as the initial concrete implementation of a structure creation module. This will serve as a proof-of-concept for the generator factory pattern, which is a key architectural choice that will allow for the seamless addition of more complex generation logic for different material types (e.g., ionic, covalent) in later cycles. The complex and computationally demanding `Exploration` and `Sampling` stages, which will involve running molecular dynamics simulations and applying sophisticated data selection algorithms, are explicitly out of scope for this cycle. They will be implemented in Cycle 2, where they can be built upon the stable and well-tested foundation laid here. By the end of this cycle, the development team will have a runnable application that demonstrates the core data flow, from configuration input to database output. This will serve as a critical milestone, validating the architectural design and providing a solid, stable platform for the advanced scientific features to come.

## 2. System Architecture

The architecture for Cycle 1 focuses on creating the foundational modules and defining their interactions. The system will be composed of a CLI entry point, a configuration loading module, a central pipeline orchestrator, a generator module with a factory pattern, a storage engine, and a database wrapper. This design emphasizes a strong separation of concerns, ensuring that each component has a single, well-defined responsibility.

**File Structure (Cycle 1 Focus):**
The files to be created or modified in this cycle are marked in bold, indicating that they are the primary deliverables for this phase of the project.

```
.
├── **pyproject.toml**
└── src
    └── mlip_autopipec
        ├── **__init__.py**
        ├── **cli.py**                # **[CREATE]** Main CLI entry point using Typer.
        ├── **config.py**             # **[CREATE]** Pydantic models defining the configuration schema.
        ├── **database.py**           # **[CREATE]** ASE DB wrapper for safe database operations.
        │
        ├── pipeline
        │   ├── **__init__.py**
        │   └── **runner.py**         # **[CREATE]** PipelineRunner orchestrator for the workflow.
        │
        ├── generators
        │   ├── **__init__.py**
        │   ├── **base.py**           # **[CREATE]** Abstract BaseGenerator class defining the interface.
        │   ├── **alloy.py**          # **[CREATE]** Concrete implementation for alloy structures.
        │   └── **factory.py**        # **[CREATE]** Factory to select the appropriate generator.
        │
        ├── exploration             # (Directory created as a placeholder, no implementation in C1)
        │   └── __init__.py
        │
        ├── sampling                # (Directory created as a placeholder, no implementation in C1)
        │   └── __init__.py
        │
        └── storage
            ├── **__init__.py**
            └── **engine.py**         # **[CREATE]** StorageEngine to handle data archival.
```

**Component Interaction Blueprint:**

1.  **`cli.py` (Typer App):**
    *   This module serves as the sole user-facing entry point for the application in this cycle.
    *   It will be implemented using the `Typer` library to create a clean, modern command-line interface.
    *   It will define a main command, for example `run-pipeline`, that accepts a single mandatory argument: `--config-path`, which is the path to the user's YAML configuration file.
    *   Its primary responsibility is to orchestrate the application's startup and shutdown. It will parse the command-line arguments, call the configuration loading function, initialize the `PipelineRunner`, and invoke its `run()` method.
    *   It will implement robust, top-level error handling. Any exception that propagates up from the lower layers of the application will be caught here, and a user-friendly error message will be printed to the console using the `rich` library for better formatting. It will also be responsible for setting the exit code of the application appropriately (0 for success, non-zero for failure).
    *   **Code Blueprint:**
        ```python
        import typer
        from rich.console import Console
        from .pipeline.runner import PipelineRunner
        from .config import load_config, ConfigurationError

        app = typer.Typer()
        console = Console()

        @app.command()
        def run_pipeline(config_path: Path = typer.Option(..., "--config-path", help="Path to the YAML configuration file.")):
            """Runs the MLIP-AutoPipe generation pipeline."""
            try:
                console.print("[bold green]Starting pipeline...[/bold green]")
                config = load_config(config_path)
                runner = PipelineRunner(config)
                runner.run()
                console.print("[bold green]Pipeline completed successfully.[/bold green]")
            except ConfigurationError as e:
                console.print(f"[bold red]Configuration Error:[/bold red] {e}")
                raise typer.Exit(code=1)
            except Exception as e:
                console.print(f"[bold red]An unexpected error occurred:[/bold red] {e}")
                raise typer.Exit(code=1)
        ```

2.  **`pipeline/runner.py` (PipelineRunner):**
    *   This class is the central orchestrator of the entire pipeline. It contains the core business logic that dictates the sequence of operations.
    *   Its constructor will accept the fully validated Pydantic configuration object, ensuring that it only ever operates on valid parameters. This is an example of dependency injection.
    *   The `run()` method will implement the simplified workflow for Cycle 1: **Generation -> Storage**.
    *   It will instantiate the `GeneratorFactory` and use it to get the correct generator instance based on the configuration. It will then call the generator's `generate()` method to create the initial structures.
    *   Finally, it will pass the generated structures to the `StorageEngine` to be saved to the database.
    *   **Code Blueprint:**
        ```python
        from ..config import FullConfig
        from ..generators.factory import GeneratorFactory
        from ..storage.engine import StorageEngine
        from ..database import AseDBWrapper

        class PipelineRunner:
            def __init__(self, config: FullConfig):
                self.config = config
                self.db = AseDBWrapper(path=self.config.database_path)
                self.generator = GeneratorFactory.get_generator(self.config.system)
                self.storage_engine = StorageEngine(self.db)

            def run(self):
                # Stage 1: Generation
                print("Generating initial structures...")
                initial_structures = self.generator.generate()
                print(f"Generated {len(initial_structures)} structures.")

                # Stage 2: Storage
                print("Saving structures to database...")
                self.storage_engine.save(initial_structures, group_name="initial_structures")
                print("Structures saved.")
        ```

## 3. Design Architecture

This cycle establishes a schema-first design paradigm, where the Pydantic models defined in `config.py` serve as the canonical, unambiguous, and single source of truth for all of the application's data structures and configuration parameters. This approach is fundamental to the project's goal of creating a robust and maintainable system, as it ensures that all data flowing through the application is rigorously validated against a strict, well-defined contract. This catches a large class of potential errors at the earliest possible stage, before any significant computation has been performed.

**`config.py` - Pydantic Schemas:**
This file is the cornerstone of the system's design. It provides a declarative and easily readable definition of the expected structure and constraints of the user's YAML configuration file.

*   **`BaseModel` Configuration:** All models will inherit from `pydantic.BaseModel`. They will be strictly configured with `model_config = ConfigDict(extra="forbid")`. This is a critical design choice that prevents users from passing unknown or misspelled configuration keys. If the user provides a key that is not defined in the schema, Pydantic will immediately raise a `ValidationError`, providing clear and actionable feedback. This simple configuration choice eliminates a common and frustrating source of user error.
*   **`SystemConfig`:** This model will define the schema for the `system` section of the YAML file.
    *   `generator_type`: A string field that will be constrained to an `Enum` (`"alloy"`, `"ionic"`, etc.) in a future iteration. For Cycle 1, a simple string is sufficient, and the factory will validate its value. This ensures only supported generator types can be requested.
    *   `elements`: A `list[str]`, which will be validated to contain valid chemical symbols from the periodic table using a custom validator.
    *   `composition`: A `dict[str, float]`. The keys of this dictionary must be a subset of the `elements` list, and the values must sum to 1.0. A custom `@model_validator` will be implemented to enforce this crucial physical constraint, ensuring that the user provides a chemically sensible composition.
    *   `num_structures`: An integer field constrained to be positive using `pydantic.PositiveInt`. This prevents nonsensical requests for zero or a negative number of structures.
*   **`FullConfig`:** This is the top-level model that composes all other configuration models (e.g., `system: SystemConfig`). It provides a single, unified object that represents the entire configuration state of the application. It will also define global parameters, such as `database_path: Path`, which will be automatically cast to a `pathlib.Path` object for robust and platform-agnostic path manipulation.

**Producers and Consumers of Data:**
*   **Producer:** The primary producer of the configuration data is the user, who authors a YAML file. The `pyyaml` library will act as the parser, converting this text file into a Python dictionary. Pydantic's model validation (`FullConfig.model_validate`) then consumes this dictionary and produces the validated `FullConfig` object.
*   **Consumers:** Every major component within the application will be a consumer of the `FullConfig` object or its sub-models. For example, the `PipelineRunner` consumes the entire object, while the `GeneratorFactory` specifically consumes the `SystemConfig` sub-model. This pattern of dependency injection, where components are given the configuration they need, makes the system highly modular and easy to test, as each component explicitly declares its dependencies.

**Extensibility:**
The design is explicitly built for future extensibility. The `GeneratorFactory` uses the `generator_type` string from the configuration to decide which concrete generator class to instantiate. To add a new generator in the future (e.g., an `IonicGenerator` in Cycle 2), a developer would simply need to create the new generator class, ensuring it inherits from the `BaseGenerator` abstract class, and then update the factory's simple `if/elif/else` logic to recognize the new `generator_type`. This powerful design pattern decouples the high-level pipeline runner from the low-level implementation details of the various generators, making the system easy to extend without requiring modifications to the core pipeline logic.

## 4. Implementation Approach

The implementation will proceed in a logical, dependency-first order. We will start by building the core data structures and foundational components, and then progressively build the higher-level logic on top of them. This ensures that at each step, we are building upon a solid, tested foundation.

1.  **Project Scaffolding:** Create the complete directory structure as defined in the System Architecture section. This includes creating all sub-packages (`pipeline`, `generators`, `storage`, etc.) and initializing them with `__init__.py` files to ensure they are recognized as Python packages.
2.  **Pydantic Models (`config.py`):** The very first implementation task is to write the Pydantic models. This is the "schema-first" approach in action. We will define the `SystemConfig` and `FullConfig` models with all their fields, type hints, and validation constraints. This includes implementing the custom `@model_validator` for the composition sum. This step is complete only when the schema is fully defined and can correctly validate both correct and incorrect example configuration dictionaries.
3.  **Database Wrapper (`database.py`):** Next, implement the `AseDBWrapper` class. This class will be a high-level abstraction over the `ase.db` library. It will handle the details of connecting to the database, starting transactions, and writing `Atoms` objects. The goal is to create a simple, safe interface (e.g., `add_structures(atoms_list, group_name)`) so that the rest of the application does not need to know the specific details of the `ase.db` API. This encapsulation is key for maintainability.
4.  **Base Generator (`generators/base.py`):** Define the abstract base class `BaseGenerator` using Python's `abc` module. It will define the public interface for all generators by declaring an abstract method `generate() -> list[Atoms]`. It can also provide concrete implementations of common utility or validation methods (e.g., an `_is_valid_structure` method that checks for atomic overlaps) that can be inherited and reused by all subclasses, reducing code duplication.
5.  **Alloy Generator (`generators/alloy.py`):** Create the concrete `AlloyGenerator` class, ensuring it inherits from `BaseGenerator`. Implement the `generate()` method according to the contract defined by the base class. The initial implementation will use `ase.build.bulk` to create a primitive cell, `ase.build.make_supercell` to expand it, and then custom logic to randomly replace atomic symbols based on the requested composition. It must use the validation methods inherited from its base class to ensure it only outputs physically plausible structures.
6.  **Generator Factory (`generators/factory.py`):** Implement the `GeneratorFactory` class. It will have a single static method, `get_generator(config: SystemConfig) -> BaseGenerator`. This method will contain a simple `if/elif/else` control structure that inspects `config.generator_type` and returns an instance of the corresponding concrete generator class. It will raise an error if the type is not recognized.
7.  **Storage Engine (`storage/engine.py`):** Implement the `StorageEngine` class. Its role is simple in this cycle. Its constructor will take an `AseDBWrapper` instance (dependency injection), and it will have a single public method, `save(atoms_list, group_name)`, which will simply delegate the call to the database wrapper. This provides a nice abstraction layer for the concept of "storage".
8.  **Pipeline Runner (`pipeline/runner.py`):** With all the backend components now in place and unit-tested, we can implement the `PipelineRunner`. It will tie everything together, instantiating the components and calling their methods in the correct sequence as described in the blueprint in the System Architecture section.
9.  **CLI (`cli.py`):** The final implementation step is the user-facing CLI. With all the underlying logic complete, this module becomes a thin layer that is responsible only for parsing command-line arguments and triggering the pipeline.
10. **Testing:** Throughout this entire process, unit tests will be written concurrently with the implementation of each component, following the detailed strategy outlined below. This Test-Driven Development (TDD) -like approach ensures that each component is verified before it is integrated into the larger system.

## 5. Test Strategy

Testing in Cycle 1 is of paramount importance as it establishes the foundation for the quality of the entire application. The strategy is to have comprehensive unit tests for every component and a small set of integration tests to verify that the components work together correctly for the primary use case.

**Unit Testing Approach (Min 300 words):**
Unit tests will be created in the `tests/` directory, with a file structure that mirrors the source code (e.g., `tests/generators/test_alloy.py`). The `pytest` framework will be used as the test runner, and `pytest-mock` will be used for mocking dependencies to ensure tests are isolated.
*   **`test_config.py`:** This will be one of the most extensive and critical test suites. It will rigorously test the Pydantic models, as they are the gateway for all data into the system.
    *   A "happy path" test will load a dictionary representing a perfectly valid configuration and assert that the `FullConfig` model is created successfully without raising any exceptions.
    *   Multiple failure-case tests will be written to verify the custom validators. For example, one test will provide a composition that does not sum to 1.0 and assert that a `pydantic.ValidationError` is raised. Another will provide an element in the composition that is not in the `elements` list.
    *   A test for the `extra="forbid"` configuration will provide a dictionary with an extra, undefined key (e.g., `colour: "red"`) and assert that a `ValidationError` is raised, confirming that we protect against user typos.
    *   Tests for Pydantic's built-in constraints will be included, such as providing a negative value for `num_structures` and asserting that a `ValidationError` is raised.
*   **`test_database.py`:** To test the `AseDBWrapper` in isolation from the filesystem, the `ase.db.connect` function will be mocked using `mocker.patch`. We will create a mock connection object that mimics the real object's interface. The test will then call `AseDBWrapper.add_structures()` and assert that the `connection.write()` method on the mock object was called with the correct `Atoms` objects and any associated metadata. This allows us to verify the wrapper's logic without the overhead and non-determinism of disk I/O.
*   **`test_generators.py`:** The `AlloyGenerator` will be tested to ensure its output is both syntactically and semantically correct. We will initialize it with a simple but deterministic configuration (e.g., generate 5 SiC structures in a 2x2x2 supercell). We will then assert that the `generate()` method returns a list containing exactly 5 `ase.Atoms` objects. We will then perform a detailed inspection of one of these objects to assert that it contains the correct total number of atoms and that the ratio of Si to C atoms is correct. We will also write a specific test for the overlap validation logic by manually creating an `Atoms` object with two atoms at the same coordinate and asserting that the generator's internal validation method correctly identifies and flags this invalid structure.

**Integration Testing Approach (Min 300 words):**
While unit tests are essential for verifying components in isolation, integration tests are necessary to ensure that they communicate and work together as a cohesive system. For Cycle 1, we will focus on a single, critical integration test for the main CLI workflow.
*   **`test_cli.py`:** We will use the `typer.testing.CliRunner` utility to perform an end-to-end test of the simplified pipeline, from command-line invocation to database creation.
    *   The test function will be structured using the "Arrange, Act, Assert" pattern.
    *   **Arrange:** The test will first create a temporary directory using `pytest`'s built-in `tmp_path` fixture. This ensures the test is hermetic and does not leave artifacts on the filesystem. Inside this directory, it will programmatically create a simple, valid `config.yaml` file. This configuration will specify creating a small number of structures (e.g., 3) and will set the `database_path` to a file within the same temporary directory (e.g., `tmp_path / "test.db"`).
    *   **Act:** The `CliRunner.invoke()` method will be called, passing in the arguments to run the pipeline, including the path to the temporary config file.
    *   **Assert:** The test will perform several assertions to verify the outcome. First, it will assert that the CLI command exited with a success code of 0 and that the output printed to the console contains a success message. Second, it will check the filesystem to assert that the database file was actually created at the specified path. Finally, and most importantly, it will use the real `ase.db.connect` to connect to this newly created database. It will query the database and assert that it contains exactly the number of structures that were specified in the configuration file. This final check verifies that the entire data flow—from CLI parsing, through configuration, generation, and storage—is working correctly.
    *   A second integration test will be created to verify the failure path. It will use an invalid config file and assert that the CLI exits with a non-zero exit code and prints an informative error message.
