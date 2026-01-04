# CYCLE01/SPEC.md

## 1. Summary

Cycle 1 marks the foundational stage of the MLIP-AutoPipe project. The primary objective of this cycle is to deliver a robust, end-to-end, command-line-driven workflow for generating and storing initial atomic structures. This initial version will not include the advanced simulation or sampling capabilities; instead, it focuses on establishing the core architecture, data models, and essential components that all subsequent features will rely upon. The key deliverables for this cycle are a fully functional Command-Line Interface (CLI) built with Typer for a modern user experience, a comprehensive set of Pydantic models for type-safe configuration, a factory for creating various structure generators, concrete implementations for Alloy and Ionic structure generators, and a dedicated wrapper for managing the ASE database. The design philosophy for this cycle is to build a "walking skeleton" of the application. This approach, common in agile development, involves creating a minimal yet functional end-to-end implementation of the workflow. In this case, it means a user can define a material in a configuration file, run a single command, and see a valid database of atomic structures as the output.

This methodology is crucial for several reasons. Firstly, it forces the early integration of all major architectural layers—from the user-facing CLI and configuration parsing down to the business logic of the generators and the persistence layer of the database. This early integration helps to identify and resolve potential architectural flaws or incorrect assumptions about how components will interact, which is far less costly than discovering such issues later in the development process. Secondly, it provides a tangible and testable artifact early on. This allows for immediate feedback and validation, ensuring that the project is on the right track. By the end of this cycle, we will have a solid, command-line-driven tool that, while limited in features, is complete in its end-to-end execution. This provides immense value, as it can already be used for basic structure generation tasks. More importantly, it provides a stable and well-tested foundation upon which the more complex and scientifically innovative features of Cycle 2—such as the MD/MC exploration engine—can be confidently built. The successful completion of Cycle 1 will result in a reliable, extensible, and immediately useful piece of software that perfectly sets the stage for the advanced capabilities to come.

## 2. System Architecture

The architecture for Cycle 1 is intentionally focused and streamlined. It involves creating the essential file structure and components for the core pipeline, specifically the Generation and Storage stages. The modules for the subsequent Exploration and Sampling stages will be created as placeholder files with minimal content, clearly delineating where future functionality will be integrated in Cycle 2. This approach ensures a clean and well-organized codebase from the outset.

**File Structure for Cycle 1:**

The ASCII tree below provides a detailed blueprint of the files to be created or modified during this cycle. The files marked in **bold** are the primary targets for implementation. This structure adheres to modern Python packaging standards, with application code residing in `src/mlip_autopipec` and tests in the `tests` directory. This separation is crucial for maintainability and for building a distributable package.

```
src/mlip_autopipec/
├── **__init__.py**              # Makes the directory a Python package.
├── **cli.py**                     # **Implementation Target**: Typer-based CLI application. Will contain the `run-pipeline` command.
├── **config.py**                  # **Implementation Target**: All Pydantic configuration models (e.g., `SystemConfig`, `FullConfig`).
├── **database.py**                # **Implementation Target**: The `AseDBWrapper` class for all database operations.
├── **pipeline.py**                # **Implementation Target**: The initial `PipelineRunner` class, orchestrating only generation and storage.
├── **generators/**                # Package for all structure generation logic.
│   ├── **__init__.py**
│   ├── **base.py**                # **Implementation Target**: Defines the `BaseStructureGenerator` Abstract Base Class.
│   ├── **alloy.py**               # **Implementation Target**: Concrete implementation for generating alloy structures.
│   └── **ionic.py**               # **Implementation Target**: Concrete implementation for generating ionic crystal structures.
└── explorers/                     # Placeholder for Cycle 2.
    ├── __init__.py                # Will be created as an empty file.
    └── md_engine.py               # Will be created as an empty file.
└── samplers/                      # Placeholder for Cycle 2.
    ├── __init__.py                # Will be created as an empty file.
    ├── base.py                    # Will be created as an empty file.
    └── fps.py                     # Will be created as an empty file.
tests/
├── **conftest.py**                # **Implementation Target**: Shared pytest fixtures, such as a temporary database fixture.
├── **test_cli.py**                # **Implementation Target**: Integration tests for the Typer CLI application.
├── **test_config.py**             # **Implementation Target**: Unit tests for the Pydantic configuration models.
├── **test_database.py**           # **Implementation Target**: Unit tests for the `AseDBWrapper`.
└── **test_generators.py**         # **Implementation Target**: Unit tests for the `AlloyGenerator` and `IonicGenerator`.
```

The core of this architecture is its modularity. The `generators` are developed as a distinct package, completely decoupled from the `database` and `pipeline` logic. This is enforced through the use of the `BaseStructureGenerator` ABC, which acts as a formal contract. The `PipelineRunner` in `pipeline.py` will be designed to work with any object that fulfills this contract, meaning it will be able to run any new generator added in the future without modification. Similarly, the `AseDBWrapper` encapsulates all database-specific code. If, in the future, we decided to switch from SQLite to a different database backend (like PostgreSQL), only the `database.py` file would need to be changed; the rest of the application, which interacts only with the wrapper, would remain untouched. This strict separation of concerns is fundamental to creating a maintainable and extensible system. The placeholder files for `explorers` and `samplers` serve as a clear architectural roadmap, showing exactly where the new functionality of Cycle 2 will be integrated, ensuring a smooth transition between development phases. This level of planning and structure is essential for the long-term health and success of the project.

## 3. Design Architecture

This cycle emphasizes a rigorous schema-first design approach, with Pydantic models serving as the "source of truth" for all data structures within the system. This methodology is central to the project's commitment to robustness and reliability. All components, from the CLI to the database wrapper, will be designed to consume these validated models, ensuring data integrity and consistency throughout the entire pipeline. This prevents a wide class of common runtime errors and makes the system's behavior far more predictable.

**Pydantic-Based Schema Design (`config.py`):**

The configuration will be meticulously structured into a hierarchy of nested Pydantic models. This provides a clear, self-documenting, and type-safe schema for all user-configurable parameters. This is vastly superior to using untyped dictionaries, which are prone to key errors and type mismatches.

*   `SystemConfig`: This model will define the core physical system.
    *   `elements`: A `list[str]` of chemical symbols (e.g., `['Fe', 'Pt']`).
    *   `composition`: A `dict[str, float]` mapping each element to its fractional composition. The model will include a custom `@model_validator` to ensure the sum of the composition values is approximately 1.0, a critical physical constraint.
    *   `crystal_structure`: An `Enum` of allowed crystal structures (e.g., 'fcc', 'bcc', 'hcp'), preventing users from entering invalid strings.
    *   `n_structures`: A `PositiveInt` from Pydantic, which automatically validates that the number is an integer greater than zero.

*   `GeneratorConfig`: A container for generator-specific settings.
    *   `type`: An `Enum` to select the generator type (e.g., 'alloy', 'ionic'). This will be used by the factory to instantiate the correct generator class.
    *   It will also support additional fields for generator-specific parameters, using Pydantic's `Extra.allow` if necessary.

*   `DBConfig`: Contains settings for the database.
    *   `path`: A `Path` object from `pathlib`, ensuring the path is handled correctly across different operating systems. It will have a default value (e.g., `'./structures.db'`).

*   `FullConfig`: The top-level model that aggregates all other configuration models into a single, cohesive object (`system: SystemConfig`, `generator: GeneratorConfig`, `db: DBConfig`). This model will be the primary data object passed between high-level components.

All models will be configured with `model_config = ConfigDict(extra="forbid")` to strictly prevent users from providing unknown or misspelled configuration fields, providing immediate and clear feedback on typos. This rigorous, upfront validation is a cornerstone of the system's design.

**Component Interfaces (Contracts):**

*   **`generators.base.BaseStructureGenerator` (ABC):**
    This Abstract Base Class defines the formal contract for all structure generators, ensuring interoperability.
    ```python
    from abc import ABC, abstractmethod
    from ase import Atoms
    from ..config import SystemConfig

    class BaseStructureGenerator(ABC):
        @abstractmethod
        def generate(self, config: SystemConfig) -> list[Atoms]:
            """Generates a list of Atoms objects based on the system config."""
            pass
    ```
    This design is a textbook example of the Dependency Inversion Principle. The high-level `PipelineRunner` depends on this abstraction (`BaseStructureGenerator`), not on the concrete low-level implementations (`AlloyGenerator`). This decouples the components and is the key to the system's extensibility.

*   **`database.AseDBWrapper`:**
    This class will provide a clean, high-level API for all database operations, completely encapsulating the implementation details of the underlying `ase.db` library.
    ```python
    from ase import Atoms

    class AseDBWrapper:
        def __init__(self, db_path: str): ...
        def write_atoms(self, atoms_list: list[Atoms], metadata: dict) -> None: ...
        def get_natoms(self) -> int: ...
        def get_all_atoms(self) -> list[Atoms]: ...
    ```
    This design isolates the rest of the application from the persistence layer. The `PipelineRunner` does not need to know how the atoms are being stored, only that it can call the `write_atoms` method on the wrapper. This makes the system far easier to test, as the `AseDBWrapper` can be easily replaced with a mock object during unit tests of the pipeline.

## 4. Implementation Approach

The implementation for Cycle 1 will proceed in a logical, bottom-up sequence. This approach ensures that foundational components are built and tested before the components that depend on them are created. This systematic process minimizes integration issues and facilitates a smooth development workflow.

1.  **Implement Pydantic Models (`config.py`):** The first step is to lay the foundation by defining the data structures. The `SystemConfig`, `GeneratorConfig`, `DBConfig`, and `FullConfig` models will be implemented with strict type hints and validation rules as designed in the previous section. This includes writing custom validators, for example, to ensure the chemical compositions sum to 1.0. Crucially, a corresponding test file, `tests/test_config.py`, will be created in parallel. This test file will contain unit tests that verify both valid and invalid configuration scenarios, ensuring the validation logic is correct and robust from the very beginning.

2.  **Implement Database Wrapper (`database.py`):** With the data models in place, the next step is to build the persistence layer. The `AseDBWrapper` class will be created. Its `__init__` method will establish a connection to the SQLite database file specified by the path. The core `write_atoms` method will be implemented to iterate through a list of `ase.Atoms` objects and save each one to the database using the `ase.db` library's API. A corresponding `tests/test_database.py` will be written to unit test this wrapper, using pytest's `tmp_path` fixture to ensure that tests are isolated and do not leave artifacts on the filesystem.

3.  **Implement Generator Abstraction (`generators/base.py`):** The next step is to define the contract for all generators. The `BaseStructureGenerator` Abstract Base Class will be created with its single abstract method, `generate`. This step is simple but architecturally critical, as it enforces the design contract for all subsequent generator implementations.

4.  **Implement Concrete Generators (`generators/alloy.py`, `generators/ionic.py`):** Now, we will create the concrete generator classes, `AlloyGenerator` and `IonicGenerator`, ensuring they inherit from `BaseStructureGenerator`. The `generate` method for the `AlloyGenerator` will contain the logic for creating random alloy structures. This will involve using the `ase.build.bulk` function to create a pristine supercell and then iterating through the atoms to randomly assign chemical symbols based on the specified composition. The `IonicGenerator` will be more complex, requiring logic to ensure charge neutrality. A corresponding `tests/test_generators.py` will be created to unit test each generator, asserting that they produce the correct number of physically plausible structures.

5.  **Implement Initial Pipeline (`pipeline.py`):** With the core components built and tested, it is time to orchestrate them. The `PipelineRunner` class will be created. Its `__init__` method will accept a `FullConfig` object. Its `run` method will implement the simple, two-stage workflow for this cycle. It will first use a simple factory function (or a match-case block) to instantiate the correct generator based on the `config.generator.type` string. It will then call the generator's `generate` method to produce the list of `ase.Atoms` objects. Finally, it will instantiate the `AseDBWrapper` and call `write_atoms` to save the generated structures to the database.

6.  **Implement CLI (`cli.py`):** The final step is to create the user-facing entry point. A Typer application will be created with a single command, `run-pipeline`. This command will use `typer.Argument` to accept the path to a YAML configuration file. The function body will be responsible for loading and parsing the YAML file, then using the dictionary to instantiate the `FullConfig` Pydantic model. This automatically triggers all the validation logic. If validation succeeds, it will instantiate the `PipelineRunner` with the config object and call its `run` method. It will also include robust error handling, using a `try...except` block to catch `FileNotFoundError` or `pydantic.ValidationError` and print a user-friendly error message to the console before exiting gracefully. `tests/test_cli.py` will be written to perform an end-to-end integration test of this command.

## 5. Test Strategy

The testing strategy for Cycle 1 is designed to build a pyramid of tests, with a broad base of fast, isolated unit tests and a narrow peak of slower, comprehensive integration tests. This approach ensures high code coverage and confidence in the system's correctness while keeping the test suite fast and maintainable. The `pytest` framework will be used for all tests.

**Unit Testing Approach (Min 300 words):**
Unit tests are the foundation of our testing strategy. They are designed to be fast, deterministic, and to test each piece of the application in complete isolation. We will use `pytest` for test running and `pytest-mock` for creating mock objects to isolate components from their dependencies.
*   **Configuration Models (`test_config.py`):** Each Pydantic model will be subjected to rigorous testing. We will create test cases for the "happy path," ensuring that a valid dictionary is correctly parsed into a model instance with the correct types and default values. More importantly, we will write numerous tests for failure scenarios. For example, we will test what happens when a required field is missing, when a field has the wrong data type (e.g., a string for `n_structures`), or when a value violates a constraint (e.g., `n_structures=0` or a composition that does not sum to 1.0). For each of these invalid cases, we will use `pytest.raises(pydantic.ValidationError)` to assert that the expected exception is raised. This ensures that our configuration system is robust and provides clear feedback to the user on invalid input.
*   **Generators (`test_generators.py`):** The concrete generator classes (`AlloyGenerator`, `IonicGenerator`) will be tested to confirm their correctness. The tests will instantiate a generator and call its `generate` method with a test `SystemConfig` object. The assertions will be manifold: we will check that the returned list contains the correct number of `ase.Atoms` objects. We will then inspect individual `Atoms` objects, asserting that they contain the correct chemical elements, that the composition is as expected, and that no two atoms are unrealistically close to each other (a common failure mode in structure generation). We will use `pytest.mark.parametrize` to efficiently test each generator with a variety of different configurations (e.g., different elements, compositions, and crystal structures).
*   **Database Wrapper (`test_database.py`):** The `AseDBWrapper` will be tested in isolation from the rest of the application. We will use pytest's built-in `tmp_path` fixture to create a new, empty directory for each test function. The test will initialize the `AseDBWrapper` with a path to a database file inside this temporary directory. We will then create a known list of `ase.Atoms` objects and call `write_atoms`. The core of the test will be to then create a *new* `AseDBWrapper` instance connected to the same file, read the atoms back, and assert that the retrieved structures are identical to the ones we originally wrote. This round-trip testing ensures the integrity and correctness of our data persistence layer.

**Integration Testing Approach (Min 300 words):**
While unit tests verify the components in isolation, integration tests are crucial for verifying that these components work together correctly as a system. For Cycle 1, the primary integration test will focus on the main user workflow: running the pipeline from the command line.
*   **CLI to Database (`test_cli.py`):** This test suite will use the `click.testing.CliRunner`, a utility specifically designed for testing command-line applications built with Click (which Typer is based on). The test function will first programmatically create a temporary directory using the `tmp_path` fixture. Inside this directory, it will create a valid YAML configuration file (e.g., `config.yml`). It will then use the `CliRunner.invoke()` method to execute the `run-pipeline` command, passing the path to the temporary configuration file as an argument. The first set of assertions will be on the result of the invocation itself: we will check that the command's exit code is 0, indicating success, and that there is no unexpected error output in the captured stdout/stderr. The main verification, however, happens after the command has finished. The test will construct the expected path to the output database file (e.g., `tmp_path / 'structures.db'`) and assert that the file exists. It will then use the `AseDBWrapper` to connect to this database and perform assertions on its contents. We will verify that the database contains the correct number of structures, as specified in our temporary `config.yml`. This end-to-end test provides a very high degree of confidence, as it validates the entire application stack, from the parsing of command-line arguments and the loading of the configuration file, through the instantiation of the pipeline and the execution of the generator, to the final, successful storage of the results in the database. We will also write separate integration tests for failure cases, such as providing a path to a non-existent config file, and assert that the CLI prints a user-friendly error message and exits with a non-zero status code.