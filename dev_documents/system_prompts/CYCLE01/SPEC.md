# SPEC.md - CYCLE 01

## 1. Summary

This document provides the detailed technical specifications for the first development cycle of the **MLIP-AutoPipe** project. The primary objective of Cycle 01 is to establish a solid and functional foundation for the application. This cycle focuses on implementing the core architectural components, including the project's scaffolding, the configuration and data modeling framework using Pydantic, a robust database wrapper for data persistence, and a basic, two-stage pipeline capable of generating initial alloy structures and storing them.

The scope of this cycle is intentionally limited to the essentials required to prove the viability of the core concept. We will develop a simplified pipeline that executes only the **Generation** and **Storage** stages. The more complex Exploration and Sampling stages are explicitly deferred to Cycle 02. By the end of this cycle, a developer or user will be able to define a simple alloy system in a YAML configuration file, run the command-line interface (CLI), and have a set of physically-validated seed structures saved into an ASE-compatible SQLite database. This deliverable will serve as a stable baseline upon which the more advanced simulation and sampling features will be built in the subsequent cycle. The emphasis is on code quality, clear interfaces, and robust testing of the foundational components.

Key deliverables for this cycle include:
- A well-defined directory structure and `pyproject.toml`.
- Pydantic models for parsing user configurations.
- An `AseDBWrapper` class that cleanly abstracts database interactions.
- A `BaseGenerator` interface and a concrete `AlloyGenerator` implementation.
- A minimal `WorkflowOrchestrator` service to run the Generation -> Storage pipeline.
- A basic Typer-based CLI to launch the process.
- Comprehensive unit and integration tests for all implemented components.

## 2. System Architecture

The architecture for Cycle 01 is a subset of the final system architecture, focusing on the components necessary for the initial two-stage pipeline. The system will be built within the `src/mlip_autopipec` package.

**File Structure (Cycle 01 Create/Modify):**

```
.
├── **pyproject.toml**              # To be updated with dependencies
└── src
    └── mlip_autopipec
        ├── **__init__.py**
        ├── **cli.py**                # CLI entry point
        ├── **config.py**             # Pydantic configuration models
        ├── **db.py**                 # AseDBWrapper for database logic
        ├── **domain**
        │   ├── **__init__.py**
        │   ├── **models.py**         # (Optional) Internal domain models
        │   └── **exceptions.py**     # Custom exception classes
        ├── **generators**
        │   ├── **__init__.py**
        │   ├── **base.py**           # BaseGenerator ABC
        │   └── **alloy.py**          # AlloyGenerator implementation
        ├── **pipeline.py**           # Simplified PipelineRunner
        └── **services**
            ├── **__init__.py**
            └── **workflow.py**       # WorkflowOrchestrator service
```

**Component Blueprint:**

*   **`cli.py`**: This module will contain the Typer application. It will define a single command, `run`, that accepts two arguments: `--config-path` (a string path to the user's YAML file) and `--db-path` (a string path for the output database). Its sole responsibility is to parse these inputs, instantiate the `WorkflowOrchestrator`, and invoke its `execute` method within a try/except block for graceful error handling.

*   **`config.py`**: This file is central to the project's flexibility. It will define a set of Pydantic `BaseModel` classes that mirror the structure of the required YAML configuration file.
    *   `SystemConfig`: Will contain fields like `elements` (a list of strings, e.g., `['Fe', 'Pt']`), `composition` (a dictionary mapping element to ratio, e.g., `{'Fe': 0.5, 'Pt': 0.5}`), and `num_structures` (an integer).
    *   `GeneratorConfig`: A container for generator-specific settings.
    *   `FullConfig`: The top-level model that nests the other configuration models. It will have a method `from_yaml(filepath: str)` to load and parse a configuration file.

*   **`db.py`**: This module implements the `AseDBWrapper`.
    *   The class will be initialized with the `db_path`.
    *   It will have a primary method, `write_structures(structures: list[ase.Atoms])`. This method will handle the connection to the SQLite database (creating it if it doesn't exist) and iterate through the list of `ase.Atoms` objects, writing each one to the database. It will use a context manager (`with connect(self.db_path) as db:`) to ensure the database connection is properly handled.

*   **`generators/base.py`**: Will define an abstract base class `BaseGenerator` using Python's `abc` module. It will have one abstract method, `generate() -> list[ase.Atoms]`, which all concrete generator classes must implement.

*   **`generators/alloy.py`**: Will contain the `AlloyGenerator` class, which inherits from `BaseGenerator`.
    *   Its constructor will accept the relevant `SystemConfig`.
    *   The `generate` method will contain the logic for creating random alloy structures. This involves creating a base crystal lattice, making a supercell, and then randomly assigning atomic species to the sites according to the specified composition.
    *   It must include a private method `_is_valid(atoms: ase.Atoms) -> bool` that checks for physical constraints, specifically ensuring that no two atoms are closer than a defined threshold. This method will be called for every generated structure before it is returned.

*   **`services/workflow.py`**: This module contains the `WorkflowOrchestrator`.
    *   Its `execute` method will be the main entry point for the business logic.
    *   It will first instantiate the `AlloyGenerator` based on the configuration.
    *   It will then call the generator's `generate` method to get the list of seed structures.
    *   Finally, it will instantiate the `AseDBWrapper` and call `write_structures` to save the generated data. This orchestrates the simple Generation -> Storage flow.

## 3. Design Architecture

The design for Cycle 01 is centered around establishing a robust, schema-first development process using Pydantic. This ensures that all data, whether from user input or internal processing, is strictly validated, preventing a wide class of runtime errors and making the system's behavior predictable and reliable.

**Pydantic Schema Design (`config.py`):**

The configuration schema is the primary contract between the user and the application. It must be designed to be both intuitive for the user and easy to parse and validate for the application.

*   **`FullConfig(BaseModel)`**: The root of our configuration.
    *   `system: SystemConfig`: Defines the physical system to be generated.
    *   `generator: GeneratorConfig`: Namespace for all generator-related settings.
    *   **Invariants**: The `model_config` will be set to `ConfigDict(extra="forbid")` to prevent users from passing unknown configuration keys, which helps catch typos.
    *   **Consumers**: The `WorkflowOrchestrator` is the primary consumer. The `cli.py` module produces it by parsing the YAML file.
    *   **Extensibility**: By nesting configuration into categories (`system`, `generator`), we make it easy to add new top-level categories in the future (e.g., `exploration`, `sampling`) without breaking backward compatibility for existing YAML files.

*   **`SystemConfig(BaseModel)`**:
    *   `elements: list[str]`: A list of elemental symbols. Must contain at least one element.
    *   `composition: dict[str, float]`: Maps elements to their fractional composition.
    *   `num_structures: int = Field(gt=0)`: The number of structures to generate. Must be a positive integer.
    *   **Constraints and Validation**:
        *   A `@model_validator` will be implemented to ensure that the keys in `composition` are a subset of the `elements` list.
        *   Another `@model_validator` will check that the sum of the values in `composition` is approximately equal to 1.0 (within a small tolerance). This enforces a valid physical constraint.

*   **`GeneratorConfig(BaseModel)`**:
    *   `name: Literal["alloy"]`: For Cycle 01, this is a fixed string. In the future, it will be a literal that allows for selecting different generator types. This acts as a factory key.
    *   `min_atomic_distance: float = Field(gt=0)`: The minimum allowed distance between any two atoms, used for physical validation.

This schema-first approach means we define the "shape" of our data before writing the logic that manipulates it. This serves as a precise pre-implementation design document.

**Internal Data Flow:**

The data flow within the application is designed to be a simple, linear progression of strongly-typed objects.
1.  **`str` (filepath)**: The CLI receives a raw string path to the config file.
2.  **`FullConfig` (Pydantic Model)**: `cli.py` passes the filepath to a `FullConfig.from_yaml()` classmethod, which parses and validates the file, producing a rich `FullConfig` object.
3.  **`list[ase.Atoms]`**: This validated `FullConfig` object is passed to the `WorkflowOrchestrator`. The orchestrator uses it to configure and instantiate the `AlloyGenerator`. The generator's `generate()` method is called, which returns a list of `ase.Atoms` objects. The `ase.Atoms` object is a well-defined data structure from an external library, acting as our primary internal domain model for an atomic structure.
4.  **Database Records**: The `WorkflowOrchestrator` passes the `list[ase.Atoms]` to the `AseDBWrapper`'s `write_structures` method, which serializes these objects into the SQLite database.

This typed, sequential flow minimizes the chance of errors and makes the logic easy to follow and debug.

## 4. Implementation Approach

The implementation will proceed in a logical, bottom-up order, starting with the data structures and moving up to the user-facing CLI.

1.  **Project Setup**:
    *   Create the directory structure as outlined in the System Architecture section.
    *   Initialize `pyproject.toml`. Add the necessary dependencies for Cycle 01: `typer`, `pydantic`, `pyyaml`, and `ase`. Also, include development dependencies like `pytest`.

2.  **Pydantic Models (`config.py`)**:
    *   Implement the `SystemConfig`, `GeneratorConfig`, and `FullConfig` Pydantic models.
    *   Add the field validators (`@model_validator`) to `SystemConfig` to ensure the consistency of `elements` and `composition` and to check that the composition sums to 1.0.
    *   Implement the `from_yaml` classmethod on `FullConfig` to handle loading the configuration from a file path.

3.  **Database Wrapper (`db.py`)**:
    *   Create the `AseDBWrapper` class.
    *   The `__init__` method will simply store the `db_path`.
    *   The `write_structures` method will implement the core logic using `ase.db.connect`. The method will look like:
        ```python
        from ase.db import connect

        def write_structures(self, structures: list[ase.Atoms]):
            with connect(self.db_path) as db:
                for atoms in structures:
                    db.write(atoms)
        ```

4.  **Generator Implementation (`generators/`)**:
    *   In `base.py`, define the `BaseGenerator` ABC with the single abstract method `generate`.
    *   In `alloy.py`, implement the `AlloyGenerator`. The `generate` method will use `ase.build.bulk` to create a primitive cell, `ase.build.make_supercell` to expand it, and then iterate through the atoms, randomly changing their symbols based on the composition.
    *   Implement the `_is_valid` method using `atoms.get_all_distances()` and checking if the minimum distance (ignoring the diagonal) is above the threshold from the config.

5.  **Workflow Service (`services/workflow.py`)**:
    *   Create the `WorkflowOrchestrator` class.
    *   The `execute(config: FullConfig, db_path: str)` method will tie everything together. It will:
        1.  Create an `AlloyGenerator` instance using `config.system` and `config.generator`.
        2.  Call `generator.generate()` to get the structures.
        3.  Create an `AseDBWrapper` instance with `db_path`.
        4.  Call `db_wrapper.write_structures()` with the generated structures.
        5.  (Optional) Add print statements or logging to indicate progress.

6.  **CLI (`cli.py`)**:
    *   Set up a Typer `app`.
    *   Define the `run` command function, decorated with `@app.command()`.
    *   Use `typer.Option` to define the `--config-path` and `--db-path` arguments, making them required. Add help text.
    *   Inside the `run` function, call `FullConfig.from_yaml` to parse the config.
    *   Instantiate `WorkflowOrchestrator` and call `execute`.
    *   Wrap the execution in a `try...except` block to catch potential errors (e.g., file not found, invalid config) and print a user-friendly error message.

## 5. Test Strategy

The testing for Cycle 01 is crucial as it validates the foundation of the entire application.

**Unit Testing Approach (Min 300 words):**

Unit tests will be created to validate each component in isolation, using mocking to remove external dependencies like the filesystem.

*   **`test_config.py`**: We will test the Pydantic models thoroughly.
    *   **Success Case**: A test will be written that uses a valid sample YAML file (stored as a string or in a test assets folder) to create a `FullConfig` object. We will then assert that the parsed fields (e.g., `config.system.elements`) match the expected values.
    *   **Failure Cases**: Multiple tests will be written to ensure validation logic works. We will test for `ValidationError` when:
        *   A required field like `elements` is missing.
        *   A field has the wrong type (e.g., `num_structures` is a string).
        *   The composition keys do not match the elements list.
        *   The composition values do not sum to 1.0.
        *   An unknown key is provided in the YAML file.
    These tests ensure our configuration contract is robust.

*   **`test_db.py`**: The `AseDBWrapper` will be tested using an in-memory database to avoid filesystem side effects.
    *   We will create a list of mock `ase.Atoms` objects.
    *   We will instantiate `AseDBWrapper` with a path to a temporary database file.
    *   After calling `write_structures`, we will use `ase.db.connect` directly in the test to read from the temporary database and assert that the number of rows matches the number of structures we wrote. We will also read one of the rows back and check that its properties (e.g., chemical symbols) match the original `ase.Atoms` object.

*   **`test_generators.py`**: The `AlloyGenerator` will be tested to verify its core logic.
    *   We will create a sample `SystemConfig` and `GeneratorConfig`.
    *   We will call the `generate` method and make several assertions on the result:
        *   The length of the returned list is equal to `num_structures`.
        *   The composition of elements in each generated `ase.Atoms` object is correct.
        *   The `_is_valid` method is correctly filtering out structures with overlapping atoms. This can be tested by temporarily setting `min_atomic_distance` to a very large value, which should result in the generator returning an empty list.

**Integration Testing Approach (Min 300 words):**

A single, comprehensive integration test will ensure that all the units work together as a cohesive whole.

*   **`test_pipeline.py`**: This test will simulate a real user running the application from the command line.
    *   **Setup**: The test will use pytest's `tmp_path` fixture to create a temporary directory for all artifacts. It will create a sample `config.yaml` file inside this directory with a simple configuration (e.g., generate 5 structures of SiC).
    *   **Execution**: The test will use `typer.testing.CliRunner` to invoke the CLI command. The command will be run like this: `runner.invoke(app, ["--config-path", str(config_file), "--db-path", str(db_file)])`. We will check that the command exits with a success code (`result.exit_code == 0`).
    *   **Validation**: After the CLI command has finished, the test will perform assertions on the side effects. It will check that the database file has been created at the specified path.
    *   **Verification**: Crucially, the test will then use `ase.db.connect` to open the output database. It will query the database to verify that it contains exactly 5 rows (the `num_structures` from our sample config). It will also read the first row and confirm that the atoms present are indeed Silicon and Carbon. This end-to-end test provides high confidence that the core workflow of Cycle 01 is functioning correctly, from the user's input (config file) to the final output (database).
