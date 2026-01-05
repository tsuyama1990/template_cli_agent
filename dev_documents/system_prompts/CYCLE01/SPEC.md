# Cycle 1 Specification: Core Pipeline and Basic Generation/Storage

## 1. Summary

Cycle 1 marks the foundational stage of the MLIP-AutoPipe project. The primary objective of this cycle is to establish a robust and operational skeleton for the entire pipeline. This involves creating a minimal, yet fully functional, end-to-end workflow that can be executed from a single command-line interface (CLI) command. The scope is intentionally limited to the most critical components to ensure a solid base is built before adding complexity.

Specifically, this cycle will deliver:
- A complete project directory structure with initial dependency management via `pyproject.toml`.
- A Typer-based CLI that serves as the main user entry point, responsible for parsing arguments and initiating the workflow.
- A core `WorkflowOrchestrator` that controls the sequence of operations.
- A dedicated database wrapper (`AseDBWrapper`) to handle all interactions with the ASE-managed SQLite database, ensuring that data persistence logic is centralized and abstracted.
- The first, and simplest, structure generator: an `AlloyGenerator` capable of creating random binary alloy structures based on user-defined composition and lattice parameters. This includes implementing a `BaseGenerator` to enforce common validation rules, such as checking for overlapping atoms.
- A `GenerationEngine` to manage the logic of calling the appropriate generator.
- A `StorageEngine` to take the generated structures and persist them using the `AseDBWrapper`.

At the conclusion of Cycle 1, the system will not perform any complex simulations or sampling; the Exploration and Sampling engines will be implemented as simple pass-through stubs. However, a user will be able to define a simple binary alloy system in a configuration file and, with one command, generate a set of physically plausible seed structures and have them correctly stored in a queryable database. This provides the essential groundwork and a tangible, testable product upon which all future cycles will build.

## 2. System Architecture

The focus of Cycle 1 is to create the core files that constitute the pipeline's backbone. The file structure will be established as follows, with the new or modified files for this cycle marked in bold.

```
src/mlip_autopipec/
├── __init__.py
├── **cli.py**                  # Typer-based command-line interface
├── config.py               # Pydantic models for Hydra configuration (stubbed)
├── database/
│   ├── __init__.py
│   └── **ase_db_wrapper.py**   # Wrapper for ASE DB interactions
├── domain/
│   ├── __init__.py
│   ├── models.py           # Core Pydantic data models (stubbed)
│   └── interfaces.py       # Abstract interfaces (stubbed)
├── engines/
│   ├── __init__.py
│   ├── **generation_engine.py**  # Orchestrates structure generation
│   ├── exploration_engine.py # (stubbed)
│   ├── sampling_engine.py    # (stubbed)
│   └── **storage_engine.py**     # Orchestrates final database storage
├── generators/
│   ├── __init__.py
│   ├── **base_generator.py**     # Base class for all generators
│   └── **alloy_generator.py**    # Implementation for alloys
├── utils/
│   └── __init__.py
└── **workflow_orchestrator.py**  # Main orchestrator class
```

**File Blueprints:**

*   **`cli.py`**:
    *   This file will contain the main entry point for the application.
    *   It will use the `Typer` library to define a simple CLI command.
    *   The command will accept a single argument: the path to the Hydra configuration file.
    *   Its primary responsibility is to instantiate the `WorkflowOrchestrator` and trigger its `run()` method.

*   **`database/ase_db_wrapper.py`**:
    *   Will contain the `AseDBWrapper` class.
    *   This class will manage the connection to the SQLite database using `ase.db`.
    *   It will provide a high-level API for the rest of the application, including methods like `connect(db_path)`, `write_atoms(atoms_list)`. This abstracts away the low-level ASE DB calls.

*   **`generators/base_generator.py`**:
    *   Will define an abstract base class `BaseGenerator` using Python's `abc` module.
    *   It will have an abstract method `generate(self) -> list[ase.Atoms]`.
    *   It will contain concrete validation methods that are common to all generators, such as `_check_overlap(atoms)`, which ensures no two atoms are too close, preventing physically impossible structures.

*   **`generators/alloy_generator.py`**:
    *   Will contain the `AlloyGenerator` class, which inherits from `BaseGenerator`.
    *   It will be initialized with a configuration object specifying elements, composition, lattice type, etc.
    *   The `generate()` method will implement the logic to create a supercell, populate it with the specified elements in the correct ratio at random positions, and then run the validation checks from the base class.

*   **`engines/generation_engine.py`**:
    *   Will contain the `GenerationEngine` class.
    *   Its role is to select and run the correct generator based on the input configuration. For this cycle, it will be hardcoded to use the `AlloyGenerator`.
    *   It will return a list of the generated `ase.Atoms` objects.

*   **`engines/storage_engine.py`**:
    *   Will contain the `StorageEngine` class.
    *   It will take a list of `ase.Atoms` objects and an instance of the `AseDBWrapper`.
    *   Its `run()` method will simply call the `write_atoms()` method on the database wrapper.

*   **`workflow_orchestrator.py`**:
    *   Will contain the `WorkflowOrchestrator` class.
    *   Its constructor will initialize all the required engine components.
    *   Its `run()` method will execute the pipeline stages in the correct order:
        1.  Call `GenerationEngine` to create the structures.
        2.  (Call stubbed Exploration and Sampling engines).
        3.  Call `StorageEngine` to save the final structures.

## 3. Design Architecture

This cycle establishes the Pydantic-based design for configuration management, ensuring that all inputs are validated and strongly typed from the very beginning. This "schema-first" approach is critical for building a robust system.

**Pydantic Models (`config.py`):**

The initial configuration schema will be simple, focusing only on the parameters needed for the `AlloyGenerator`.

```python
from pydantic import BaseModel, Field
from typing import list, dict

class AlloyGeneratorConfig(BaseModel):
    """Configuration for the AlloyGenerator."""
    elements: list[str] = Field(..., min_length=2, description="List of element symbols.")
    composition: dict[str, float] = Field(..., description="Dictionary mapping element to its fraction.")
    lattice_constant: float = Field(..., gt=0, description="Lattice constant in Angstroms.")
    crystal_structure: str = Field("fcc", description="Crystal structure (e.g., 'fcc', 'bcc').")
    supercell_size: list[int] = Field([5, 5, 5], description="Size of the supercell.")
    num_structures: int = Field(10, gt=0, description="Number of structures to generate.")

class MainConfig(BaseModel):
    """Main configuration model."""
    generator: AlloyGeneratorConfig
    db_path: str = "structures.db"

```

**Design Principles:**

*   **Data Integrity:** By defining `AlloyGeneratorConfig` with Pydantic, we ensure that any configuration used to run the pipeline is automatically validated. For example, `lattice_constant` must be a positive number (`gt=0`), and `elements` must be a list containing at least two strings. This prevents a large class of runtime errors caused by invalid configuration.
*   **Producer/Consumer:** The primary producer of this configuration is the user, who writes a YAML file. The Hydra library, in conjunction with our Pydantic models, acts as the consumer and validator, loading the YAML and parsing it into a `MainConfig` object. This typed object is then passed down through the `WorkflowOrchestrator` to the relevant components, ensuring type safety throughout the application.
*   **Extensibility:** This structure is highly extensible. In future cycles, we can add new configuration models (e.g., `IonicGeneratorConfig`) and use a discriminated union in `MainConfig` to allow the user to select the generator type directly from the configuration file. For now, we establish the pattern of a single, top-level configuration object that contains nested, component-specific configurations.

## 4. Implementation Approach

The implementation will proceed in a logical, step-by-step manner to ensure each component is built and tested correctly.

1.  **Project Setup:**
    *   Create the directory structure as outlined in the System Architecture.
    *   Initialize a `pyproject.toml` file.
    *   Add initial dependencies: `typer`, `pydantic`, `ase`, `hydra-core`, `ruff`, `pytest`. Run `uv sync` to create the virtual environment.

2.  **Database Wrapper (`AseDBWrapper`):**
    *   Create the `ase_db_wrapper.py` file.
    *   Implement the `AseDBWrapper` class.
    *   The `connect` method will take a file path and use `ase.db.connect(path)` to establish a connection.
    *   The `write_atoms` method will iterate through a list of `ase.Atoms` objects and use the `db.write(atoms)` method for each one.

3.  **Generators (`BaseGenerator`, `AlloyGenerator`):**
    *   Create `base_generator.py` and define the `BaseGenerator` ABC.
    *   Implement the `_check_overlap` method, which calculates the distance matrix of an atoms object and raises a `ValueError` if any distance is smaller than a reasonable threshold (e.g., 1.0 Angstrom).
    *   Create `alloy_generator.py`. The `AlloyGenerator` will inherit from `BaseGenerator`.
    *   The constructor will accept an `AlloyGeneratorConfig` object.
    *   The `generate` method will use `ase.build.bulk` to create a primitive cell, `ase.build.make_supercell` to expand it, and then randomly assign elements to atomic positions based on the specified composition. It will call `_check_overlap` before returning the final structure.

4.  **Engines (`GenerationEngine`, `StorageEngine`):**
    *   Implement `GenerationEngine`. It will instantiate `AlloyGenerator` with the configuration and call its `generate` method.
    *   Implement `StorageEngine`. It will take the list of atoms from the previous step and the database wrapper instance and call the wrapper's `write_atoms` method.

5.  **Orchestrator and CLI (`WorkflowOrchestrator`, `cli.py`):**
    *   Implement `WorkflowOrchestrator`. Its `run` method will instantiate and call the engines in sequence. It will handle the instantiation of the `AseDBWrapper` and pass it to the `StorageEngine`.
    *   Implement `cli.py`. Use Typer to create a main function that accepts the config path.
    *   Inside the main function, use Hydra to load the YAML configuration and parse it into the `MainConfig` Pydantic model.
    *   Instantiate `WorkflowOrchestrator` with the loaded config and call `run()`.

## 5. Test Strategy

Testing in Cycle 1 is crucial to validate the core architecture and ensure a stable foundation.

**Unit Testing Approach (Min 300 words):**

Unit tests will focus on isolating each new component to verify its correctness independently.
*   **`AseDBWrapper`:** We will use `pytest` and Python's `unittest.mock` library. The tests will not write to the actual filesystem. We can use an in-memory SQLite database (`db_path=":memory:"`) for testing. We will create a test that calls `write_atoms` with a list of mock `ase.Atoms` objects and then uses the underlying `ase.db` API to assert that the correct number of rows were written to the database.
*   **`AlloyGenerator`:** The generator will be tested by instantiating it with a fixed `AlloyGeneratorConfig` and calling `generate()`. The assertions will be comprehensive:
    1.  Check that the returned object is a list of `ase.Atoms` objects.
    2.  Verify that the number of returned structures matches `num_structures` in the config.
    3.  For each structure, verify the total number of atoms is correct based on the supercell size.
    4.  Check that the chemical symbols in each structure match the `elements` in the config.
    5.  Verify that the composition of the generated structure is approximately correct (within statistical fluctuations for a random assignment).
    6.  Crucially, test the validation logic by creating a configuration that is guaranteed to produce overlapping atoms (e.g., a very small lattice constant) and assert that the `generate` method raises the expected `ValueError`.
*   **CLI:** The Typer CLI will be tested using `typer.testing.CliRunner`. We will test that the command runs successfully with a valid config path and that it fails with an appropriate error message if the config path is invalid or does not exist.

**Integration Testing Approach (Min 300 words):**

Integration testing will verify that the components work together as a complete system. The primary integration test will simulate the user's end-to-end workflow.
*   **Test Setup:** The test will use a `pytest` fixture to create a temporary directory using `tmp_path`. A simple, valid Hydra YAML configuration file for the `AlloyGenerator` will be created programmatically within this temporary directory.
*   **Execution:** The test will invoke the CLI command using `CliRunner`, passing the path to the newly created configuration file. The database path in the config will also point to a file within the temporary directory.
*   **Assertions:** After the CLI command finishes, the test will perform a series of assertions to verify the entire process was successful:
    1.  Check the command's exit code is 0, indicating success.
    2.  Assert that the specified SQLite database file now exists in the temporary directory.
    3.  Connect to the generated database using `ase.db.connect()`.
    4.  Query the database to get the total number of rows (structures) and assert that it matches the `num_structures` specified in the input configuration file.
    5.  Retrieve one or two structures from the database and perform basic checks on them (e.g., verify atom counts and symbols) to ensure the data was written correctly. This test validates the entire flow: CLI parsing, Hydra configuration loading, Pydantic validation, orchestrator execution, generator logic, and database writing.
