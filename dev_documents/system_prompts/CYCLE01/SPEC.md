# CYCLE 01: SPECIFICATION - Core Command-Line Pipeline

## 1. Summary

This document provides the detailed technical specification for Cycle 1 of the MLIP-AutoPipe project. The primary goal of this cycle is to build the foundational backbone of the application: a robust, command-line-driven pipeline capable of executing an end-to-end workflow for generating training data for a specific, well-defined use case—binary and ternary alloys. This cycle will focus on implementing the core architectural patterns, including the Pydantic-based configuration, the decoupled orchestrator-engine-database design, and the essential components for structure generation, simplified exploration, and database storage.

The scope of Cycle 1 is intentionally limited to establish a solid, testable foundation before introducing more complex scientific features. The exploration stage will be implemented using a basic Molecular Dynamics (MD) simulation with a standard, non-machine-learning potential (like Effective Medium Theory, EMT) to validate the parallel processing and data handling logic. Advanced features such as the hybrid MD/MC engine, Farthest Point Sampling, and support for diverse material systems (ionic, covalent, etc.) are explicitly deferred to Cycle 2. The key deliverables for this cycle are a functional Command-Line Interface (CLI), a well-defined database schema wrapper, the `AlloyGenerator`, a simplified `LabelingEngine` for exploration, and the `WorkflowOrchestrator` that ties them all together. At the end of this cycle, a user will be able to write a simple YAML configuration file, run a single command, and produce a valid ASE database of atomic structures for an alloy system.

## 2. System Architecture

The architecture for Cycle 1 implements the core components of the system as defined in the main `SYSTEM_ARCHITECTURE.md`. The focus is on establishing the foundational file structure and class interactions.

**File Structure (Cycle 1 Focus):**
The files and directories to be created or modified in this cycle are marked in **bold**.

```
src/mlip_autopipec/
├── cli/
│   └── **main.py**              # CLI entry point (Click-based)
├── pipeline/
│   └── **orchestrator.py**      # Contains the main WorkflowOrchestrator class
│   └── **interfaces.py**        # Defines abstract interfaces (e.g., ILabelingEngine)
│   └── **factories.py**         # Defines factories for creating concrete classes
├── generators/
│   ├── __init__.py
│   ├── **base.py**              # Abstract BaseStructureGenerator class
│   ├── **factory.py**           # GeneratorFactory
│   └── **alloy.py**             # AlloyGenerator implementation
├── exploration/
│   ├── __init__.py
│   └── **engine.py**            # Core MD engine (Simplified for Cycle 1)
├── storage/
│   ├── __init__.py
│   └── **database.py**          # ASE Database wrapper
├── config/
│   ├── __init__.py
│   └── **schema.py**            # Pydantic models for configuration validation
└── common/
    └── **atoms_utils.py**       # Utility functions for manipulating ASE Atoms objects
tests/
├── unit/
│   ├── **test_config_schema.py**
│   ├── **test_database.py**
│   ├── **test_alloy_generator.py**
│   └── **test_orchestrator.py**
└── integration/
    └── **test_cli_pipeline.py**
```

**Architectural Blueprint:**

The system's design is rooted in the Dependency Inversion Principle. The high-level `WorkflowOrchestrator` does not depend on concrete implementations of engines or generators. Instead, it depends on abstractions (`interfaces.py`). This allows for a loosely coupled system where components can be easily swapped or replaced.

1.  **Configuration (`config/schema.py`):** A comprehensive set of Pydantic models will be created to validate the user's configuration file. This includes a top-level `FullConfig` model that nests `SystemConfig`, `ExplorationConfig`, etc. For Cycle 1, the `SystemConfig` will include fields specific to alloys, like `alloy_type` (e.g., 'random') and `lattice_constant`. This strict, schema-first approach is critical for preventing configuration errors at runtime.

2.  **Database Wrapper (`storage/database.py`):** The `AseDBWrapper` class will be implemented to provide a clean, high-level API for all database operations. It will abstract away the low-level details of the `ase.db` library. Key methods will include:
    *   `__init__(self, db_path: str)`: Connects to or creates the database.
    *   `add_initial_structures(self, atoms_list: list[Atoms])`: Adds the seed structures, marking them as 'unlabeled'.
    *   `get_unlabeled_structure_ids() -> list[int]`: Retrieves the IDs of structures that need to be processed by the exploration engine.
    *   `get_atoms_by_id(self, id: int) -> Atoms`: Fetches a single structure.
    *   `update_with_labeled_result(self, original_id: int, result: DFTResult)`: Updates the database with the results from the exploration, marking the entry as 'labeled'. The `DFTResult` will be a Pydantic model containing energy, forces, and stress.

3.  **Generator (`generators/`):**
    *   `base.py`: The `BaseStructureGenerator` abstract base class will define the public interface: `generate(config: SystemConfig) -> list[Atoms]`. It will also contain concrete implementations of essential validation logic like checking for atomic overlap, which is common to all generators.
    *   `alloy.py`: The `AlloyGenerator` will implement the `generate` method. It will use ASE's `ase.build.bulk` to create a primitive cell, construct a supercell, and then randomly replace atomic species to achieve the target composition. It will rigorously apply the validation methods from its parent class.
    *   `factory.py`: The `GeneratorFactory` will be a simple function that takes the system configuration and returns an `AlloyGenerator` instance.

4.  **Labeling Engine (`exploration/engine.py`):**
    *   The `ILabelingEngine` interface in `pipeline/interfaces.py` will define a single method: `run(atoms: Atoms) -> DFTResult`.
    *   The concrete `LabelingEngine` will implement this interface. Its `run` method will attach a simple ASE calculator (e.g., `ase.calculators.emt.EMT`) to the `Atoms` object and run a short MD simulation (`ase.md.langevin.Langevin`). After the run, it will extract the final energy, forces, and stress, package them into a `DFTResult` Pydantic model, and return it. This isolates the "doing" of the simulation from the "managing" of the workflow.

5.  **Orchestrator (`pipeline/orchestrator.py`):**
    *   The `WorkflowOrchestrator` is the central coordinator. Its constructor will receive instances of the database wrapper, the generator, and the labeling engine (dependency injection).
    *   Its main `run()` method will implement the core logic:
        1.  Check if initial structures exist in the DB. If not, call the generator and save the results to the DB.
        2.  Get the list of all unlabeled structures from the DB.
        3.  In a loop, for each structure:
            a. Fetch the `Atoms` object from the DB.
            b. Pass the `Atoms` object to the `LabelingEngine.run()` method.
            c. Save the returned `DFTResult` back to the DB using the `update` method.

6.  **CLI (`cli/main.py`):**
    *   A `click`-based CLI will be created.
    *   It will have a single command, `run-pipeline`, which takes one argument: the path to a YAML configuration file.
    *   The command will first parse the YAML file and validate it against the `FullConfig` Pydantic model.
    *   It will then use the factory functions from `pipeline/factories.py` to instantiate the `AseDBWrapper`, `AlloyGenerator`, and `LabelingEngine`.
    *   Finally, it will create an instance of `WorkflowOrchestrator` with these dependencies and call its `run()` method.

## 3. Design Architecture

This cycle establishes the project's commitment to a schema-driven design, where Pydantic models are the single source of truth for data structures passing between modules. This ensures data integrity and provides clear, machine-readable contracts.

**Pydantic-based Schema Design:**

1.  **`config.schema.FullConfig`:** This is the top-level model, representing the entire user-provided configuration. It enforces the overall structure of the input file.
    *   **Consumers:** The CLI (`main.py`) is the primary consumer. It deserializes and validates the user's YAML file into this model at the very start of the program.
    *   **Producers:** The user is the producer, writing a YAML file.
    *   **Constraints:** `extra='forbid'` will be used to prevent users from providing unknown configuration keys, catching typos early.

2.  **`storage.database.DFTResult`:** This model represents the structured output of a single simulation run.
    *   **Producers:** The `LabelingEngine` is the sole producer. After a simulation, it must package its findings (energy, forces, stress) into this strict format.
    *   **Consumers:** The `WorkflowOrchestrator` consumes this model and passes it to the `AseDBWrapper` for persistence.
    *   **Constraints:**
        *   `energy: float`: A standard floating-point number.
        *   `forces: list[list[float]]`: A nested list representing a NumPy array of shape (N, 3). A custom validator will ensure the list of lists has the correct shape.
        *   `stress: list[float]`: A list of 6 floats representing the Voigt-notation stress tensor. A validator will ensure the list has exactly 6 elements.
    *   **Extensibility:** This model is designed for future extension. For example, a `trajectory: list[Atoms]` field could be added later without breaking existing components.

3.  **Dependency Injection and Interfaces (`pipeline/interfaces.py`, `pipeline/factories.py`):**
    *   The design explicitly decouples the `WorkflowOrchestrator` from concrete engines and generators.
    *   `interfaces.py` will define the contracts, e.g., `class ILabelingEngine(Protocol): def run(self, atoms: Atoms) -> DFTResult: ...`.
    *   The `WorkflowOrchestrator`'s constructor will be type-hinted against these interfaces: `__init__(self, ..., engine: ILabelingEngine)`.
    *   `factories.py` will contain functions like `create_engine(config: ExplorationConfig) -> ILabelingEngine`. This function will read the config and decide which concrete engine to instantiate (for Cycle 1, always the `LabelingEngine`).
    *   This design is crucial for testability, as it allows unit tests for the orchestrator to inject mock engines. It is also key for future extensibility, as a new engine (e.g., a Quantum Espresso engine) could be added just by updating the factory, with no changes to the orchestrator itself.

## 4. Implementation Approach

The implementation will proceed in a logical order, starting from the foundational data structures and moving up to the high-level orchestration.

1.  **Project Scaffolding:** Create the directory structure as outlined in the "System Architecture" section using `mkdir -p`. Create empty `__init__.py` files to make them Python packages.
2.  **Pydantic Schemas (`config/schema.py`):** Implement `FullConfig` and all nested configuration models first. This defines the "shape" of the data for the entire application.
3.  **Database Wrapper (`storage/database.py`):** Implement the `AseDBWrapper` class. This component is self-contained and can be developed and unit-tested in isolation. It will depend only on the `ase` library.
4.  **Generator Implementation (`generators/`):**
    a. Create the `BaseStructureGenerator` ABC in `base.py` with its validation methods.
    b. Create the `AlloyGenerator` in `alloy.py`, inheriting from the base class.
    c. Implement the `GeneratorFactory` in `factory.py`.
    d. Write unit tests for the generator, ensuring it produces valid structures and rejects invalid ones.
5.  **Labeling Engine Implementation (`exploration/engine.py`):**
    a. Define the `ILabelingEngine` protocol in `pipeline/interfaces.py`.
    b. Implement the `DFTResult` Pydantic model.
    c. Create the `LabelingEngine` class. For now, its `run` method will use the simple `ase.calculators.emt.EMT`.
    d. Unit test the engine in isolation. A mock `Atoms` object can be used to test that the engine correctly calls the calculator and returns a properly formatted `DFTResult`.
6.  **Orchestrator Implementation (`pipeline/orchestrator.py`):**
    a. Implement the `WorkflowOrchestrator` class. Its constructor will accept the database wrapper, generator, and engine.
    b. Implement the `run` method following the logic defined in the blueprint above.
    c. Write unit tests for the orchestrator using `unittest.mock`. Mock objects for the database, generator, and engine will be injected to verify that the orchestrator calls them in the correct sequence and with the correct data.
7.  **CLI Implementation (`cli/main.py`):**
    a. Use the `click` library to create the `run-pipeline` command.
    b. The command logic will handle file I/O (reading the YAML), configuration validation with Pydantic, and the setup of the dependency injection container (i.e., calling the factories).
8.  **Integration Testing (`tests/integration/test_cli_pipeline.py`):**
    a. Create a final integration test that simulates a user running the command.
    b. Use `click.testing.CliRunner` to invoke the CLI command within the test.
    c. The test will create a temporary directory, write a minimal `config.yml` file, and run the pipeline.
    d. Assertions will be made against the final state of the temporary database file to ensure the pipeline ran successfully and produced the expected output.

## 5. Test Strategy

**Unit Testing Approach (Min 300 words):**
The unit testing strategy for Cycle 1 is focused on component isolation. Each class will be tested independently of its collaborators to ensure its internal logic is correct. This is made possible by the dependency injection architecture. We will use `pytest` as the test runner and `unittest.mock` for creating mock objects.

*   **Configuration (`test_config_schema.py`):** We will test the Pydantic models directly. Tests will ensure that valid data successfully creates a model instance, while invalid data (e.g., `temperature_K = -10`) raises a `pydantic.ValidationError`. We will also test that providing an unknown key in the configuration raises an error, enforcing our strict schema.
*   **Database (`test_database.py`):** The `AseDBWrapper` will be tested by patching the `ase.db.connect` function at the module level. This means our tests will not touch the actual filesystem. We will create a mock connection object and assert that our wrapper's methods (e.g., `add_initial_structures`) result in the correct calls to the mock's methods (e.g., `mock_connection.write()`). We will verify that the arguments passed are correctly formatted.
*   **Generator (`test_alloy_generator.py`):** The `AlloyGenerator` will be tested by calling its `generate` method with a sample `SystemConfig`. We will assert that the returned list of `ase.Atoms` objects has the correct length, that each object has the right number of atoms, and that the chemical composition matches the request. We will also write tests that are *expected* to fail, for example by requesting a configuration that would result in overlapping atoms, and assert that a custom `PhysicsViolationError` is raised.
*   **Orchestrator (`test_orchestrator.py`):** This is a key unit test. We will create mock objects for `AseDBWrapper`, `AlloyGenerator`, and `LabelingEngine`. We will then instantiate the `WorkflowOrchestrator` with these mocks. The test will call the `run()` method and then make assertions on the mocks, for example: `mock_generator.generate.assert_called_once()`, `mock_engine.run.assert_called_with(test_atoms)`, and `mock_db.update_with_labeled_result.assert_called()`. This verifies the orchestration logic without needing to perform any real computation or I/O.

**Integration Testing Approach (Min 300 words):**
The integration testing strategy for Cycle 1 is designed to verify that the isolated components, once assembled, work together correctly. The primary integration test will focus on the highest level of integration: the command-line interface.

*   **CLI End-to-End Test (`test_cli_pipeline.py`):** We will use the `click.testing.CliRunner` and `pytest`'s `tmp_path` fixture. The test case will perform the following steps:
    1.  **Setup:** Create a temporary directory (`tmp_path`). Inside this directory, create a simple YAML configuration file (`config.yml`) for a small binary alloy (e.g., CuAu). The configuration will specify a small number of initial structures and a very short MD run.
    2.  **Execution:** Use `CliRunner.invoke()` to run the main CLI command, passing the path to the `config.yml` file. We will use `catch_exceptions=False` to get a full traceback if the application fails.
    3.  **Verification:** After the command finishes, the test will inspect the state of the temporary directory.
        *   It will check that an SQLite database file (e.g., `results.db`) has been created.
        *   It will connect to this database using the real `ase.db.connect`.
        *   It will query the database to assert that the correct number of initial structures were created and that all of them have been successfully 'labeled' (i.e., processed by the engine).
        *   It will check the row count of the database and verify it matches the number of generated structures.
This test covers the entire application stack, from CLI parsing and configuration validation to the final database write, providing a high degree of confidence that all the components are correctly integrated. It acts as a safety net to catch issues that unit tests, with their heavy use of mocking, might miss, such as mismatched data formats between components.
