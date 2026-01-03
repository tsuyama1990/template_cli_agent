# SPEC.md for CYCLE01: Core Data Pipeline & Manual Workflow

## 1. Summary

This document provides the detailed technical specification for CYCLE01 of the MLIP-AutoPipe project. The primary, strategic objective of this cycle is to construct the foundational infrastructure and the core data processing pipeline of the application. This initial phase is deliberately constrained and does not concern itself with the complex automation or active learning features that will ultimately define the project. Instead, its focus is squarely on creating a robust, reliable, and manually-operable workflow. This workflow will enable a user to take a single, well-defined atomic structure, command the system to perform a high-fidelity DFT calculation to "label" it with its corresponding energy, forces, and stresses, store these results persistently and safely in a structured database, and finally, train a basic Machine Learning Interatomic Potential (MLIP) from the accumulated data. This cycle is analogous to building the factory floor, installing the heavy machinery, and verifying that each machine works correctly in isolation before attempting to automate the entire assembly line.

At the successful conclusion of CYCLE01, the system will possess several key, well-tested components that form the backbone of the entire application. These include: strictly-defined Pantic models for all configuration management, ensuring type safety and validation; a dedicated database wrapper providing a clean, high-level API for all data persistence, abstracting away the underlying database implementation; a reliable `LabelingEngine` capable of robustly interfacing with Quantum Espresso to perform single-point calculations; and a functional `TrainingEngine` for generating valid ACE potentials from the labeled data. The interaction between these disparate components will be managed by a simple, sequential `WorkflowOrchestrator`, which in turn will be exposed to the end-user via a clear and simple command-line interface. This cycle is arguably the most critical, as it lays the essential groundwork and validates the core architectural decision of a modular, database-mediated workflow. The success of this cycle will be measured by the system's ability to reliably and reproducibly produce a valid, loadable MLIP from a manually curated set of input structures, thereby providing a solid, trustworthy foundation upon which the sophisticated automation logic of CYCLE02 can be built.

## 2. System Architecture

The architecture for CYCLE01 is a carefully selected subset of the full system architecture, focusing exclusively on the components required to stand up the manual data pipeline. The system will be developed as a modern Python package named `mlip_autopipec`, adhering to current best practices for packaging and dependency management as defined by `pyproject.toml`.

**File Structure for CYCLE01:**
(Bold files and directories are those that will be created or significantly modified during this cycle)

```
.
├── **pyproject.toml**
├── README.md
├── **src**
│   └── **mlip_autopipec**
│       ├── **__init__.py**
│       ├── **cli.py**                # User-facing CLI application logic
│       ├── **config.py**             # Pydantic data schemas for configuration
│       ├── **orchestrator.py**       # Core workflow logic
│       ├── **database.py**           # Data persistence layer
│       └── **modules**
│           ├── **__init__.py**
│           ├── **c_labeling_engine.py** # Quantum Espresso interface
│           └── **d_training_engine.py** # ACE model training interface
└── **tests**
    ├── **__init__.py**
    ├── **test_config.py**
    ├── **test_orchestrator.py**
    └── **modules**
        ├── **__init__.py**
        ├── **test_c_labeling_engine.py**
        └── **test_d_training_engine.py**
```

**Component Blueprint:**

*   **`pyproject.toml`**: This file is the cornerstone of the project's packaging and dependency story. It will be initialized using `uv` and will clearly define the project's metadata, such as its name (`mlip-autopipec`), version, and authors. Most importantly, it will list all direct dependencies required for this cycle, including `ase` for atomic structure manipulation, `pydantic` for configuration modeling, `typer` for the command-line interface, and `pyyaml` for parsing configuration files. It will also define a separate group for development dependencies, such as `pytest` and `pytest-mock`, and will be configured to run tests from the `tests` directory.

*   **`src/mlip_autopipec`**: This directory will contain all the source code for the installable Python package.

*   **`config.py`**: This file is of critical importance as it embodies the "Schema-First" design principle. It will contain the Pydantic models that define the canonical structure of all configuration parameters. For CYCLE01, it will primarily define the `FullConfig` model, which will be loaded directly from a user-created YAML file. This model will be composed of nested sub-models for DFT settings (`DFTSettings`) and MLIP training settings (`MLIPSettings`), ensuring that every parameter is validated at runtime before any computation begins. This prevents a large class of errors related to misconfiguration.

*   **`database.py`**: This module will implement the `AseDBWrapper` class. The purpose of this class is to provide a high-level abstraction layer over the `ase.db` SQLite database, decoupling the rest of the application from the specific implementation details of the database. It must provide clear, intention-revealing methods like `connect(db_path)`, `add_atoms(atoms_obj, metadata)`, `get_all_labeled_atoms()`, and `get_atom_by_id(id)`. It will be responsible for managing the database connection lifecycle and ensuring that all data is written and read in a consistent, predictable format.

*   **`modules/c_labeling_engine.py`**: This file will contain the `LabelingEngine` class. Its sole responsibility is to serve as a robust and reliable interface to the Quantum Espresso DFT code. It will be initialized with the `FullConfig` object. Its main method will accept an `ase.Atoms` object, generate a valid Quantum Espresso input file in a temporary directory, execute the `pw.x` binary as a subprocess, and then meticulously parse the text-based output file to extract the key results: total energy, atomic forces, and the virial stress tensor. In this cycle, it will implement basic error handling, primarily by checking the exit code of the subprocess and raising a specific, informative exception if the DFT calculation fails to complete successfully.

*   **`modules/d_training_engine.py`**: This file will contain the `TrainingEngine` class. This class is responsible for the machine learning aspect of the pipeline. It will query the database, via the `AseDBWrapper`, to fetch the complete set of labeled structures. It will then implement the "Delta Learning" strategy by calculating the energies and forces from a simple baseline potential (e.g., Lennard-Jones, with parameters derived from heuristics) and subtracting them from the raw DFT values. It will then use these computed residuals to train a new ACE potential by interfacing with a suitable external library (e.g., `pace-neutrons` or a similar package). The final trained potential file (e.g., `model.ace`) will be saved to a user-specified output directory.

*   **`orchestrator.py`**: This file will house the `WorkflowOrchestrator` class. For CYCLE01, its logic will be simple and sequential. It will be initialized with a `FullConfig` object and an instance of the `AseDBWrapper`. It will have two primary public methods: `run_labeling(atoms_file)`, which will read an `ase.Atoms` object from a file, coordinate with the `LabelingEngine` to perform the DFT calculation, and then use the `AseDBWrapper` to store the results. The second method, `run_training()`, will coordinate with the `TrainingEngine` to query the database and produce the final potential file.

*   **`cli.py`**: This module provides the primary user entry point for the application. It will be built using the `Typer` library to ensure a clean, modern, and well-documented command-line experience. It will define two distinct commands for this cycle, reflecting the manual nature of the workflow:
    *   `mlip-pipe label --config-file <path> --atoms-file <path> --db-file <path>`: This command will be responsible for creating the necessary objects (config, database wrapper), instantiating the orchestrator, and invoking the `run_labeling` method.
    *   `mlip-pipe train --config-file <path> --db-file <path>`: This command will similarly instantiate the orchestrator and invoke the `run_training` method.

## 3. Design Architecture

This cycle places a heavy emphasis on a Schema-First design, with Pydantic models serving as the unshakeable foundation of the entire architecture. All data that flows through the system, from the initial user configuration to the final DFT results, will be strongly typed, rigorously validated, and self-documenting. This approach is instrumental in building a robust and maintainable system.

**Pydantic Schema Design (`config.py`):**

The entire design is centered around a single, comprehensive, and deeply nested `FullConfig` model. This strategic choice avoids the pitfalls of implicit state and scattered configuration parameters, making the entire workflow perfectly reproducible from a single, explicit configuration file. Every parameter that can affect the outcome of a calculation will be defined here.

```python
#
# Pseudocode for config.py illustrating detailed structure
#
from pydantic import BaseModel, Field, PositiveFloat, FilePath
from typing import Dict, List, Literal

class SystemSettings(BaseModel):
    """Defines the chemical system to be modeled."""
    elements: List[str] = Field(..., min_items=1)
    composition: Optional[Dict[str, float]] = None

class DFTSettings(BaseModel):
    """Detailed settings for the Quantum Espresso DFT Calculation."""
    code: Literal["quantum_espresso"] = "quantum_espresso"
    command: str = "pw.x"
    pseudopotentials_protocol: str = "SSSP_1.3_PBE_precision"
    ecutwfc: PositiveFloat = Field(..., description="Wavefunction cutoff in Ry")
    ecutrho: PositiveFloat = Field(..., description="Density cutoff in Ry")
    kpoints_density: PositiveFloat = Field(..., description="K-point mesh density in 1/Angstrom, used for automatic mesh generation.")
    smearing: Literal["gaussian", "mv", "fd"] = "mv"
    degauss: PositiveFloat = 0.02
    magnetism: Literal["none", "ferromagnetic", "antiferromagnetic"] = "none"
    # Additional fields for convergence thresholds, mixing params etc. can be added.

class BaselinePotentialSettings(BaseModel):
    """Settings for the baseline potential in Delta Learning."""
    type: Literal["none", "lj_auto"] = "lj_auto"

class MLIPSettings(BaseModel):
    """Settings for the ACE Model Training."""
    model_type: Literal["ace"] = "ace"
    r_cut: PositiveFloat = Field(..., description="Cutoff radius in Angstrom")
    delta_learning: bool = True
    baseline: BaselinePotentialSettings = BaselinePotentialSettings()
    loss_weights: Dict[Literal["energy", "force", "stress"], PositiveFloat] = {
        "energy": 1.0,
        "force": 100.0,
        "stress": 10.0
    }

class FullConfig(BaseModel):
    """The complete, validated configuration for the MLIP-AutoPipe workflow."""
    system: SystemSettings
    dft_compute: DFTSettings
    mlip_training: MLIPSettings
```

**Key Architectural Invariants and Constraints:**
*   **Configuration Immutability**: The `FullConfig` object, once loaded and validated at the start of a command, must be treated as strictly immutable by all downstream components. Modules are allowed to read from it to guide their behavior, but they are forbidden from modifying it. This ensures that the behavior of the system is consistent and reproducible throughout a single run.
*   **Stateless Modules**: The `LabelingEngine` and `TrainingEngine` must be designed to be completely stateless. All the information they require to operate must be provided to them either during their initialization (via the `FullConfig` object) or as arguments to their primary execution methods. They must not hold any state between calls. All stateful information, such as the list of processed structures, is managed exclusively by the database and the orchestrator.
*   **Database as the Single Source of Truth**: The ASE database file is designated as the single, canonical source of truth for all scientific data, including all atomic structures and their computed properties. The local filesystem should only be used for storing configuration files, final output artifacts (like the `.ace` potential), and temporary files required by external programs (like QE input/output). This centralization of data is key to ensuring data integrity and consistency.

**Data Consumers and Producers:**
*   **User**: The primary producer of the initial `exec_config.yaml` file and one or more atomic structure files (e.g., in `.cif` or `.xyz` format). The user is the final consumer of the trained `.ace` potential file.
*   **`cli.py`**: Acts as the initial producer of the in-memory `FullConfig` object by parsing and validating the user's YAML file.
*   **`LabelingEngine`**: Consumes an `ase.Atoms` object and the `FullConfig` object. It produces the core scientific data: the calculated DFT energy, forces, and stress tensor.
*   **`AseDBWrapper`**: The primary intermediary for data. It is both a producer (when reading from the database) and a consumer (when writing to it) of `ase.Atoms` objects and their associated scientific data.
*   **`TrainingEngine`**: Consumes the `FullConfig` object and the complete set of labeled data retrieved from the database. It produces the final, tangible artifact of the workflow: the trained `.ace` model file.

## 4. Implementation Approach

The implementation for CYCLE01 will proceed in a logical, bottom-up order, starting from the foundational data structures and moving progressively outwards to the application and command-line logic. This ensures that each layer is built upon a stable, well-tested foundation.

1.  **Setup `pyproject.toml`**: The very first step is to initialize the project using `uv init`. This command creates the `pyproject.toml` file. We will then edit this file to add the core dependencies for this cycle: `pydantic`, `ase`, `typer[all]`, and `pyyaml`. We will also add a dedicated group for development dependencies: `pytest` and `pytest-mock`. Finally, we will configure the project structure, telling the build system that our package resides in the `src/mlip_autopipec` directory.
2.  **Implement `config.py`**: With the project structure in place, the next step is to create the Pydantic models as designed in the section above. This is a critical early step as these models define the "data contracts" that all other components in the system will adhere to. We will write comprehensive unit tests in `tests/test_config.py` to ensure that the validation logic is working as expected. For example, a test will confirm that attempting to create a `DFTSettings` object with a negative `ecutwfc` value correctly raises a `pydantic.ValidationError`.
3.  **Implement `database.py`**: Next, we will create the `AseDBWrapper` class. It should have a simple, clean API that hides the underlying `ase.db` implementation details. The implementation will primarily use the standard `ase.db.connect()` function. We will write dedicated tests in `tests/modules/test_database.py` that create a temporary, in-memory SQLite database, add `ase.Atoms` objects to it, retrieve them, and assert that the data written to and read from the database is perfectly consistent.
4.  **Implement `modules/c_labeling_engine.py`**: This is one of the most complex steps in this cycle. We will create the `LabelingEngine` class.
    *   The input file generation logic will use Python's f-strings or a simple templating mechanism to create the Quantum Espresso input text from the `ase.Atoms` object and the `FullConfig` settings.
    *   The `subprocess.run` call will be used to execute the `pw.x` binary. It is critical to configure this call correctly to capture `stdout`, `stderr`, and the process's return code to allow for robust error checking.
    *   The output parser will be a separate, private method within the class. It will need to be robust, using precise regular expressions to find and parse the specific lines in the output file that contain the required data, such as `!    total energy`, `Forces acting on atoms`, and `total   stress`.
    *   **Testing (`test_c_labeling_engine.py`)**: Testing this component is absolutely crucial. We will use `pytest-mock` to patch the `subprocess.run` function, completely isolating our tests from any external Quantum Espresso dependency. The tests will verify several key behaviors:
        *   Given a specific `ase.Atoms` object and config, is the generated input file text correct?
        *   Does the engine correctly parse a sample, pre-saved, valid Quantum Espresso output file?
        *   Does the engine raise an appropriate, specific exception if `subprocess.run` indicates a failure (e.g., returns a non-zero exit code) or if the output file is malformed and cannot be parsed?
5.  **Implement `modules/d_training_engine.py`**: We will then create the `TrainingEngine` class. This class will likely act as a wrapper around an external library for the actual ACE model training. The engine's primary job is to prepare the data in the precise format required by that library. The delta learning logic (subtracting the baseline potential) will be implemented here. Tests will involve creating a mock database with a few atoms and ensuring that the external training function is called with the correctly transformed data.
6.  **Implement `orchestrator.py`**: The `WorkflowOrchestrator` will be implemented to tie all the previously built components together. It will instantiate the engines and the database wrapper. Its methods will contain the high-level application logic, for example: "read atoms from this file, pass them to the labeling engine, and then store the results in the database."
7.  **Implement `cli.py`**: The final implementation step is to create the Typer CLI. This will involve parsing file paths and other arguments from the command line, using `pyyaml` to load the `FullConfig` from the specified file, and then initializing and running the orchestrator. We will write integration tests for the CLI using Typer's built-in test runner to invoke the commands and check for the expected outcomes, such as the successful creation of a database file or a final potential file.

## 5. Test Strategy

The test strategy for CYCLE01 is designed to be comprehensive, ensuring the reliability of each individual component in isolation and the integrity of the entire data pipeline when the components are connected. The strategy is divided into two main pillars: unit testing and integration testing.

**Unit Testing Approach:**

The unit testing approach for this cycle is centered on the principle of strict isolation. Each class, and each public method within it, will be tested independently of its collaborators and external dependencies. The `pytest` framework, in combination with the `pytest-mock` library, will be our primary tool for achieving this. For the `config.py` Pydantic models, we will write a suite of tests that attempt to create instances with both valid and invalid data. For example, we'll assert that providing a negative value for `r_cut` to the `MLIPSettings` model correctly and predictably raises a `pydantic.ValidationError`. This ensures that our data contracts are rigorously enforced at the application's boundaries, preventing invalid data from propagating through the system.

The `AseDBWrapper` will be tested against a temporary, in-memory SQLite database to ensure that our tests are fast and do not have any side effects on the filesystem. We will test the full data lifecycle: we will add an `ase.Atoms` object, read it back, and verify that all its properties (atomic numbers, positions, cell) are perfectly intact. We will also test our ability to update its metadata and, finally, to delete it.

The most critical and detailed unit tests will be for the `LabelingEngine`. We will adopt a "mockist" approach to its testing. The external dependency on the Quantum Espresso binary will be completely mocked out using `mocker.patch('subprocess.run')`. This allows us to test the engine's internal logic in a controlled and reproducible environment. We will have a primary "happy path" test case where the mocked `subprocess.run` returns successfully. We will provide a pre-saved, valid Quantum Espresso output file as a string. The test will then assert that the engine's parser correctly extracts the floating-point numbers for energy, forces, and stress from this sample file. Another critical set of tests will focus on failure modes. One test will simulate a Quantum Espresso crash by having the mock return a non-zero exit code, and we will assert that our `LabelingEngine` catches this and raises a custom, informative exception, such as `DFTCalculationError`. This ensures that our error handling is correct and that the application will fail gracefully. We will also thoroughly test the input file generation logic, asserting that for a given set of atoms and configuration parameters, the generated Quantum Espresso input string contains the correct keywords, values, and formatting.

**Integration Testing Approach:**

While unit tests are essential for verifying that the individual components work correctly in isolation, integration tests are necessary to verify that they communicate and work together as intended. The focus of integration testing in CYCLE01 is on the "happy path" of the entire manual data pipeline. We will create a small-scale test that programmatically mimics the exact workflow a user would follow. This test will not mock the core internal components (`Orchestrator`, `LabelingEngine`, `TrainingEngine`, `AseDBWrapper`), but it will strategically mock the heaviest and most environment-dependent external dependency: the actual DFT calculation itself.

The primary integration test will reside in a file like `tests/test_pipeline.py`. This test will perform the following steps:
1.  Define a simple, valid `FullConfig` object, either directly in the test code or by loading it from a dedicated test YAML file.
2.  Use `pytest`'s `tmp_path` fixture to create a temporary directory and a temporary SQLite database for a clean, isolated test run.
3.  Create a simple `ase.Atoms` object in code (e.g., an H2 molecule or a 2-atom Si cell).
4.  Instantiate the `WorkflowOrchestrator` with the test configuration and a path to the temporary database.
5.  Crucially, we will patch the `subprocess.run` call inside the `LabelingEngine` module. This patch will be designed to bypass the real Quantum Espresso calculation but return a valid, pre-canned output string that mimics a successful run. This is a key point: we are not mocking the engine or the orchestrator, only the expensive external call that the engine makes.
6.  Invoke the orchestrator's `run_labeling` method, passing the path to the test structure.
7.  Assert that the temporary database file was created and now contains exactly one entry for the test structure, and that this entry has the corresponding (mocked) energy and forces.
8.  Invoke the orchestrator's `run_training` method.
9.  Assert that a new potential file (e.g., `model.ace`) has been successfully created in the temporary directory.

This single, powerful test verifies the entire chain of communication: CLI -> Orchestrator -> LabelingEngine -> Database -> TrainingEngine. It confirms that data flows correctly between all the components, that file paths and database connections are handled properly, and that the end-to-end process successfully produces the expected final artifact without any errors. This provides a very high degree of confidence in the foundational stability of the core pipeline before we move on to the significantly more complex automation and active learning features of Cycle 2.
