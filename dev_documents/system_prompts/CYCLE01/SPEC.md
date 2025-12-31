# CYCLE01/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 1 of the MLIP-AutoPipe project. The primary objective of this foundational cycle is to construct the core engine of the entire pipeline, establishing a robust workflow for automated DFT calculations and subsequent MLIP training. This cycle is critical as it lays the architectural groundwork, project structure, and essential tooling upon which all future functionality will be built. The core philosophy of "removing the human expert from the loop" begins here, with a focus on creating a system that can programmatically take a set of atomic structures, calculate their physical properties using a powerful DFT engine (Quantum Espresso), and then use this high-fidelity data to train a Machine Learning Interatomic Potential (MLIP) based on the Atomic Cluster Expansion (ACE) framework.

Success in this cycle means delivering a functional, command-line-driven tool that, while not yet fully autonomous in its data generation, is completely automated in its data processing and model training phases. We will establish the modern Python project structure using `pyproject.toml` and the `uv` package manager, ensuring a reproducible and high-performance development environment. The two key modules to be delivered are the **Labeling Engine (Module C)**, which will act as a sophisticated wrapper around Quantum Espresso, and the **Training Engine (Module D)**, which will implement the delta learning methodology with the ACE model. Data persistence will be handled by a dedicated `AseDBWrapper`, ensuring that all computed data is stored in a structured, queryable format. The orchestrator will be kept simple, responsible only for executing the linear `Label -> Train` workflow. By the end of this cycle, we will have a tangible, working system that proves the viability of the core components and provides a solid foundation for the more advanced features planned in subsequent cycles, such as automated structure generation and active learning.

## 2. System Architecture

The system architecture for Cycle 1 is focused on establishing the core components and their interactions. The file structure is designed to be modular and scalable, clearly separating concerns between the command-line interface, data management, orchestration, and the individual processing modules.

**File Structure:**

The following ASCII tree shows the files to be created or modified in this cycle. Bolded files are the primary focus of implementation for Cycle 1.

```
.
├── dev_documents/
│   └── ...
├── src/
│   └── mlip_autopipec/
│       ├── **__init__.py**
│       ├── **cli.py**          # Main entry point (Click-based CLI)
│       ├── data/
│       │   ├── **__init__.py**
│       │   ├── **database.py**   # AseDB wrapper class
│       │   └── **models.py**     # Pydantic data models for config
│       ├── modules/
│       │   ├── **__init__.py**
│       │   ├── **c_labeling_engine.py**
│       │   └── **d_training_engine.py**
│       ├── **orchestrator.py** # Main workflow controller
│       └── utils/
│           ├── **__init__.py**
│           └── **dft_utils.py**  # Helpers for QE integration
├── tests/
│   ├── unit/
│   │   ├── data/
│   │   │   └── **test_database.py**
│   │   ├── modules/
│   │   │   └── **test_c_labeling_engine.py**
│   │   │   └── **test_d_training_engine.py**
│   │   └── utils/
│   │       └── **test_dft_utils.py**
│   └── e2e/
│       └── **test_cycle01_workflow.py**
├── **pyproject.toml**          # Project definition and dependencies
```

**Component Blueprint:**

*   **`pyproject.toml`**: This file will define all project metadata, dependencies (e.g., `click`, `ase`, `pydantic`, `mace-torch`), and tool configurations (e.g., `pytest`, `ruff`). It is the single source of truth for setting up the project's environment using `uv`.
*   **`cli.py`**: This will implement the main user-facing interface using `click`. A command `mlip-pipe run-cycle1` will be created to trigger the workflow. It will handle parsing command-line arguments, such as the path to the configuration file and the ASE database.
*   **`orchestrator.py`**: The `Orchestrator` class will be the central coordinator. Its constructor will accept the configuration and a database path. A public method, `run_label_and_train_workflow()`, will execute the Cycle 1 pipeline by sequentially initializing and calling the `LabelingEngine` and `TrainingEngine`.
*   **`data/database.py`**: The `AseDBWrapper` class will abstract all interactions with the `ase.db` SQLite database. It will provide methods like `connect()`, `disconnect()`, `get_rows_to_label()`, and `update_row_with_dft_results()`. This encapsulation is crucial for ensuring database transactions are handled correctly and for simplifying testing.
*   **`data/models.py`**: Pydantic models will be defined here to represent the structure of the configuration file. This ensures that all settings are type-checked and validated at runtime. `DFTConfig` will model settings for Quantum Espresso (e.g., `ecutwfc`, `kpoints_density`), and `TrainingConfig` will model settings for the ACE model (e.g., `r_cut`, `delta_learning`).
*   **`modules/c_labeling_engine.py`**: The `LabelingEngine` class will be responsible for the entire DFT calculation process. It will query the database for unlabeled structures, and for each one, it will:
    1.  Generate a valid Quantum Espresso input file using helper functions from `dft_utils.py`.
    2.  Execute `pw.x` as a subprocess.
    3.  Monitor the subprocess and handle potential errors or timeouts.
    4.  Parse the output file to extract forces, energy, and stress.
    5.  Update the corresponding record in the database with the results.
*   **`modules/d_training_engine.py`**: The `TrainingEngine` class will handle MLIP training. It will query the database for all successfully labeled structures, convert them into the format expected by the training library (e.g., a list of ASE `Atoms` objects), configure the ACE model and trainer based on the settings, and execute the training process. It will also implement the "delta learning" logic by calculating and subtracting the baseline potential's contribution before training.
*   **`utils/dft_utils.py`**: This module will contain pure, stateless functions for handling Quantum Espresso's file formats. Key functions will include `create_qe_input_from_atoms()` and `parse_qe_output()`. This separation isolates the complex logic of file format conversion and parsing, making it independently testable and reusable.

## 3. Design Architecture

The design architecture for Cycle 1 is centered around creating a robust, type-safe, and decoupled system. Pydantic-based schemas are the cornerstone of this design, ensuring that data, especially configuration data, is strictly validated and clearly defined throughout the application.

**Pydantic Schema Design:**

The configuration system will be governed by a hierarchy of Pydantic models defined in `src/mlip_autopipec/data/models.py`.

*   **`DFTCompute` (Model)**: This model will encapsulate all parameters required by the `LabelingEngine`.
    *   `code`: A string literal, initially fixed to `"quantum_espresso"`.
    *   `command`: A string for the execution command (e.g., `"mpirun -np 32 pw.x"`).
    *   `pseudopotentials`: A string representing the selected protocol (e.g., `"SSSP_1.3_PBE_precision"`).
    *   `ecutwfc`, `ecutrho`: Floats for wavefunction and density cutoffs.
    *   `kpoints_density`: A float defining the k-point mesh density.
    *   *Invariants*: Validation rules will ensure that `ecutrho` is always greater than or equal to `ecutwfc`.

*   **`MLIPTraining` (Model)**: This model will contain all settings for the `TrainingEngine`.
    *   `model_type`: A string literal, fixed to `"ace"` for this cycle.
    *   `r_cut`: A float for the model's cutoff radius.
    *   `delta_learning`: A boolean to enable or disable the delta learning strategy.
    *   `base_potential`: An optional string specifying the baseline potential (e.g., `"lj_auto"`).
    *   `loss_weights`: A nested model or dictionary mapping "energy", "force", and "stress" to their respective float weights in the loss function.
    *   *Invariants*: A validator will ensure that if `delta_learning` is `True`, `base_potential` must not be `None`.

*   **`Cycle01Config` (Top-level Model)**: This model will compose the other models into a single configuration object.
    *   `dft_compute`: An instance of the `DFTCompute` model.
    *   `mlip_training`: An instance of the `MLIPTraining` model.
    *   `database_path`: A `FilePath` type from Pydantic, ensuring the file path exists.

**Data Flow and Consumers/Producers:**

*   **Producer (`cli.py`)**: The CLI will read a YAML configuration file from disk. It will then use Pydantic's parsing capabilities to instantiate the `Cycle01Config` model. This validated model object is the primary data passed to the `Orchestrator`.
*   **Consumer (`Orchestrator`)**: The `Orchestrator` consumes the `Cycle01Config` object. It acts as a router, passing the `dft_compute` part of the config to the `LabelingEngine` and the `mlip_training` part to the `TrainingEngine`.
*   **Consumer (`LabelingEngine`, `TrainingEngine`)**: These modules are the ultimate consumers of the configuration models. They use the validated parameters within these models to configure their internal logic, such as setting up the DFT calculation or initializing the ACE model trainer.

This Pydantic-driven design ensures that any invalid configuration is caught at the earliest possible moment (at startup), preventing runtime errors deep within the calculation or training loops. It also serves as a form of self-documentation, as the Pydantic models provide a clear and unambiguous definition of all available settings. This architecture is designed for future extensibility; adding new configuration options will be as simple as adding a new field to a Pydantic model.

## 4. Implementation Approach

The implementation will proceed in a logical, step-by-step manner, starting with the foundational data layer and progressively building up to the full workflow.

1.  **Project Setup**:
    *   Initialize the project directory structure as outlined in the architecture.
    *   Create the `pyproject.toml` file, defining the project name (`mlip-autopipec`) and adding initial dependencies: `click`, `ase`, `pydantic`, `pyyaml`, `mace-torch`, and `pytest`.
    *   Set up a virtual environment using `uv venv` and install the dependencies with `uv pip install -e .`.

2.  **Data Layer Implementation**:
    *   Define the Pydantic configuration models (`DFTCompute`, `MLIPTraining`, `Cycle01Config`) in `src/mlip_autopipec/data/models.py`.
    *   Implement the `AseDBWrapper` class in `src/mlip_autopipec/data/database.py`. Start with essential methods: `connect`, `disconnect`, `add_atoms`, and `get_atoms_by_id`. Use a test-driven approach, writing a corresponding test in `tests/unit/data/test_database.py` for each method.

3.  **Quantum Espresso Utilities**:
    *   In `src/mlip_autopipec/utils/dft_utils.py`, write the `create_qe_input_from_atoms` function. This function will take an ASE `Atoms` object and a `DFTCompute` config object and return a string containing a valid QE input file. This should be tested thoroughly, ensuring all parameters are correctly formatted.
    *   Next, implement the `parse_qe_output` function. This will take a string (the content of a QE output file) and extract the final energy, forces, and stress. Write unit tests in `tests/unit/utils/test_dft_utils.py` with example output files to cover success and failure cases.

4.  **Labeling Engine Implementation**:
    *   Create the `LabelingEngine` class in `src/mlip_autopipec/modules/c_labeling_engine.py`.
    *   The `execute` method will be the main entry point. It will use the `AseDBWrapper` to fetch structures that need labeling.
    *   For each structure, it will call the `dft_utils` functions to generate the input, then use Python's `subprocess` module to run `pw.x`. **Crucially, `shell=False` must be used for security and robustness.**
    *   After the process completes, it will read the output file, call the parsing utility, and update the database with the results via the `AseDBWrapper`.
    *   Unit tests for this module will mock the `subprocess.run` call and the `AseDBWrapper` to test the logic in isolation.

5.  **Training Engine Implementation**:
    *   Create the `TrainingEngine` class in `src/mlip_autopipec/modules/d_training_engine.py`.
    *   The `execute` method will fetch all labeled data from the database.
    *   It will implement the delta learning logic: for each structure, calculate the energy/forces from a simple reference potential (e.g., Lennard-Jones) and subtract this from the DFT values.
    *   It will then configure and run the ACE model training process using the `mace-torch` library. The trained model will be saved to a file.
    *   Unit tests will use a small, fixed dataset and mock the `mace-torch` training function to verify that the data preparation and delta learning logic are correct.

6.  **Orchestration and CLI**:
    *   Implement the `Orchestrator` class to tie everything together.
    *   Implement the `cli.py` module to create the command-line interface.
    *   Finally, write an end-to-end test in `tests/e2e/test_cycle01_workflow.py` that simulates a user running the CLI command. This test will use a real (but small and fast) YAML config file and a temporary database, and will mock the `subprocess.run` call to avoid running a real DFT calculation, instead returning a pre-canned output file. This verifies that the entire pipeline is connected correctly.

## 5. Test Strategy

The test strategy for Cycle 1 is designed to build confidence in the core components of the system from the ground up, ensuring that each part is reliable before it is integrated into the whole.

**Unit Testing Approach (Min 300 words):**
The unit testing strategy focuses on isolating each component and verifying its logic independently. We will use `pytest` as our test runner and `unittest.mock` for creating mock objects.
*   **`AseDBWrapper`**: The tests in `test_database.py` will validate the correctness of all database interactions. A temporary SQLite database file will be created during test setup and torn down afterwards. We will test a full CRUD (Create, Read, Update, Delete) cycle: adding atoms, retrieving them, updating their metadata (e.g., marking them as 'labeled'), and ensuring data integrity. This ensures that the foundation of our data persistence layer is solid.
*   **`dft_utils`**: In `test_dft_utils.py`, we will test the input file generation and output parsing logic exhaustively. For input generation, we will create various ASE `Atoms` objects (with different elements, cell sizes, and magnetic properties) and assert that the generated Quantum Espresso input string contains the correct keywords and values. For parsing, we will have static text files containing examples of QE output (including successful runs and runs that failed to converge). The tests will assert that our parsing function correctly extracts the required data in success cases and returns an appropriate error or `None` in failure cases.
*   **`LabelingEngine`**: The tests in `test_c_labeling_engine.py` will verify the orchestration logic of the DFT calculations. The key here is to mock external interactions. We will mock the `AseDBWrapper` to provide a controlled list of 'unlabeled' atoms. We will also mock `subprocess.run` to simulate the execution of `pw.x` without actually running it. The test will assert that the engine calls the subprocess with the correct command and that it correctly processes the (mocked) standard output and updates the (mocked) database.
*   **`TrainingEngine`**: In `test_d_training_engine.py`, we will focus on the data preparation logic. A test will create a small, fixed dataset of labeled `Atoms` objects. We will verify that the delta learning calculation is performed correctly by manually calculating the expected residual values and asserting that the engine's output matches. The actual call to the `mace-torch` training function will be mocked to avoid a lengthy training process, and we will simply assert that it was called with the correctly prepared data and configuration.

**Integration Testing Approach (Min 300 words):**
The integration testing for Cycle 1 will focus on verifying the data flow and interactions between the major components: CLI -> Orchestrator -> Engines -> Database.
*   **Orchestrator-Engine Integration**: We will write tests that initialize the `Orchestrator` with mock `LabelingEngine` and `TrainingEngine` objects. The purpose is to verify that the orchestrator calls the `execute` method of each engine in the correct sequence (`LabelingEngine` first, then `TrainingEngine`). This test ensures the high-level control flow is correct.
*   **Engine-Database Integration**: We will test the interaction between the engines and a real (temporary) database. For the `LabelingEngine`, the test will pre-populate a temporary `AseDB` with a few unlabeled structures. It will then run the engine (with `subprocess.run` mocked) and assert that after execution, the rows in the database have been correctly updated with the parsed DFT results. Similarly, for the `TrainingEngine`, the test will pre-populate the database with labeled structures, run the engine, and assert that it correctly queries and processes all available data. This confirms that the SQL queries and data models used by the engines are compatible with the actual database schema.
*   **End-to-End Workflow Test**: The primary integration test will be `tests/e2e/test_cycle01_workflow.py`. This will be a comprehensive test that simulates a user's action. It will use the `click.testing.CliRunner` to invoke the `mlip-pipe run-cycle1` command from within the test. A small, cycle-specific YAML configuration file will be created on the fly. The test will mock the `subprocess.run` call to the DFT code to ensure speed and determinism, providing a fixture of a successful QE output. The test will then assert on the final state: that the database contains correctly labeled structures and that a model file has been created by the `TrainingEngine`. This test serves as the ultimate validation that all components are wired together correctly and that the data flows seamlessly from configuration file to final trained model.
