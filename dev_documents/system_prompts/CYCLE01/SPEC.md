# CYCLE 01: SPEC.md - Core Engine

## 1. Summary

Cycle 01 is the foundational stage of the MLIP-AutoPipe project. The primary objective of this cycle is to establish the core, linear workflow responsible for processing pre-existing structural data, labeling it with first-principles calculations, and training a Machine Learning Interatomic Potential (MLIP). This cycle deliberately omits the complexities of structure generation and active learning to focus exclusively on creating a robust and reliable "processing engine." At the end of this cycle, the system will be capable of taking a user-provided set of atomic configurations and an execution configuration file as input, and producing a trained potential as output.

This initial phase involves implementing several key components that form the backbone of the entire application. We will develop the Pydantic-based configuration models to ensure all settings are type-safe and validated. A dedicated database wrapper for the ASE (Atomic Simulation Environment) database will be created to handle all data persistence, ensuring that every structure and its associated calculation results are stored and tracked systematically. The two main functional modules for this cycle are the `LabelingEngine` and the `TrainingEngine`. The `LabelingEngine` will serve as an automated, robust interface to the Quantum Espresso DFT code, handling input file generation, execution, output parsing, and basic error recovery. The `TrainingEngine` will consume the data produced by the `LabelingEngine` and use it to train a potential based on the Atomic Cluster Expansion (ACE) model, implementing the "Delta Learning" methodology for improved accuracy and physical realism. Finally, a `WorkflowOrchestrator` will be implemented to manage the sequential execution of these components, all driven by a simple command-line interface.

## 2. System Architecture

The focus of Cycle 01 is to build the core components and establish the fundamental project structure. The file system will be laid out as described in the main architectural document, with development concentrated on the files essential for the labeling and training workflow.

**File Structure (Cycle 01 Focus):**

```
.
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── **cli.py**              # Main entry point (Click-based CLI)
│       ├── **configs/**
│       │   ├── __init__.py
│       │   └── **models.py**       # Pydantic models for configuration files
│       ├── **data/**
│       │   ├── __init__.py
│       │   └── **database.py**     # ASE Database wrapper and data management
│       ├── **modules/**
│       │   ├── __init__.py
│       │   ├── a_structure_generator.py
│       │   ├── b_explorer_sampler.py
│       │   ├── **c_labeling_engine.py**
│       │   ├── **d_training_engine.py**
│       │   └── e_simulation_engine.py
│       ├── **utils/**
│       │   ├── __init__.py
│       │   └── **dft_utils.py**    # Helpers for DFT input/output parsing
│       └── **workflow.py**         # The Workflow Orchestrator logic
└── pyproject.toml
```

The files marked in **bold** are the primary deliverables for this cycle. The files for modules A, B, and E will be created as empty placeholders to establish the full architecture but will contain no functional code until later cycles. The `utils` directory is introduced for housing helper functions, such as the detailed parsers required for Quantum Espresso's output format, ensuring a clean separation of concerns from the main engine logic.

## 3. Design Architecture

This cycle establishes the Pydantic-based schema-driven design, which is central to the project's robustness. All configurations and data contracts will be explicitly defined in Pydantic models.

**Pydantic Schema (`configs/models.py`):**

The configuration will be defined by a series of nested Pydantic models to ensure clarity and validation.

*   `DFTComputeConfig`:
    *   `code: str` (e.g., "quantum_espresso")
    *   `command: str` (e.g., "mpirun -np 4 pw.x")
    *   `pseudopotentials: str` (e.g., "SSSP_1.3_PBE_precision")
    *   `ecutwfc: float` (Wavefunction cutoff energy)
    *   `ecutrho: float` (Density cutoff energy)
    *   `kpoints_density: float`
    *   `smearing: str` (e.g., "mv")
    *   `degauss: float`

*   `MLIPTrainingConfig`:
    *   `model_type: str` (e.g., "ace")
    *   `r_cut: float` (Cutoff radius for the potential)
    *   `delta_learning: bool`
    *   `base_potential: str` (e.g., "lj_auto")
    *   `loss_weights: Dict[str, float]` (e.g., `{'energy': 1.0, 'force': 100.0}`)

*   `MainConfig`:
    *   `system: Dict[str, Any]` (e.g., `{'elements': ['Si']}`)
    *   `dft_compute: DFTComputeConfig`
    *   `mlip_training: MLIPTrainingConfig`

**Class and Module Design:**

*   **`AseDBWrapper` (`data/database.py`):**
    *   This class will abstract away the low-level `ase.db` API.
    *   It will provide methods like `connect(path)`, `add_atoms(atoms_list)`, `get_unlabeled_rows()`, and `update_row(row_id, data, key_value_pairs)`.
    *   Consumers of this class (the engines) will not need to know the underlying database implementation details. They will interact with `Atoms` objects and simple dictionaries.

*   **`LabelingEngine` (`modules/c_labeling_engine.py`):**
    *   The constructor will accept a `DFTComputeConfig` object and an `AseDBWrapper` instance (Dependency Injection).
    *   Its main public method, `run()`, will fetch unlabeled structures from the database.
    *   For each structure, it will:
        1.  Use helper functions in `ase.io` and `dft_utils.py` to generate a Quantum Espresso input file string from the `Atoms` object and `DFTComputeConfig`.
        2.  Execute `pw.x` using `subprocess.run`. It must capture stdout and stderr for parsing and error checking.
        3.  Call a parser function (from `dft_utils.py`) to extract the final energy, forces, and stress from the stdout string.
        4.  Update the corresponding row in the database with the results and set its status to 'labeled'.
    *   It will include basic error handling to detect SCF convergence failures from the output.

*   **`TrainingEngine` (`modules/d_training_engine.py`):**
    *   The constructor will accept an `MLIPTrainingConfig` and an `AseDBWrapper` instance.
    *   Its main public method, `run()`, will fetch all labeled structures from the database.
    *   It will implement the Delta Learning logic: if `delta_learning` is true, it will calculate the energy/forces from the specified `base_potential` and subtract them from the DFT values.
    *   It will then prepare the data in the format required by the ACE training library and invoke the training process.
    *   After training, it will save the final model artifact (e.g., `model.ace`) to the working directory.

## 4. Implementation Approach

The implementation will proceed in a logical, bottom-up fashion, starting with data structures and moving to functional components.

1.  **Project Setup:** Initialize the project with `uv` and create the directory structure as outlined above. Add initial dependencies (`click`, `pydantic`, `ase`, `pyyaml`) to `pyproject.toml`.

2.  **Configuration Models:** Implement all Pydantic models in `src/mlip_autopipec/configs/models.py`. Include validation logic where necessary (e.g., ensuring loss weights are positive).

3.  **Database Wrapper:** Implement the `AseDBWrapper` in `src/mlip_autopipec/data/database.py`. Focus on creating a clean and simple API for the rest of the application.

4.  **DFT Utilities:** In `src/mlip_autopipec/utils/dft_utils.py`, write the core function for parsing Quantum Espresso output files. This function should be robust and able to handle slight variations in output format. Use regular expressions to find and extract the total energy, and the blocks for atomic forces and stress.

5.  **Labeling Engine:** Implement the `LabelingEngine` class. Inject its dependencies. Write the logic to generate input files using ASE's built-in writers, run the `subprocess` command, and call the parsing utility. The focus is on robustly handling the external process call.

6.  **Training Engine:** Implement the `TrainingEngine` class. For this initial cycle, the interaction with the actual ACE training library can be simplified. The main focus is on correctly retrieving data from the database, processing it according to the Delta Learning configuration, and ensuring the data flow is correct.

7.  **Orchestrator and CLI:** Implement the `WorkflowOrchestrator` in `src/mlip_autopipec/workflow.py`. This class will instantiate the database wrapper and the two engines, passing the relevant configuration objects and the database instance to them. It will call `labeling_engine.run()` followed by `training_engine.run()`. Finally, implement the `cli.py` file using `click` to parse the command-line arguments, load the YAML configuration, instantiate the orchestrator, and run the workflow.

## 5. Test Strategy

Testing in Cycle 01 is crucial for building a reliable foundation. The strategy is divided into isolated unit tests and tests that verify the connection between components.

**Unit Testing Approach:**
(Located in `tests/unit/`)
*   **`configs/models.py`:** Test Pydantic model validation. For instance, provide invalid data (e.g., a negative loss weight) and assert that a `ValidationError` is raised.
*   **`data/database.py`:** Test the `AseDBWrapper`. Use `pytest-mock` to patch `ase.db.connect` and `ase.db.Connection` objects. This allows testing the wrapper's logic (e.g., how it formats queries) without any actual file I/O, making the tests extremely fast and reliable.
*   **`utils/dft_utils.py`:** Create a dedicated test file with sample Quantum Espresso output stored as a string multiline literal. The test will call the parsing function with this string and assert that the correct floating-point values for energy, forces, and stress are extracted. This is a critical test for the robustness of the `LabelingEngine`.
*   **`modules/c_labeling_engine.py`:** Test the `LabelingEngine` class in isolation. Mock the `AseDBWrapper` dependency to provide a known `Atoms` object. Mock the `subprocess.run` call to return a `CompletedProcess` object containing the sample QE output string from the previous test. The test will then assert that the engine calls the `db.update_row` method with the correctly parsed data.
*   **`modules/d_training_engine.py`:** Test the `TrainingEngine`. Mock the `AseDBWrapper` to provide a list of labeled `Atoms` objects. The test will verify that the Delta Learning logic is applied correctly (or skipped if disabled). The actual call to the ACE library's training function will be mocked to avoid a lengthy training process, asserting only that it was called with correctly prepared data.

**Integration Testing Approach:**
(Located in `tests/e2e/`)
*   **Orchestrator Test:** A single, powerful integration test will be created to verify the entire Cycle 01 workflow.
    *   **Setup:** The test will create a temporary directory. It will programmatically create a temporary ASE database and populate it with one or two unlabeled `ase.Atoms` objects (e.g., a slightly distorted Si atom). It will also create a valid `exec_config_dump.yaml` in the directory.
    *   **Execution:** The test will mock `subprocess.run` to prevent actual DFT execution, making it return a valid, pre-computed output string whenever it's called. It will also mock the final ACE training call in the `TrainingEngine`. The test will then instantiate and run the `WorkflowOrchestrator`.
    *   **Verification:** After the orchestrator finishes, the test will connect to the temporary database and assert that the status of the atoms has been changed to 'labeled' and that the energy/force data has been correctly written. It will also assert that the mocked training function was called, confirming that the data flowed correctly from one engine to the next.
