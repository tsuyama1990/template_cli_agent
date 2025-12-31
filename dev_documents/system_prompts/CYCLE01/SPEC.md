# SPEC.md: Cycle 01 - The Core Engine

## 1. Summary

Cycle 01 represents the foundational pillar of the MLIP-AutoPipe project, a crucial first step that prioritises stability, reliability, and architectural soundness over feature richness. The primary objective of this cycle is to construct the essential, non-interactive backbone of the entire pipeline. This involves developing and integrating the core components required to perform the most fundamental task: taking a single, predefined atomic structure, accurately calculating its physical properties (energy, forces, and stress) using a first-principles Density Functional Theory (DFT) engine, and subsequently training a functional Machine Learning Interatomic Potential (MLIP) from this single data point. To maintain focus and mitigate risk, this cycle deliberately and completely omits user-facing features such as configuration files, automated structure generation, and the complexities of active learning. The exclusive focus is on establishing a robust, verifiable, and architecturally clean workflow that can serve as a solid foundation for all subsequent development cycles.

The two principal software components to be designed, implemented, and tested in this cycle are the **Labelling Engine (Module C)** and the **Training Engine (Module D)**. The Labelling Engine will be far more than a simple script; it will be a sophisticated and resilient wrapper around the Quantum Espresso (QE) software package. Its purpose is to abstract away the significant complexities and potential pitfalls of setting up and running DFT calculations. A key design requirement is robustness, which will be achieved by incorporating automated error detection and a multi-step recovery strategy. This will enable the engine to handle common calculation failures, such as lack of self-consistent field (SCF) convergence, without requiring manual intervention. The Training Engine will serve as a clean, high-level interface to the chosen MLIP framework (initially, the Atomic Cluster Expansion or ACE model). A critical feature of this engine will be the implementation of the "delta learning" technique. This method, which involves training the model on the *residual* between a simple, computationally cheap physical potential (like Lennard-Jones) and the true DFT values, is of paramount importance for ensuring the physical plausibility and numerical stability of the final MLIP, especially in high-energy or short-range configurations.

To ensure data integrity and provenance from the very beginning, all data persistence will be handled by a dedicated `AseDB` wrapper class. This class will manage all interactions with an underlying SQLite database, providing a structured and reliable way to store, query, and manage the state of atomic structures and their associated DFT results. By the successful conclusion of Cycle 01, the project will possess a functional, albeit minimal, end-to-end workflow. This will not only prove the viability of the core architectural concept—the tight integration of automated DFT labelling with MLIP training—but will also provide a stable, testable platform upon which the more complex features of subsequent cycles can be confidently built.

## 2. System Architecture

For Cycle 01, the system architecture is deliberately and strategically streamlined to focus exclusively on the core data processing pipeline. The workflow is entirely linear, non-interactive, and deterministic, designed to prove the fundamental mechanics of the system. It is initiated not by a user, but by a developer-run script that specifies a single, hardcoded input structure. This minimalist approach allows for focused development and testing of the most critical components: the Orchestrator, the Labelling Engine, the Training Engine, and the Database interface. The interaction between these components is simple and direct, forming a single, clear data flow path that is easy to debug and verify.

The architectural workflow is as follows:
1.  **Initiation**: The process begins when a developer executes a test script. This script instantiates the `Orchestrator` and invokes its main workflow method, passing it a predefined ASE Atoms object (e.g., a simple molecule or a small crystal unit cell). This object is the sole input to the entire system for this cycle.
2.  **Labelling Request and State Management**: The Orchestrator immediately communicates with the `AseDB` wrapper to save the incoming structure and assign it an initial state, such as 'awaiting_labelling'. It then passes the Atoms object to the `Labelling Engine` for processing.
3.  **DFT Calculation**: The Labelling Engine encapsulates all the logic for interacting with Quantum Espresso. It dynamically generates a QE input file from the Atoms object's properties, executes the `pw.x` binary as a managed subprocess, and carefully monitors the execution. Upon completion, it parses the voluminous text output to extract the key physical quantities. This module also contains the crucial error handling logic to gracefully manage issues like SCF non-convergence.
4.  **Data Persistence (Labels)**: The extracted energy, forces, and stress are packaged into a structured `DFTResult` data object. The Orchestrator receives this object and passes it to the `AseDB`, which persists the results and updates the status of the original structure to 'labelled_success' or 'labelled_failure'.
5.  **Training Data Retrieval**: The Orchestrator then queries the `AseDB` for a complete batch of training data. In this simplified cycle, this batch will consist of the single data point that was just generated and stored.
6.  **MLIP Training**: The retrieved data pair (Atoms object and DFTResult) is passed to the `Training Engine`.
7.  **Delta Learning Calculation**: The Training Engine first calculates the energy and forces for the Atoms object using a simple, hardcoded baseline physical potential (e.g., Lennard-Jones). It then computes the residual by subtracting these baseline values from the DFT ground-truth values.
8.  **Model Training Execution**: The engine then uses the chosen MLIP library (ACE) to train a new model, using the calculated residuals as the training target.
9.  **Data Persistence (Model)**: The resulting trained model is serialised and saved to a designated file in the filesystem. The Orchestrator also instructs the `AseDB` to store metadata about this model, such as its creation timestamp and a link to the training data used to generate it.

This architectural design enforces a strong separation of concerns, which is critical for a maintainable and extensible system. The Orchestrator is the "brain," concerned only with the high-level sequence of operations. The Engines are the "muscles," containing the detailed, complex implementation logic for their specific domains (DFT calculation and MLIP training). The database acts as a passive, transactional data store, completely decoupling the data generation process (by the Labelling Engine) from the data consumption process (by the Training Engine). This modularity is a key strategic investment, as it will allow for the seamless integration of more sophisticated orchestration logic, such as the active learning loops planned for later cycles, without requiring any disruptive modifications to the core, validated engines developed in this foundational cycle.

## 3. Design Architecture

This cycle will establish the foundational code structure within the `src/mlip_autopipec/` directory, which will serve as the blueprint for all future development. The design strongly emphasises object-oriented principles, including encapsulation, clear public APIs, and the use of data models (Pydantic) to ensure type safety and data integrity throughout the application. The goal is to create a codebase that is not only functional but also clean, readable, and highly testable.

**Key Classes and APIs:**

1.  **`ac_cdd_core.services.AseDB`** (Located in `data/database.py`)
    *   **Purpose**: To provide a high-level, abstracted interface for all database operations, insulating the rest of the application from the specific syntax of the underlying `ase.db` library.
    *   **Public API**:
        *   `__init__(db_path: str)`: Connects to the specified SQLite database file, creating it if it does not exist.
        *   `add_atoms(atoms: ase.Atoms, state: str = 'initial') -> int`: Adds a new ASE Atoms object to the database. It returns the unique integer ID assigned to the new row. The `state` parameter is crucial for tracking the structure's progress through the workflow.
        *   `get_atoms_to_label(limit: int = 10) -> List[ase.Atoms]`: Fetches a list of structures from the database that are in a state indicating they need to be processed by the Labelling Engine.
        *   `write_dft_result(atoms_id: int, result: DFTResult)`: Writes the complete outcome of a DFT calculation (the `DFTResult` object) to the database. This method is responsible for linking the result to the original atoms row and updating that row's state to 'labelled'.
        *   `update_state(atoms_id: int, new_state: str)`: Provides a generic mechanism to update the state of a structure, which will be useful for marking failures or other statuses.
        *   `get_training_data() -> Tuple[List[ase.Atoms], List[DFTResult]]`: Retrieves a complete dataset of all successfully labelled atoms and their corresponding results, formatted conveniently for the Training Engine.

2.  **`ac_cdd_core.domain_models.DFTResult`** (Located in `data/models.py`)
    *   **Purpose**: To serve as a robust, type-safe data transfer object for the results of a DFT calculation. Using a Pydantic model here prevents common errors associated with passing around unstructured dictionaries.
    *   **Attributes**:
        *   `energy: float`: The total potential energy calculated by DFT.
        *   `forces: np.ndarray`: A NumPy array of shape (N, 3), where N is the number of atoms, representing the forces on each atom.
        *   `stress: np.ndarray`: A NumPy array of shape (6,) representing the Voigt-ordered stress tensor.
        *   `status: Literal['success', 'failed']`: A mandatory status flag to clearly indicate the outcome of the calculation.
        *   `metadata: Dict`: An open-ended dictionary to store supplementary information, such as DFT wall time, convergence parameters, or error messages.

3.  **`mlip_autopipec.modules.LabellingEngine`** (Located in `modules/labelling_engine.py`)
    *   **Purpose**: To encapsulate and manage the entire lifecycle of running a Quantum Espresso calculation.
    *   **Public API**:
        *   `__init__(qe_command: str, pseudo_dir: str)`: The constructor is initialised with the system-specific command needed to run QE (e.g., "mpirun -np 4 pw.x") and the path to the directory containing pseudopotential files.
        *   `run(atoms: ase.Atoms) -> DFTResult`: The primary public method. It orchestrates the entire process for a single structure and returns a comprehensive `DFTResult` object.
    *   **Key Internal Methods**:
        *   `_generate_input_file(atoms: ase.Atoms) -> str`: A private method that performs the critical translation from an `ase.Atoms` object into the highly specific, formatted string required for a QE input file.
        *   `_execute_qe(input_content: str) -> Tuple[str, str, int]`: A private method that handles the execution of the `pw.x` binary in a secure subprocess, safely capturing its stdout, stderr, and return code for later analysis.
        *   `_parse_output(output_content: str) -> Dict`: A private method containing robust regular expressions and parsing logic to extract the required physical quantities from the QE text output.
        *   `_handle_error(atoms: ase.Atoms, stdout: str, stderr: str) -> DFTResult`: A private method that implements the recovery logic. For Cycle 01, it will focus on diagnosing common errors and returning a `DFTResult` with a 'failed' status and informative metadata.

4.  **`mlip_autopipec.modules.TrainingEngine`** (Located in `modules/training_engine.py`)
    *   **Purpose**: To abstract the complexities of the underlying MLIP framework and implement the delta learning logic.
    *   **Public API**:
        *   `__init__(model_type: str = 'ace')`: Initialises the trainer for a specific model type.
        *   `train(data: List[Tuple[ase.Atoms, DFTResult]]) -> Any`: The main public method. It takes the full labelled dataset and returns a trained, serialisable model object.
    *   **Key Internal Methods**:
        *   `_get_baseline_potential(elements: List[str]) -> ase.calculators.Calculator`: A private method that returns a simple, physics-based ASE calculator (e.g., `ase.calculators.lj.LennardJones`) to be used for the delta learning baseline.
        *   `_compute_delta(...) -> Dict`: A private method that calculates the residual between the DFT ground truth and the baseline potential for energy, forces, and stress.
        *   `_prepare_training_dataset(...)`: A private method that takes the calculated deltas and transforms them into the specific data format required by the underlying ACE training library's API.

This detailed, object-oriented design ensures that each class has a single, well-defined responsibility, which is paramount for building a system that is testable, maintainable, and ready for the addition of more complex features in subsequent development cycles.

## 4. Implementation Approach

The implementation for Cycle 01 will follow a logical, bottom-up progression. We will start with the fundamental data structures and the data persistence layer, then build the core processing engines, and finally tie them together with a simple orchestrator. This layered approach ensures that each component can be independently unit-tested before being integrated, which significantly simplifies debugging and validation.

1.  **Project Setup and Environment**:
    *   The first step is to create the complete project directory structure as outlined in the Design Architecture section. This includes creating all the necessary sub-packages (`config`, `data`, `modules`, `utils`) with their `__init__.py` files.
    *   Next, the `pyproject.toml` file will be created. This file is central to the project's dependency management and packaging. It will be populated with the initial set of dependencies required for this cycle, including `ase`, `numpy`, `pydantic`, `click`, and the specific MLIP training library that will be used (e.g., `mace-torch`).
    *   The development environment will be set up using `uv venv` to create a virtual environment, and the project will be installed in editable mode (`uv pip install -e .`) to ensure that changes to the source code are immediately reflected in the installed package.

2.  **Data Layer Implementation**:
    *   The development of the code itself will begin with the data layer, as it is a dependency for almost all other components. The `DFTResult` Pydantic model will be implemented first in `data/models.py` to define the data contract for DFT calculations.
    *   Following this, the `AseDB` wrapper class will be implemented in `data/database.py`. The initial version will focus on the essential methods required for the Cycle 01 workflow: connecting to a local SQLite file (`mlip.db`), adding new atoms with an initial state, and writing the `DFTResult` upon completion of a calculation.

3.  **Labelling Engine Development (Module C)**:
    *   The `LabellingEngine` is arguably the most complex component of this cycle. Development will start with the `_generate_input_file` method. This is a critical and detailed task that involves correctly translating the properties of an `ase.Atoms` object (cell vectors, atomic positions, chemical symbols) into the strict, column-formatted text required by Quantum Espresso. This includes setting default values for k-points, pseudopotentials, and calculation type (`scf`).
    *   Next, the `_parse_output` method will be implemented. This requires writing robust regular expressions to reliably find and extract the total energy, the block of atomic forces, and the stress tensor from the potentially very long QE output file.
    *   The `_execute_qe` method will then be implemented using Python's `subprocess.run`. Care will be taken to correctly handle the standard input, output, and error streams.
    *   These private methods will be combined into the main public `run` function.
    *   Finally, the initial error handling logic will be developed. This will involve inspecting the captured stderr and stdout for key error messages (like "SCF not converged") and ensuring that a `DFTResult` with a 'failed' status is returned in these cases.

4.  **Training Engine Development (Module D)**:
    *   Development of the `TrainingEngine` will begin with the implementation of the baseline potential logic. This will likely involve a simple heuristic to generate reasonable Lennard-Jones parameters based on atomic properties available in ASE.
    *   The `_compute_delta` method will then be implemented. This is a relatively straightforward step that involves using the baseline calculator to get its energy and forces, and then subtracting them from the DFT values stored in the `DFTResult` object.
    *   The main `train` method will be written last, as it depends on the delta computation. This will involve the specific API calls required by the ACE library to configure a training job, pass the dataset, and execute the training loop.
    *   Finally, the logic for serialising the trained model object to a file (e.g., using `torch.save`) will be implemented.

5.  **Orchestrator and Temporary CLI**:
    *   A very basic `Orchestrator` class will be created in `orchestrator.py`. For this cycle, it will contain a single method, `run_cycle01_workflow(atoms: ase.Atoms)`, which will contain the simple, linear sequence of calls to the database and the engines.
    *   To provide a simple entry point for testing, a minimal `click` command will be created in `main.py`. This command will do nothing more than instantiate the `Orchestrator` and call the workflow method with a hardcoded test case (e.g., a single Argon atom or an H2 molecule), printing a success or failure message to the console. This serves as the final integration point for all the components developed in the cycle.

This methodical, step-by-step approach ensures that each component is built on a solid, previously-tested foundation, which is essential for managing complexity and ensuring the robustness of the final integrated system.

## 5. Test Strategy

The test strategy for Cycle 01 is of paramount importance, as it will validate the core functionality upon which the entire project will be built. The strategy is designed to be comprehensive, focusing on isolated unit tests for each component to ensure its internal logic is correct, and culminating in a single, crucial integration test to verify that the components work together as expected. All tests will be automated and designed to run quickly without any external dependencies like a real DFT code.

**Unit Testing Approach (Min 600 words):**

The unit testing for Cycle 01 will be extremely fine-grained and thorough, adhering to the principle that every logical path in the code should be tested. All external dependencies—particularly the filesystem, subprocesses, and the MLIP training library—will be rigorously mocked using Python's `unittest.mock` library. This ensures that the tests are fast, deterministic, and focus solely on the logic of the code being tested.

*   **`AseDB`**: The database wrapper will be tested against a temporary, in-memory SQLite database (`":memory:"`) to avoid filesystem I/O. The test suite will cover all public methods. For `add_atoms`, we will add several atoms and assert that the returned IDs are unique and incrementing. We will then query the database directly to confirm the data was inserted correctly. For `write_dft_result`, we will first add an atom, then call the method with a mock `DFTResult`, and finally assert that the corresponding row in the database was updated with the correct result data and that its state was correctly changed from 'initial' to 'labelled'. The `get_training_data` method will be tested by populating the database with a mix of labelled and unlabelled atoms and asserting that the method returns only the correctly labelled ones.

*   **`LabellingEngine`**: This component, with its external process interaction, requires the most extensive and careful mocking.
    *   **`_generate_input_file`**: This method's correctness is critical. It will be tested by creating a variety of `ase.Atoms` objects representing different physical systems (e.g., a molecule with no periodic boundaries, a bulk crystal with PBC, systems with different elements) and calling the method for each. The test will then assert that the generated QE input string is character-for-character identical to a pre-validated, "golden" input file stored as a test asset. This ensures that our generated inputs are always valid.
    *   **`_execute_qe`**: The `subprocess.run` call within this method will be completely mocked using `unittest.mock.patch`. The mock will be configured to return different `subprocess.CompletedProcess` objects to simulate various outcomes. We will have test cases that simulate a successful run (return code 0, valid stdout), a run that produces an error on stderr, a run that returns a non-zero exit code, and a run that "hangs" (by raising a timeout exception). This validates our process handling logic.
    *   **`_parse_output`**: This method will be tested independently by feeding it a variety of strings containing sample QE output text. We will have test cases for standard successful outputs, as well as outputs that might be missing certain optional blocks (like the stress tensor for a molecule). For each input, we will assert that the method extracts the correct floating-point values for energy, forces, and stress, and handles missing data gracefully.

*   **`TrainingEngine`**: The MLIP training library itself (e.g., `mace-torch`) will be completely mocked at the module level. We are testing our integration logic, not the library's ability to train.
    *   **`_compute_delta`**: This method will be tested with a mock baseline calculator whose `get_potential_energy` and `get_forces` methods are configured to return known, fixed values. We will provide a specific `DFTResult` object and assert that the calculated delta values for energy and forces are arithmetically correct.
    *   **`train`**: The primary purpose of this test is to verify the data pipeline into the training library. We will provide a sample dataset and then assert that the `train` method of the mocked library was called exactly once, and that the data passed to it was in the correct format (e.g., a list of specific data objects) and contained the correct delta values calculated in the previous step.

**Integration Testing Approach (Min 300 words):**

While unit tests are essential for verifying components in isolation, a single integration test is required to ensure they communicate and work together correctly in the prescribed sequence. The goal of the Cycle 01 integration test is to verify the "happy path" of the entire workflow, from a starting Atoms object to a saved model file. Crucially, this test will still mock the slow external processes (DFT and model training) to ensure it can run quickly as part of an automated test suite.

The integration test will target the top-level `Orchestrator.run_cycle01_workflow` method. The test setup will be carefully orchestrated:
1.  **Mocking the `LabellingEngine`**: The entire `LabellingEngine` class will be replaced with a mock object. Its `run` method will be configured to return a predefined, successful `DFTResult` object immediately, without actually calling Quantum Espresso. This allows us to control the data that flows from the labelling to the training stage.
2.  **Mocking the `TrainingEngine`**: Similarly, the `TrainingEngine`'s `train` method will be mocked. It will be configured to return a simple placeholder object (e.g., the string "trained_model") without performing any actual training. This allows us to verify that the training step was triggered and to check what happened with its output.
3.  **Using a temporary `AseDB`**: The test will instantiate the `AseDB` to use a temporary database file on the filesystem, which will be automatically cleaned up after the test completes.

The test will then execute the orchestrator's workflow method with a sample `ase.Atoms` object. Finally, a series of assertions will verify the entire chain of events:
*   Assert that the temporary database file was created.
*   Connect to the database and confirm that the initial atom was added and that its state was correctly updated from 'initial' to 'labelled'.
*   Assert that the `TrainingEngine.train` mock was called with the exact `DFTResult` that the mocked `LabellingEngine` was configured to produce.
*   Assert that an attempt was made to save the return value of the mocked `train` method ("trained_model") to a file. The file I/O for saving the model can also be mocked using `mock_open` to avoid filesystem writes.

This single, powerful integration test provides a high degree of confidence that the core data flow, state management, and inter-component communication of the pipeline are all functioning correctly, providing a solid, verified foundation for Cycle 02.
