# Specification: Cycle 1 - Core Engine & Automation Foundation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

Cycle 1 represents the foundational phase of the MLIP-AutoPipe project. The primary objective of this cycle is to construct the essential backbone of the automated workflow, focusing on the most critical components: the DFT labelling and MLIP training engines. This cycle will deliver a system capable of performing the core task of converting a given set of atomic structures into a trained machine learning potential without manual intervention. The central theme is automation. By the end of this cycle, we will have eliminated the most time-consuming manual steps in a typical MLIP development workflow: running DFT calculations and orchestrating the model training process. The successful completion of this cycle will yield a robust, scriptable pipeline that can be executed from a single command.

The scope of this cycle is deliberately focused. We will not be addressing the generation of initial structures or the sophisticated active learning loop. Instead, the cycle assumes that a set of input structures is provided. The system will then process these structures through two main modules. First, `Module C: Labeling Engine`, will take the structures and automatically manage the execution of Quantum Espresso calculations to obtain the necessary energy, force, and stress labels. This involves programmatically generating input files with physically sound parameters, managing the execution of the DFT code, parsing the output, and storing the results in a structured database. Second, `Module D: Training Engine`, will consume the labelled data from the database and conduct the training of a delta-learning MLIP. This module will handle all aspects of the training process, from data preparation to the final serialization of the trained model.

A crucial deliverable of this cycle is the implementation of the "Two-Tier Configuration" strategy. This innovative approach is key to achieving our goal of a user-friendly system. We will develop the `ConfigExpander`, a heuristic engine that takes a minimal user-provided configuration file and intelligently populates it with all the detailed parameters required for the DFT and training steps. This will abstract away the complexity of the underlying tools, allowing users to interact with the system at a high level. By the conclusion of Cycle 1, we will have a functional, albeit not yet fully autonomous, pipeline that demonstrates the viability of our automated approach and provides a solid, extensible foundation upon which the subsequent cycles will build. The successful implementation of this cycle will result in a command-line tool that can be pointed at a set of structures and a minimal configuration file, and will output a trained potential without any further user interaction. This will represent a major step forward in the automation of MLIP development and will provide a solid foundation for the more advanced features to be added in later cycles. The importance of this cycle cannot be overstated, as it will serve as the bedrock for all future development. A robust and reliable core engine is essential for the success of the project as a whole.

## 2. System Architecture

The system architecture for Cycle 1 is a linear, two-stage pipeline consisting of the Labeling Engine (Module C) and the Training Engine (Module D), orchestrated by a central workflow manager that is configured via the Two-Tier configuration system. This architecture is a minimalist but essential slice of the final system, designed to prove the core concept of automated labelling and training.

**Workflow Orchestration and Configuration:**
The process begins with the user providing a minimal `input.yaml`. This file is the sole point of interaction for the user in this cycle. The main script, which acts as the workflow orchestrator, first invokes the `ConfigExpander`. This component reads the user's input, identifies the material system, and applies a set of built-in rules and heuristics to generate a fully specified configuration. For example, it will select the appropriate SSSP pseudopotentials and plane-wave cutoffs for Quantum Espresso based on the elements present. It will also define the parameters for the delta-learning model in the training engine. This fully expanded configuration, saved as `exec_config_dump.yaml`, becomes the immutable set of instructions for the rest of the workflow, ensuring reproducibility. The ConfigExpander will be designed to be highly modular, with a clear separation between the rules for different DFT codes and MLIP models. This will allow for easy extension to support other codes in the future.

**Data Flow:**
The data flow is strictly sequential. The orchestrator takes a user-provided file containing the initial atomic structures (e.g., an ASE trajectory file). It then passes these structures to `Module C: Labeling Engine`.

*   **Module C (Labeling Engine):** This module is architected as a robust wrapper around the Quantum Espresso command-line interface. It is designed to be completely non-interactive. For each atomic structure it receives, it performs the following steps:
    1.  **Input Generation:** It uses the parameters from `exec_config_dump.yaml` to dynamically generate a QE input file (`.pwi`). This includes setting the correct cell parameters, atomic positions, k-point mesh density, and SCF convergence parameters. Logic will be included to handle complexities like magnetism, automatically enabling spin-polarised calculations for magnetic elements.
    2.  **Execution Management:** It spawns a subprocess to execute `pw.x`. The command, including MPI parallelisation, will be configurable. The module will monitor the subprocess, capturing standard output and error streams.
    3.  **Output Parsing and Error Handling:** Upon completion, the module parses the QE output file to extract the final energy, atomic forces, and system stress tensor. It will be designed to detect common failure modes, such as SCF non-convergence. An initial, simple error recovery mechanism will be implemented, such as attempting a re-calculation with more conservative mixing parameters.
    4.  **Data Persistence:** All extracted data, along with the original structure and key calculation parameters (metadata), are written to a central ASE database. This creates a traceable record of every calculation performed.

*   **Module D (Training Engine):** Once the labelling process is complete, the orchestrator invokes `Module D`. This module is responsible for the machine learning aspect of the pipeline.
    1.  **Data Loading:** It queries the ASE database to retrieve the full set of labelled structures.
    2.  **Data Preparation:** It prepares the data for training. This involves splitting the data into training and validation sets and configuring the delta-learning framework. This means it will calculate the energy and forces from a baseline reference potential (e.g., ZBL) and compute the residuals that the MLIP needs to learn.
    3.  **Model Training:** It interfaces with the chosen MLIP library (e.g., the Atomic Cluster Expansion - ACE - framework). It passes the prepared data to the library's training routine and manages the training process.
    4.  **Model Serialization:** After training is complete, the module saves the trained potential to a file. This file can then be used by a simulation engine in later cycles.

This two-module architecture, governed by the configuration system, forms a solid, testable, and extensible foundation for the project. The clear separation of concerns between the two modules will allow for independent development and testing, and the use of a central database will ensure that the data is managed in a consistent and reproducible manner. The modular design will also facilitate future enhancements, such as the addition of support for other DFT codes or MLIP models.

## 3. Design Architecture

The design architecture for Cycle 1 will be implemented in Python, adhering to principles of modularity and separation of concerns. The project structure will be established using `pyproject.toml` and `uv` for dependency management.

**Project Structure:**
```
src/mlip_autoflow/
├── __init__.py
├── main.py
├── config/
│   ├── expander.py
│   └── models.py      # Pydantic models for configuration
├── data/
│   └── ase_db_manager.py
└── modules/
    ├── c_labeling_engine.py
    └── d_training_engine.py
```

**Class and Method Definitions:**

*   **`config.models.py`**: This will contain Pydantic models. `UserConfig` will define the schema for `input.yaml`. `DFTConfig`, `TrainingConfig`, and a top-level `FullConfig` will define the schema for the expanded `exec_config_dump.yaml`. This provides strong typing and validation.

*   **`config.expander.py`**: A `ConfigExpander` class will be implemented.
    *   `__init__(self, user_config: UserConfig)`: Takes the parsed user config.
    *   `expand(self) -> FullConfig`: Contains the heuristic logic. It will include methods like `_get_dft_params(elements)` and `_get_training_defaults()`. It will return a fully populated `FullConfig` object.

*   **`data.ase_db_manager.py`**: A `AseDBManager` class will handle all database interactions.
    *   `__init__(self, db_path: str)`: Connects to the ASE database file.
    *   `write_calculation(self, atoms: ase.Atoms, results: dict, metadata: dict)`: Writes a completed DFT calculation to the database.
    *   `get_all_labeled_structures(self) -> List[ase.Atoms]`: Queries the database and returns all structures that have successfully been labeled.

*   **`modules.c_labeling_engine.py`**: A `LabelingEngine` class will manage DFT calculations.
    *   `__init__(self, config: DFTConfig, db_manager: AseDBManager)`: Initializes with the DFT configuration and a database manager instance.
    *   `run_labeling_workflow(self, structures: List[ase.Atoms])`: The main public method. It iterates through the list of `ase.Atoms` objects and calls a private method for each.
    *   `_run_single_calculation(self, atoms: ase.Atoms) -> None`: A private method that generates the QE input, executes `pw.x` via `subprocess.run`, parses the output, and uses the `db_manager` to save the results.

*   **`modules.d_training_engine.py`**: A `TrainingEngine` class will manage MLIP training.
    *   `__init__(self, config: TrainingConfig, db_manager: AseDBManager)`: Initializes with the training configuration and database manager.
    *   `run_training_workflow(self)`: The main public method. It retrieves data from the `db_manager`, prepares it, trains the model, and saves it.
    *   `_prepare_delta_dataset(self, atoms_list: List[ase.Atoms])`: A private method to calculate the residuals for delta learning.
    *   `_train_model(self, dataset)`: A private method that interfaces with the chosen MLIP library (e.g., `pyace`).

*   **`main.py`**: This will be the main CLI entry point, using a library like `Typer`. It will orchestrate the entire process: parse the user config, instantiate and run the `ConfigExpander`, instantiate the `AseDBManager`, `LabelingEngine`, and `TrainingEngine`, and then call their main workflow methods in sequence.

This object-oriented design ensures that each class has a single responsibility, making the system easier to understand, test, and maintain. The use of dependency injection (e.g., passing the `db_manager` to the engines) decouples the modules from the specific data storage implementation.

## 4. Implementation Approach

The implementation of Cycle 1 will be a systematic process, starting with the project setup and progressively building up the functionality of each module. The development will be guided by a test-driven approach where possible.

**Step 1: Project Scaffolding and Dependency Management**
The first step is to set up the project structure as defined in the Design Architecture. We will initialize a Git repository and create the `src/mlip_autoflow` directory structure. A `pyproject.toml` file will be created to define project metadata and dependencies. Key dependencies to be included are:
*   `uv`: For environment and package management.
*   `typer`: For the command-line interface.
*   `pydantic`: For configuration and data modeling.
*   `ase`: The core library for atomic structures and database management.
*   `numpy`: For numerical operations.
*   `pyace` (or a similar ACE/MACE framework): The MLIP training library.
*   `pytest`: For testing.

**Step 2: Configuration System (`config` module)**
We will begin by implementing the Pydantic models in `config/models.py`. This will define the clear data contracts that the rest of the system will use. Following this, we will develop the `ConfigExpander` class. The initial version of the heuristic engine will be simple, perhaps supporting only a few elements and providing sensible default parameters for DFT and training. The logic for selecting pseudopotentials and cutoffs based on SSSP recommendations will be a key part of this implementation. We will write unit tests to ensure that a minimal `input.yaml` is correctly expanded into a valid `FullConfig` object.

**Step 3: Data Persistence (`data` module)**
Next, we will implement the `AseDBManager` class. This class will encapsulate all interactions with the `ase.db` SQLite database. We will implement the methods for writing new calculation results and for querying the database to retrieve all labeled data. The interface will be simple and focused on the needs of the labelling and training engines.

**Step 4: DFT Automation (`modules.c_labeling_engine`)**
This is a critical module. We will develop the `LabelingEngine` class. The implementation will start with the input file generation. We will create a template for the Quantum Espresso input file and use string formatting or a templating engine to populate it with parameters from the `DFTConfig` object. Next, we will implement the execution logic using Python's `subprocess` module. Careful attention will be paid to capturing the output and handling potential errors. The output parsing logic will be implemented to robustly extract the required data (energy, forces, stress) from the QE output file. We will integrate this module with the `AseDBManager` to persist the results. The entire workflow will be wrapped in the `run_labeling_workflow` method.

**Step 5: MLIP Training (`modules.d_training_engine`)**
With the labelling engine in place, we will implement the `TrainingEngine`. This class will first use the `AseDBManager` to load the data. We will then implement the logic for delta learning. This involves selecting a baseline potential (e.g., ZBL) and calculating the difference between its predictions and the DFT labels. The core of this module will be the interface to the MLIP library (e.g., `pyace`). We will write a wrapper function that takes our prepared dataset and passes it to the training function of the library. Finally, we will implement the logic to save the trained model to a file.

**Step 6: CLI Orchestration (`main.py`)**
The final step is to tie everything together in the main CLI script. We will use `Typer` to create a simple command, e.g., `mlip-pipe run-cycle1 <input.yaml> <structures.traj>`. This script will be responsible for the high-level orchestration:
1.  Load and parse `input.yaml`.
2.  Instantiate and run the `ConfigExpander`.
3.  Load the initial structures from the provided file.
4.  Instantiate the `AseDBManager`.
5.  Instantiate the `LabelingEngine` and run its workflow.
6.  Instantiate the `TrainingEngine` and run its workflow.
7.  Print a message indicating successful completion.

Throughout this process, we will add logging to provide visibility into the pipeline's execution and aid in debugging. This systematic approach will ensure that each component is well-tested and that the final integrated system is robust and reliable.

## 5. Test Strategy

The test strategy for Cycle 1 is designed to ensure the correctness and robustness of the core pipeline components. We will use a combination of unit and integration testing, leveraging the `pytest` framework.

**Unit Testing Approach:**
Unit tests will focus on testing individual classes and their methods in isolation, using mocking to remove external dependencies.
*   **`ConfigExpander`**: We will create several example `input.yaml` files (e.g., for a simple metal, a magnetic system) and write tests to assert that the `expand()` method produces a `FullConfig` object with the expected parameters. For instance, we will check that the correct pseudopotentials are selected and that spin-polarised DFT settings are enabled for magnetic elements. This ensures the heuristic engine is behaving as expected.
*   **`AseDBManager`**: We will test the database manager against a temporary, in-memory SQLite database. Tests will cover writing a new entry, reading it back, and querying for all entries. We will verify that the data is stored and retrieved correctly without corruption.
*   **`LabelingEngine`**: This requires mocking the execution of Quantum Espresso. We will replace the `subprocess.run` call with a mock object. We will have several tests:
    1.  A "happy path" test where we provide a mock successful QE output file and verify that the engine parses it correctly and writes the correct data to the (mocked) database manager.
    2.  An "error path" test where we provide a mock QE output file indicating an SCF convergence error. We will test that the engine correctly identifies the error and, in the future, that it triggers the error recovery logic.
    3.  A test to verify that the generated QE input file string is correct for a given set of parameters.
*   **`TrainingEngine`**: We will create a small, fixed dataset (a dummy `ase.db` file) for testing. We will mock the underlying MLIP training library (`pyace`). The tests will ensure that the `TrainingEngine` correctly loads the data, prepares the delta-learning dataset by subtracting the baseline potential, and calls the training library with the correctly formatted data.

**Integration Testing Approach:**
Integration tests will verify the interaction and data flow between the modules. For Cycle 1, the primary integration test will be a "mini-pipeline" that connects the `LabelingEngine` and `TrainingEngine`.
*   **Test Scenario:** The test will orchestrate the following workflow:
    1.  **Setup:** Programmatically create a small set of `ase.Atoms` objects (e.g., 2-3 simple structures). Create a valid `FullConfig` object.
    2.  **Execution:**
        *   Instantiate a real `AseDBManager` with a temporary database file.
        *   Instantiate the `LabelingEngine`. Instead of calling the real QE, we will configure it to use a "dummy" QE wrapper script that we provide. This script will take the generated input file and produce a valid, deterministic output file (pre-prepared for the test).
        *   Run the `run_labeling_workflow`.
        *   Instantiate the `TrainingEngine` using the same temporary database.
        *   Run the `run_training_workflow`. This time, we will use the *real* MLIP library but configure it for a very short training run (e.g., 1-2 epochs) on the small dataset generated by the first step.
    3.  **Verification:**
        *   Assert that the temporary database contains the expected number of labeled structures.
        *   Assert that the `TrainingEngine` successfully produces a model file.
        *   Optionally, load the model and perform a simple prediction to ensure it is a valid potential file.

This integration test will validate that the data contract between the modules (the database schema) is correct and that the components function correctly when connected. This provides high confidence in the end-to-end workflow of the Cycle 1 deliverables. The combination of unit and integration tests will ensure that the core engine of the MLIP-AutoPipe system is robust, reliable, and produces scientifically valid results.
