# Cycle 1 Specification: Core Engine and Automation Foundation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

This document provides the detailed technical specification for the inaugural development cycle of the MLIP-AutoPipe project. The singular focus of Cycle 1 is to forge the foundational backbone of the entire system—the core engine responsible for processing atomic structures into a trained Machine Learning Interatomic Potential (MLIP). This cycle is deliberately scoped to exclude the complexities of automated structure generation and active learning, concentrating instead on establishing a robust, reliable, and automated pipeline for the two most critical transformations: from structure to DFT data, and from DFT data to a functional MLIP. This involves the implementation of two key modules: **Module C (Labelling Engine)** and **Module D (Training Engine)**. By the end of this cycle, we will have a system that, while requiring manual input of structures, can autonomously manage the entire labelling and training workflow, thereby proving the viability of the project's core components and providing a stable, well-tested platform upon which all subsequent functionality will be built.

The successful completion of this cycle will yield a system with a clearly defined capability: a user will be able to provide a simple directory of atomic structure files (e.g., in CIF or POSCAR format) and, by running a single command, receive a trained MLIP model. The intelligence of the system will be demonstrated through its `LabellingEngine`, which will completely automate the traditionally manual and error-prone process of setting up and running Density Functional Theory (DFT) calculations. It will interface with a user-provided Quantum Espresso installation, but it will take full responsibility for generating the input files. This includes embedding expert-level logic for selecting calculation parameters according to the rigorous, community-vetted SSSP (Standard Solid State Pseudopotentials) protocol, ensuring high-quality, reproducible DFT data without requiring the user to have deep expertise in the intricacies of the simulation package. Furthermore, the engine will be designed for robustness, incorporating automated error detection and recovery mechanisms to handle common issues like Self-Consistent Field (SCF) convergence failures. All successfully computed data will be meticulously stored in a versioned, queryable database. This data will then be consumed by the `TrainingEngine`, which will implement the physically-motivated "delta learning" strategy. This approach, where the model learns to predict the residual between a simple physical potential and the true DFT forces, is critical for ensuring the final MLIP is stable and behaves correctly in regions far from the training data. The tangible output of this cycle is a non-trivial proof of concept: a functional, end-to-end data processing pipeline that serves as the essential core for the more advanced automation features planned in subsequent cycles.

## 2. System Architecture

The architecture for Cycle 1 is intentionally simple, linear, and non-iterative. It is designed to perfect the core data transformation workflow without the added complexities of the feedback loops or automated data generation that will be introduced later. The system at this stage operates as a batch processing pipeline, taking a collection of static input files and producing a single output artifact, the trained MLIP. The workflow is manually initiated and processes a fixed, user-provided dataset. At its heart, the architecture connects two primary components, Module C and Module D, via a durable data persistence layer, the ASE database. This separation is a critical design choice. By decoupling the DFT calculation phase from the machine learning phase, we gain significant flexibility. We can rerun the training on a completed dataset with different MLIP parameters without needing to reperform the expensive DFT calculations. It also makes the system more robust; a failure in the training module will not affect the integrity of the precious, computationally expensive data stored in the database.

```mermaid
graph TD
    A[User Input: Directory of Structures] --> B[Module C: Labelling Engine];
    B -- Creates & Executes --> C{Quantum Espresso Input};
    C -- Monitored Subprocess --> B;
    B -- Parses Output --> D[Parsed DFT Results (Energy, Forces)];
    D -- Stores Data & Metadata --> E[ASE Database];
    E -- Provides Complete Dataset --> F[Module D: Training Engine];
    F -- Applies Delta Learning --> G[MLIP Model Training];
    G -- Saves Final Artifact --> H[Trained MLIP Model];

    subgraph "Core Pipeline - Cycle 1"
        direction LR
        B;
        F;
    end
```

**Detailed Workflow Description:**

1.  **Manual Input:** The process is initiated by a user via a command-line call. The user must provide the path to a directory containing one or more atomic structure files in a standard format (e.g., POSCAR, CIF, XYZ) that can be parsed by the ASE (Atomic Simulation Environment) library.

2.  **Labelling (Module C):** The `LabellingEngine` is the workhorse of this cycle. It orchestrates the entire DFT calculation process.
    a.  **Structure Ingestion:** It begins by iterating through all provided structure files, parsing each one into an in-memory ASE `Atoms` object.
    b.  **Parameter Automation:** For each `Atoms` object, it determines the set of unique chemical elements. It then queries an internal database representing the SSSP protocol to retrieve the recommended pseudopotential files, plane-wave energy cutoffs (`ecutwfc`), and charge density cutoffs (`ecutrho`). This automated selection is a cornerstone of the system's "expert-in-the-box" design.
    c.  **Input Generation:** With these parameters, it dynamically generates a complete and valid input file for Quantum Espresso's `pw.x` executable.
    d.  **Execution and Monitoring:** The engine then executes `pw.x` as a managed subprocess, capturing its standard output and error streams in real-time. This allows it to monitor the progress of the DFT calculation.
    e.  **Parsing and Error Handling:** Upon completion of the subprocess, the engine parses the output file. Using robust regular expressions, it searches for keywords indicating success or failure. If successful, it extracts the final total energy, the forces on each atom, and the system's stress tensor. If it detects a specific, recoverable error, such as "SCF not converged," it will trigger its error recovery logic (e.g., generating a new input file with more conservative electronic mixing parameters) and retry the calculation a configurable number of times.
    f.  **Data Persistence:** Once a structure is successfully calculated, the engine stores the results. The original `Atoms` object is updated with the calculated energy, forces, and stress, and this entire bundle of information, along with critical metadata (like the QE version and input parameters used), is written to the central ASE database.

3.  **Training (Module D):** After the labelling for all input structures is complete, the `TrainingEngine` takes over.
    a.  **Data Aggregation:** It first connects to the ASE database and retrieves the entire labelled dataset that was just created.
    b.  **Delta Calculation:** It then applies the delta learning strategy. For each structure, it calculates the forces and energy using a simple, analytical reference potential (e.g., ZBL). It then subtracts these values from the high-accuracy DFT results to get the "residual" or "delta" values.
    c.  **Model Fitting:** The MLIP model (e.g., an Atomic Cluster Expansion model) is then trained exclusively on these delta values. This is the core machine learning step, where the model's parameters are optimized to best predict the quantum mechanical correction to the simple physical baseline.
    d.  **Artifact Serialization:** Once the training is complete, the engine saves the final, trained model object to a file. This model artifact contains all the necessary information to be used as a calculator in a future simulation.

## 3. Design Architecture

This cycle lays the architectural groundwork for the entire project, establishing the core class structures, directory layout, and data models that will be extended in subsequent cycles. The design emphasizes modularity, clear interfaces between components, and strict data validation to ensure robustness.

**Key Classes and Modules:**

-   **`src/mlip_autopipe/modules/labelling_engine.py`:**
    -   **`LabellingEngine` class:** This class encapsulates all logic related to DFT calculations.
        -   `__init__(self, config)`: Its constructor will take a configuration object that specifies crucial parameters like the path to the `pw.x` executable and the maximum number of retries for failed calculations.
        -   `run(self, structures: List[Atoms])`: The main public method that orchestrates the processing of a list of input structures.
        -   `_prepare_qe_input(self, atoms: Atoms) -> str`: A private method responsible for the complex logic of generating the text content of a Quantum Espresso input file from an `Atoms` object and the automated SSSP parameters. This method will be heavily unit-tested to ensure correctness.
        -   `_select_sssp_parameters(self, elements: List[str]) -> dict`: A critical method that contains the "expert knowledge." It will read from a bundled data file (e.g., a JSON file) containing the SSSP library information to determine the correct pseudopotentials and cutoffs for any given combination of elements.
        -   `_execute_dft(self, input_content: str, working_dir: Path) -> subprocess.CompletedProcess`: A wrapper around Python's `subprocess.run`. It will be responsible for creating a temporary directory for the calculation, writing the input file, launching `pw.x`, and returning the completed process object.
        -   `_parse_output(self, output_content: str) -> dict`: A method dedicated to parsing the output from `pw.x`. It will use a series of regular expressions to locate and extract the required numerical data. It will be designed to be robust to minor formatting changes between QE versions.
        -   `_handle_error(self, run_context) -> bool`: This method will contain the logic for recoverable errors. It will inspect the output of a failed run and decide whether to attempt a retry with modified parameters.

-   **`src/mlip_autopipe/modules/training_engine.py`:**
    -   **`TrainingEngine` class:** This class handles the machine learning aspect of the pipeline.
        -   `__init__(self, config)`: Takes a configuration object specifying the type of MLIP to train (e.g., 'ACE'), the choice of baseline potential, and other hyperparameters.
        -   `run(self, dataset: List[Atoms])`: The main public method that takes the full labelled dataset and returns a trained model object.
        -   `_prepare_delta_dataset(self, dataset)`: A private method that implements the core delta learning logic, subtracting the baseline potential's predictions from the DFT labels.
        -   `_train_model(self, delta_dataset)`: This method will act as a bridge to an external MLIP training library (like `pyace` or a similar framework). It will be responsible for converting the data into the format expected by the library and invoking its fitting function.
        -   `_save_model(self, model, filepath: Path)`: A utility method to serialize the trained model artifact to disk.

-   **`src/mlip_autopipe/common/db_interface.py`:**
    -   **`DatabaseInterface` class:** A dedicated data access object (DAO) to abstract away the details of the database.
        -   `__init__(self, db_path)`: Connects to the SQLite file that serves as the ASE database.
        -   `write_atoms(self, atoms: Atoms, metadata: dict)`: A method to write a single labelled `Atoms` object, along with its associated metadata dictionary, to the database. This ensures that provenance is always stored with the data.
        -   `get_all_atoms(self) -> List[Atoms]`: A method to query and return the entire dataset from the database.

**Data Models:** Strict data validation is enforced using Pydantic. In the `config/models.py` file, we will define Pydantic models for the relevant sections of the configuration file. This ensures that the pipeline fails fast with a clear error message if the configuration is invalid, rather than crashing midway through a long calculation.

## 4. Implementation Approach

The implementation of Cycle 1 will be methodical, starting with the foundational elements and building outwards, with a strong emphasis on unit testing at each stage before proceeding to integration.

1.  **Project Scaffolding:** The very first step is to establish the project's structure and development environment. This involves creating the `pyproject.toml` file, which serves as the definitive source for project metadata and dependencies. Initial dependencies will be declared, including `ase` (for atomic structure manipulation), `numpy` (for numerical data), `pydantic` (for configuration validation), and the chosen MLIP training framework (e.g., `pyace`). The directory structure (`src/mlip_autopipe`, `tests/`, etc.) will be created as outlined in the design architecture. This initial setup ensures a clean and standardized development process.

2.  **Database Interface:** The `DatabaseInterface` class will be implemented early on. This component is simple but absolutely crucial as it defines the contract for how data is stored and retrieved. It will be built as a thin wrapper around the `ase.db` module. We will write the first unit tests for the project to verify the functionality of this class, ensuring that we can write an `Atoms` object with attached results and metadata, and then read it back perfectly.

3.  **Labelling Engine (Module C) Implementation:** This is the most complex part of Cycle 1.
    a.  We will begin by creating the data file (e.g., `sssp_protocol.json`) that encodes the SSSP library information. The `_select_sssp_parameters` method will be implemented to read from this file.
    b.  Next, the `_prepare_qe_input` method will be developed. Its implementation will be guided by extensive unit tests that assert the correctness of the generated input file string for various atomic systems.
    c.  The core of the module, the `_execute_dft` and `_parse_output` methods, will be developed next. The parsing logic will be built and tested against a pre-saved collection of sample Quantum Espresso output files. This is a critical step, as it allows us to develop and test the parser thoroughly without the overhead of running actual DFT calculations in our CI environment.
    d.  The error handling logic (`_handle_error`) will be implemented and tested by feeding the parser with sample failed-run outputs.
    e.  Finally, the main `run` method will be implemented to tie all the private helper methods together into a coherent workflow.

4.  **Training Engine (Module D) Implementation:**
    a.  The `_prepare_delta_dataset` method will be implemented first. This will involve choosing a simple baseline potential available within a library like ASE (e.g., `LennardJones`) and applying it to a sample dataset.
    b.  The `_train_model` method will be written as a wrapper. Its main responsibility is data transformation—ensuring the dataset is in the precise format required by the external MLIP training library. The initial implementation will focus on correctness and will not involve complex hyperparameter tuning.
    c.  The `_save_model` method will be implemented to ensure the trained potential is saved in a standard, reusable format.

5.  **Integration and CLI Stub:** To tie everything together, a minimal command-line entry point will be created in `cli.py`. Initially, this can be a simple script using Python's built-in `argparse` that takes the required inputs (path to structures, path to config). The main logic will instantiate the `LabellingEngine`, run it, and then pass its results (via the database) to an instance of the `TrainingEngine`. The culmination of the cycle will be an end-to-end integration test on a simple, known system (e.g., a few structures of a Silicon unit cell) to verify the entire chain from input files to a trained model file.

## 5. Test Strategy

The test strategy for Cycle 1 is designed to build a foundation of confidence in the core data processing capabilities of the pipeline. Given that this cycle involves interfacing with an external, computationally expensive program (Quantum Espresso), the strategy is heavily weighted towards comprehensive unit tests that use mock objects and pre-saved data to isolate components and test them thoroughly without relying on the external dependency.

**Unit Testing Approach:**
The goal of the unit tests is to verify the correctness of each piece of logic in isolation.
-   **Labelling Engine:** The `LabellingEngine` is the most critical component to test at this stage. We will create an extensive test suite that covers every aspect of its operation. The `_select_sssp_parameters` method will be tested with various combinations of chemical elements to ensure it always returns the correct, expected pseudopotential names and energy cutoff values from our SSSP data file. The `_prepare_qe_input` method will be tested by feeding it different `ASE.Atoms` objects (representing bulk crystals, molecules, magnetic systems, etc.) and asserting that the generated Quantum Espresso input string is syntactically correct and contains all the necessary parameters. The most important unit tests will target the `_parse_output` method. We will create and maintain a dedicated directory within our test suite containing a collection of real QE output files. This collection will include examples of successful runs, runs that failed due to SCF non-convergence, runs that hit a walltime limit, and runs that crashed with an error. The unit tests will feed the contents of these files to the parser and assert that it either correctly extracts the numerical data (energy, forces, stress) or correctly identifies the specific type of error that occurred. The error handling logic (`_handle_error`) will be tested by simulating a failed run and asserting that the engine correctly modifies the input parameters (e.g., lowers the mixing beta) and attempts to retry the calculation.

-   **Training Engine:** Unit tests for the `TrainingEngine` will be designed to verify its data transformation and model interaction logic. We will create a small, synthetic dataset of `Atoms` objects with known, pre-calculated energies and forces. We will test the `_prepare_delta_dataset` method by providing a known baseline potential (like a simple Lennard-Jones potential) and asserting that the resulting delta values are mathematically correct down to a tight numerical tolerance. The interaction with the external MLIP training library will be mocked using `pytest-mock`. This allows us to test that our `TrainingEngine` correctly formats and passes the delta dataset to the library's fitting function, and correctly handles the model object that is returned, without needing to perform a real (and slow) model training in the unit test.

**Integration Testing Approach:**
While unit tests verify the components in isolation, integration tests are essential to ensure they work together correctly. For Cycle 1, the focus is on a single, comprehensive "mini-pipeline" test.
-   **Mini-Pipeline Test:** This test will be designed to be as lightweight as possible while still verifying the entire data flow through the Cycle 1 system. It will be run as part of the CI pipeline, using a real Quantum Espresso installation. The test will be configured to run on a computationally trivial system, such as a single hydrogen atom or a 2-atom silicon cell with a very low energy cutoff and a minimal k-point grid, ensuring it completes in seconds. The test will perform the following actions and assertions:
    1.  It will start with one or two structure files in a temporary directory.
    2.  It will invoke the main command-line entry point to run the pipeline on this data.
    3.  The test will assert that a real `pw.x` process is launched, confirming the subprocess execution logic is working.
    4.  It will check that the results are correctly written to a temporary ASE database file.
    5.  It will then verify that the `TrainingEngine` is able to successfully read this data from the database.
    6.  Finally, it will assert that a valid, non-empty MLIP model file is created in the output directory.
    This end-to-end test provides a crucial validation of the entire system implemented in Cycle 1. It confirms that file I/O, subprocess management, DFT parsing, database interaction, and model serialization are all correctly integrated and functioning as a cohesive whole.
