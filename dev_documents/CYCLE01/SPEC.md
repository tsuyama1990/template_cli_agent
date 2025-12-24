# CYCLE01 Specification: Core Engine and Workflow Foundation

## 1. Summary

This inaugural development cycle, CYCLE01, focuses on establishing the foundational infrastructure of the MLIP-AutoPipe system. The primary objective is to build a robust, end-to-end workflow that serves as the backbone for all subsequent features. This cycle will deliver a system capable of taking a user-provided set of atomic structures, performing automated DFT calculations to label them with forces, stresses, and energies, and subsequently training a functional Machine Learning Interatomic Potential (MLIP). This establishes the core data pipeline: from raw structure to trained model.

The key components to be implemented are the Labeling Engine (Module C) and the Training Engine (Module D). Module C will interface with Quantum Espresso (QE), automatically generating the necessary input files and parsing the outputs. It will incorporate basic robustness checks and parameter handling to ensure reliable calculations. Module D will take this labelled data and use it to train an MLIP based on the Atomic Cluster Expansion (ACE) framework.

A critical deliverable for this cycle is the setup of a modern, professional Python project environment. This includes the creation of a `pyproject.toml` file to manage all project dependencies and metadata, and the adoption of `uv` as the package manager to ensure fast, deterministic, and reproducible environments. We will also design and implement the initial schema for the project's central database (using the ASE DB library), which will be crucial for tracking data provenance and ensuring the traceability of every calculation. Furthermore, the foundational design for the Two-Tier Configuration strategy will be implemented, defining the Pydantic models for both `input.yaml` and `exec_config_dump.yaml`, although the automated heuristic expansion will be minimal in this first cycle.

## 2. System Architecture

The architecture for CYCLE01 is a simplified, linear version of the final system, focusing exclusively on the Labeling-to-Training pipeline. At this stage, the workflow is not yet fully autonomous; it relies on the user to provide the initial set of atomic configurations. The `Config Expander` will exist but will perform minimal logic, primarily validating and passing through user-specified parameters from a more detailed `input.yaml` to the `exec_config_dump.yaml`.

The data flow is as follows:
1.  **User Input**: The user provides a directory of initial structures (e.g., in CIF or POSCAR format) and a configuration file (`input.yaml`) specifying the elements, DFT parameters, and MLIP training settings.
2.  **Configuration**: The system parses `input.yaml` using Pydantic models and generates the `exec_config_dump.yaml`.
3.  **DFT Labeling (Module C)**: The `LabelingEngine` iterates through the provided structures. For each structure, it generates a QE input file (`.in`) based on the parameters in `exec_config_dump.yaml`. It then executes the `pw.x` binary (the path and command for which are user-provided). After the calculation completes, it parses the QE output file (`.out`) to extract the final energy, atomic forces, and system stress.
4.  **Database Storage**: The original structure and the extracted DFT results are stored together in the ASE database. Each entry is timestamped and linked to the configuration used, establishing a clear record of provenance.
5.  **Data Aggregation**: The `TrainingEngine` queries the database to retrieve all successfully labelled structures and their corresponding DFT data.
6.  **MLIP Training (Module D)**: The aggregated data is fed into the training workflow. This cycle will use the ACE framework. The engine will fit an ACE model to the data, performing a "Delta Learning" by default if a baseline potential is specified. The resulting trained potential file (e.g., a `.yace` file) is saved to disk, and its path and metadata are recorded in the database.

The system will be structured as a command-line application. The user will initiate the workflow with a command like `uv run mlip-pipe run --structures-dir /path/to/structures input.yaml`.

## 3. Design Architecture

This cycle will lay the groundwork for the entire project's codebase, located within the `src/mlip_autoprope` directory.

**Key Files and Classes:**

*   `src/mlip_autoprope/cli.py`: This file will house the main entry point for the command-line interface. It will use a library like `Typer` to define commands, arguments, and options. The initial command `run` will be implemented here.
*   `src/mlip_autoprope/config/core.py`: Pydantic models `InputConfig` and `ExecConfig` will be defined here. These models will provide strict data validation for the two-tier configuration files. `ExecConfig` will contain detailed subsections for `dft_compute` and `mlip_training`.
*   `src/mlip_autoprope/database/client.py`: A `DatabaseClient` class will be created to abstract all interactions with the ASE database. It will provide methods like `write_calculation(atoms, results)` and `get_all_calculations()`. This isolates the database logic from the main workflow.
*   `src/mlip_autoprope/modules/c_labeling_engine.py`:
    *   `LabelingEngine` class: The main orchestrator for this module.
    *   `_generate_qe_input(atoms, params)`: A private method to create the text for a QE input file.
    *   `_parse_qe_output(outfile_path)`: A private method to parse the output and extract key numerical data.
    *   `run(structures, config)`: The main public method that takes a list of ASE `Atoms` objects and the execution config, runs QE on each, and returns the results. It will use Python's `subprocess` module to call `pw.x`.
*   `src/mlip_autoprope/modules/d_training_engine.py`:
    *   `TrainingEngine` class: The orchestrator for the training process.
    *   `run(data, config)`: The main public method. It will take the labelled data, configure the ACE training parameters based on the config, execute the training process (leveraging the appropriate Python library for ACE), and save the final model. It will implement the Delta Learning logic, subtracting the baseline potential's contribution from the DFT labels before training.

The project structure will be initialized with `uv`, and `pyproject.toml` will be configured with initial dependencies, including `ase`, `typer`, `pydantic`, and the chosen ACE training library.

## 4. Implementation Approach

The implementation will proceed in a logical, step-by-step manner.

1.  **Project Scaffolding**: Initialize the project directory structure as outlined in the System Architecture. Create the `pyproject.toml` file and use `uv` to create and manage the virtual environment. Add initial dependencies.
2.  **Configuration Models**: Implement the Pydantic models in `src/mlip_autoprope/config/core.py`. Start with simple, flat structures and add nesting (`DFTConfig`, `TrainingConfig`) as needed. Ensure strict validation is enabled.
3.  **Database Interface**: Implement the `DatabaseClient` class. It will be initialized with a path to the database file. The initial implementation will focus on creating the necessary tables (if they don't exist) and writing/reading calculation data. Use the `ase.db` API for this.
4.  **Labeling Engine (Module C)**: This is the most complex part of this cycle.
    *   Begin with the `_generate_qe_input` function. It should take an ASE `Atoms` object and a configuration dictionary and return a formatted string. This can be tested independently.
    *   Implement the `_parse_qe_output` function. This will involve careful string matching and regular expressions to reliably find and convert the energy, force, and stress data from the QE output. This should also be tested independently with sample QE output files.
    *   Implement the main `run` method. It will orchestrate the process of calling the generator, running `subprocess.run` to execute `pw.x`, and then calling the parser. It must include error handling for cases where the QE calculation fails (e.g., SCF does not converge).
5.  **Training Engine (Module D)**:
    *   Implement the `TrainingEngine` class. The `run` method will first call the `DatabaseClient` to fetch all necessary data.
    *   It will then prepare this data for the ACE training library. This may involve converting units or restructuring the data.
    *   Implement the Delta Learning logic: if `delta_learning` is enabled, calculate the energy/forces from a simple reference potential (e.g., using ASE's built-in calculators) and subtract them from the DFT target values.
    *   Finally, it will call the ACE library's main fitting function and save the output model.
6.  **CLI and Integration**: Implement the `cli.py` file. The `run` command will tie everything together: it will parse the config, load the structures, instantiate and run the `LabelingEngine`, instantiate and run the `TrainingEngine`, and print a success message.

This incremental approach, with independent testing of the parsing and generation logic, will ensure a robust and maintainable foundation for the project.

## 5. Test Strategy

The testing for CYCLE01 is foundational, ensuring that the core components of the data processing pipeline are reliable.

### Unit Testing Approach

The focus of unit testing will be on isolating each component and verifying its logic with controlled inputs.
*   **Configuration**: Test the Pydantic models by attempting to parse valid and invalid `input.yaml` files. Assert that validation errors are raised for incorrect data types, missing required fields, or extra fields.
*   **QE Input Generation**: Create a reference ASE `Atoms` object for a simple system (e.g., H2 molecule). Write a unit test that calls `_generate_qe_input` and compares the resulting string to a pre-written, known-good QE input file. Test various parameter combinations (e.g., different `ecutwfc`, magnetism on/off).
*   **QE Output Parsing**: Create a collection of sample QE output files, including one for a successful run, one for a failed SCF convergence, and one with other potential errors. Write unit tests that call `_parse_qe_output` on each of these files and assert that it either extracts the correct numerical values or raises the appropriate exception.
*   **Database Client**: Use an in-memory SQLite database for testing. Write tests to verify that `write_calculation` correctly stores data and that `get_all_calculations` retrieves it in the expected format. Test edge cases like an empty database.
*   **Delta Learning Logic**: Create a simple dataset of DFT values and corresponding baseline potential values. Unit test the logic that calculates the residual, ensuring the subtraction is performed correctly.

### Integration Testing Approach

Integration tests will verify the connections and data flow between the different modules.
*   **Labeling to Database**: Create a test that uses a real, but very fast, Quantum Espresso calculation. A single H2 atom in a box is a good candidate. The test will:
    1.  Provide the structure to the `LabelingEngine`.
    2.  Run the engine, which will call the actual `pw.x` binary.
    3.  Assert that the engine completes without errors.
    4.  Query the test database and verify that the calculation results have been written correctly and match the output from QE.
*   **Database to Training**: Create a test that pre-populates a test database with a small, known dataset of labelled structures. The test will:
    1.  Run the `TrainingEngine`.
    2.  Assert that the engine runs to completion without errors.
    3.  Verify that a trained model file is created on disk.
    4.  (Optional) Load the model and perform a prediction on a known structure, asserting that the output is within a reasonable tolerance.
*   **Full End-to-End Test**: Combine the two integration tests above. This test will start with a minimal `input.yaml` and a structure, run the entire CLI command, and verify that it completes successfully, with the DFT data in the database and the final MLIP file on disk. This will be the key validation for the cycle's main objective.
