# CYCLE 01: UAT.md - Core Engine

## 1. Test Scenarios

This User Acceptance Test (UAT) is designed to verify the core functionality of the MLIP-AutoPipe system as developed in Cycle 01. The primary goal is to ensure that the system can successfully take a user-provided set of atomic structures, label them using an external DFT code (Quantum Espresso), and train a Machine Learning Interatomic Potential (MLIP) from the resulting data. This UAT acts as a tutorial for a foundational use case, demonstrating the system's capability as a robust data processing and model training engine.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C01-01** | **Linear Workflow Verification for Silicon** | **High** |

**Scenario UAT-C01-01: Linear Workflow Verification for Silicon**

This scenario focuses on a simple, well-understood, and computationally inexpensive material: crystalline Silicon (Si). The user will provide a small dataset of Si atomic structures, including the equilibrium diamond structure and a few structures with small, random atomic displacements. The goal is to confirm that the pipeline can execute its core functions—database population, DFT labeling, and MLIP training—sequentially and correctly. The user experience is centered on preparing two simple inputs (a configuration file and a file with atomic structures), running a single command, and observing the successful creation of a trained model file. This process validates that the fundamental data flow and component integrations are working as expected, providing a solid foundation for more complex features in subsequent cycles.

To facilitate this test, a Jupyter Notebook, `C01_UAT.ipynb`, will be provided. This notebook will guide the user through the following steps:
1.  **Setup:** Programmatically create the necessary input files:
    *   `exec_config_dump_c01.yaml`: The full configuration file specifying DFT and training parameters for Silicon.
    *   `initial_structures.xyz`: An XYZ file containing a few Si structures.
2.  **Execution:** Use the `!uv run mlip-pipe ...` command to execute the main CLI, pointing it to the generated configuration and an ASE database file.
3.  **Verification:** After the run completes, the notebook will contain cells to inspect the outputs:
    *   Connect to the output ASE database (`uat_c01.db`) and display its contents, showing that the structures are now marked as 'labeled' and contain energy and force data.
    *   Check for the existence of the final trained model file (`model.ace`).
This notebook-driven approach makes the UAT interactive and easy to follow, allowing the user to not only verify the outcome but also understand the inputs and outputs of the system clearly.

## 2. Behavior Definitions

**Feature:** Core Engine Workflow for MLIP Creation

**Scenario:** Successful execution of the labeling and training pipeline for a predefined set of structures.

*   **GIVEN** a clean working directory.
*   **GIVEN** a YAML configuration file (`exec_config_dump_c01.yaml`) specifying:
    *   System elements as `["Si"]`.
    *   DFT compute settings for Quantum Espresso.
    *   MLIP training settings for the ACE model with Delta Learning enabled.
*   **GIVEN** an XYZ file (`initial_structures.xyz`) containing three 8-atom silicon supercells:
    1.  The ideal diamond crystal structure.
    2.  The ideal structure with small (<0.1 Å) random displacements applied to each atom.
    3.  The ideal structure with larger (<0.3 Å) random displacements applied to each atom.
*   **GIVEN** an empty ASE database file (`uat_c01.db`).
*   **GIVEN** a mock of the Quantum Espresso executable that, when called, produces pre-calculated, valid output for each of the three silicon structures.

*   **WHEN** the user executes the command `uv run mlip-pipe --config-file exec_config_dump_c01.yaml --database-file uat_c01.db --input-file initial_structures.xyz`.

*   **THEN** the command should execute and exit with a success code (0).
*   **AND** the ASE database file `uat_c01.db` should be updated.
*   **AND** querying the `uat_c01.db` database should show exactly 3 rows.
*   **AND** each of the 3 rows in the database should have a `key_value_pair` indicating its state is 'labeled'.
*   **AND** each of the 3 rows in the database should contain non-zero values for energy, forces, and stress.
*   **AND** a new file named `model.ace` should be created in the working directory, representing the trained potential.
*   **AND** the standard output of the command should show log messages indicating the progression through the stages: "Initializing Workflow", "Starting Labeling Engine", "Labeling Complete", "Starting Training Engine", "Training Complete".
