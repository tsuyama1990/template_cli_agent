# CYCLE01/UAT.md

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for Cycle 1 of the MLIP-AutoPipe project. The focus of this cycle is to deliver a functional and reliable command-line interface (CLI) that can generate physically valid seed structures for various material systems. The UATs are designed to be run by a user (e.g., a computational materials scientist) to verify that the software meets their essential requirements and is easy to use. The primary tool for this UAT will be a Jupyter Notebook (`CYCLE01_UAT.ipynb`), which will provide a step-by-step, interactive guide for users to execute the CLI, inspect the outputs, and visualize the generated structures. This notebook-based approach will not only serve as a clear testing procedure but also as an introductory tutorial for new users.

| Scenario ID | Priority | Test Scenario Description                                                                |
|-------------|----------|------------------------------------------------------------------------------------------|
| UAT-C1-001  | High     | **Alloy Generation**: User wants to generate a set of 10 random FCC alloy structures of Gold (Au) and Copper (Cu) with a 50/50 composition. The user needs to verify that the output database contains 10 structures, each with the correct composition and no physical violations (e.g., overlapping atoms). This is a fundamental test of the core functionality. |
| UAT-C1-002  | High     | **Ionic Compound Generation**: User wants to generate 5 rutile-type Titanium Dioxide (TiO2) structures. The user needs to confirm that the generated structures have the correct 1:2 stoichiometry, are charge-neutral, and respect the minimum atomic distance constraints. This tests the generator for more complex, multi-element, charge-balanced systems. |
| UAT-C1-003  | Medium   | **Invalid Configuration Handling**: User attempts to run the pipeline with a faulty configuration file (e.g., an alloy composition that does not sum to 1.0, or a non-positive number for the structure count). The system should gracefully exit with a clear, user-friendly error message explaining what is wrong with the configuration, rather than crashing with an obscure traceback. This ensures the tool is robust and user-friendly. |
| UAT-C1-004  | Medium   | **Database Inspection and Visualization**: User wants to inspect the generated database from UAT-C1-001. Using the provided Jupyter Notebook, the user will connect to the output database, read the generated structures back into memory as ASE `Atoms` objects, and use a visualization library (like `nglview`) to visually inspect a few of the generated AuCu structures. This confirms that the output is in a standard, usable format and provides visual confirmation of success. |

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behavior for each UAT scenario. These will be implemented as guided steps within the `CYCLE01_UAT.ipynb` notebook.

---

**Scenario: UAT-C1-001 - Successful Generation of AuCu Alloy Structures**

*   **GIVEN** a user has a clean working directory.
*   **AND** the user creates a YAML configuration file named `config_aucu.yaml` specifying an FCC "alloy" system with elements `['Au', 'Cu']`, a composition of `{'Au': 0.5, 'Cu': 0.5}`, and `num_structures: 10`.
*   **AND** the configuration specifies an output database path of `aucu_structures.db`.
*   **WHEN** the user executes the command `mlip-autopipec generate config_aucu.yaml` from the terminal.
*   **THEN** the command should complete successfully with an exit code of 0.
*   **AND** a file named `aucu_structures.db` should be created in the working directory.
*   **AND** when the user inspects the database, it should contain exactly 10 atomic structures.
*   **AND** for each structure in the database, the chemical formula should be 'AuCu' (or an equivalent representation of a 50/50 mix, verifiable by checking atom counts).
*   **AND** for each structure, the distance between any two atoms must be greater than the default minimum distance specified in the tool's validation parameters.

---

**Scenario: UAT-C1-002 - Successful Generation of TiO2 Ionic Structures**

*   **GIVEN** a user has a clean working directory.
*   **AND** the user creates a YAML configuration file named `config_tio2.yaml` specifying an "ionic" system with stoichiometry `{'Ti': 1, 'O': 2}`, `num_structures: 5`, and a known crystal structure (e.g., rutile).
*   **AND** the configuration specifies an output database path of `tio2_structures.db`.
*   **WHEN** the user executes the command `mlip-autopipec generate config_tio2.yaml`.
*   **THEN** the command should complete successfully.
*   **AND** a file named `tio2_structures.db` should be created.
*   **AND** the database should contain exactly 5 structures.
*   **AND** each structure in the database must have a chemical formula of 'TiO2'.
*   **AND** each structure should be physically valid, with no overlapping atoms.

---

**Scenario: UAT-C1-003 - Graceful Failure with Invalid Configuration**

*   **GIVEN** a user has a clean working directory.
*   **AND** the user creates a YAML configuration file named `config_invalid.yaml` where the alloy composition is `{'Au': 0.5, 'Cu': 0.4}` (sums to 0.9 instead of 1.0).
*   **WHEN** the user executes the command `mlip-autopipec generate config_invalid.yaml`.
*   **THEN** the command should fail with a non-zero exit code.
*   **AND** the tool should print a clear error message to the console, such as "Configuration Error: Alloy compositions must sum to 1.0".
*   **AND** no database file should be created.

---

**Scenario: UAT-C1-004 - Inspection and Visualization of Generated Data**

*   **GIVEN** the `aucu_structures.db` file has been successfully created from scenario UAT-C1-001.
*   **AND** the user is working within the `CYCLE01_UAT.ipynb` Jupyter Notebook.
*   **WHEN** the user executes the notebook cell that contains Python code to connect to `aucu_structures.db` using ASE's `ase.db.connect`.
*   **AND** the user executes the cell to read all structures from the database into a list of `Atoms` objects.
*   **THEN** the list should contain 10 `Atoms` objects.
*   **AND** when the user executes a visualization cell (e.g., using `nglview.show_ase`), an interactive 3D view of the first AuCu structure should be displayed within the notebook, confirming the structure was written and can be read correctly.
