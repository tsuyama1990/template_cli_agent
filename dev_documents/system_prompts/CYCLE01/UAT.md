# Cycle 1 User Acceptance Testing (UAT): Core Pipeline Functionality

This document outlines the User Acceptance Testing (UAT) scenarios for the successful completion of Cycle 1. The goal of this cycle is to deliver a foundational, end-to-end command-line tool that can generate a set of simple alloy structures and store them in a database. This UAT is designed to be executed by a user to confirm that this core functionality is met and provides a stable base for future development.

## 1. Test Scenarios

The primary method for conducting this UAT will be through a Jupyter Notebook. This approach allows for an interactive and educational experience, where the user can not only run the tests but also inspect the outputs directly, providing a clear and tangible demonstration of the system's capabilities. It also serves as a practical tutorial for a new user.

| Scenario ID | Scenario Name                               | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C1-01   | Successful End-to-End Pipeline Execution    | High     |
| UAT-C1-02   | Handling of Invalid Configuration           | High     |
| UAT-C1-03   | Verifying Database Content and Integrity    | High     |

### UAT-C1-01: Successful End-to-End Pipeline Execution

**Description:**
This is the "happy path" scenario. The user will create a standard configuration file for a common binary alloy (e.g., Copper-Gold, CuAu) and run the main command-line tool. The purpose is to verify that the tool executes without errors and produces the expected output file. This test confirms that all the basic components—CLI, configuration parsing, orchestration, generation, and storage—are correctly integrated and functional.

**Execution via Jupyter Notebook (`UAT_Cycle1.ipynb`):**
The notebook will be structured with Markdown cells explaining each step and code cells to be executed.
1.  **Setup:** The first code cell will programmatically create a YAML configuration file named `config_cuau.yaml` in the notebook's directory. This file will define parameters for generating 10 structures of a Cu-Au alloy.
2.  **Execution:** The next cell will use the `!mlip-autopipec --config config_cuau.yaml` shell command to run the pipeline. The output of the command will be displayed directly in the notebook.
3.  **Verification:** The user will visually inspect the output for any error messages. A successful run should print a confirmation message, e.g., "Pipeline executed successfully. 10 structures saved to cuau_structures.db". A final code cell will check for the existence of the `cuau_structures.db` file, providing a clear "PASSED" or "FAILED" message.

### UAT-C1-02: Handling of Invalid Configuration

**Description:**
This scenario tests the system's robustness and user-friendliness when provided with incorrect or nonsensical input. The Pydantic-based configuration design should catch these errors gracefully and provide informative feedback to the user. We will test two failure modes: a syntactically incorrect YAML file and a semantically incorrect configuration (e.g., a negative lattice constant).

**Execution via Jupyter Notebook (`UAT_Cycle1.ipynb`):**
1.  **Syntactic Error:** A cell will create a deliberately broken YAML file (`config_broken.yaml`) with improper indentation. The next cell will attempt to run the pipeline with this file.
2.  **Verification:** The user will verify that the tool does *not* crash with an unhandled exception. Instead, it should display a clear error message from the Hydra framework indicating that the configuration file could not be parsed.
3.  **Semantic Error:** Another cell will create a YAML file (`config_invalid.yaml`) that is syntactically correct but contains a physically impossible value, such as `lattice_constant: -3.5`.
4.  **Verification:** When the pipeline is run with this file, the user will verify that the application exits gracefully and prints a user-friendly error message from Pydantic, explaining exactly which field is invalid (e.g., "Validation Error: `lattice_constant` must be greater than 0").

### UAT-C1-03: Verifying Database Content and Integrity

**Description:**
This scenario ensures that the data stored in the database is correct and matches the input configuration. It goes beyond simply checking for the file's existence and dives into the content itself. This test is crucial for confirming that the data persistence layer is functioning as expected.

**Execution via Jupyter Notebook (`UAT_Cycle1.ipynb`):**
This test will build upon the successful run from UAT-C1-01.
1.  **Database Connection:** A code cell will use the `ase.db` library to connect to the `cuau_structures.db` file generated in the first scenario.
2.  **Content Verification:**
    *   A query will be run to count the number of rows in the database. The user will verify that this count is exactly 10.
    *   A cell will retrieve the first record from the database as an `ase.Atoms` object.
    *   The user will inspect the properties of this object. A series of `assert` statements will programmatically check that the atoms object contains only 'Cu' and 'Au' symbols, and that the total number of atoms is correct based on the supercell size defined in the configuration.
    *   A final cell can provide a simple 3D visualization of the retrieved structure using a library like `nglview`, giving the user immediate visual confirmation that the generated structure is a reasonable-looking alloy.

## 2. Behavior Definitions

The expected behavior of the system is defined below using Gherkin-style syntax.

**Scenario: Successful End-to-End Pipeline Execution**

*   **GIVEN** a valid YAML configuration file `config_cuau.yaml` exists, specifying the generation of 10 Copper-Gold (CuAu) alloy structures and an output database path of `cuau_structures.db`.
*   **WHEN** the user executes the command `mlip-autopipec --config config_cuau.yaml` from the terminal.
*   **THEN** the command should execute without raising any errors and should complete with an exit code of 0.
*   **AND** a success message should be printed to the console, indicating that 10 structures have been successfully generated and saved.
*   **AND** a new file named `cuau_structures.db` should be present in the current working directory.

**Scenario: Handling of a Semantically Invalid Configuration**

*   **GIVEN** a YAML configuration file `config_invalid.yaml` exists, which is syntactically correct but contains a semantic error (e.g., `lattice_constant: -3.5`).
*   **WHEN** the user executes the command `mlip-autopipec --config config_invalid.yaml`.
*   **THEN** the application should terminate gracefully with a non-zero exit code.
*   **AND** a clear, user-friendly error message should be printed to the console, specifying which configuration parameter is invalid and why (e.g., "Error: `lattice_constant` must be a positive value").
*   **AND** no database file should be created.

**Scenario: Verifying Database Content and Integrity**

*   **GIVEN** the pipeline has been successfully executed with a configuration to generate 10 CuAu structures in a 5x5x5 supercell of an FCC lattice.
*   **AND** the output database `cuau_structures.db` has been created.
*   **WHEN** a user or script connects to `cuau_structures.db` using the ASE database library.
*   **THEN** a query for the total number of structures in the database should return exactly 10.
*   **AND** when the first structure is retrieved from the database, it should be a valid `ase.Atoms` object.
*   **AND** this `ase.Atoms` object should contain a total of 500 atoms (5x5x5 supercell * 4 atoms per FCC unit cell).
*   **AND** the chemical symbols of all atoms in the object must be either 'Cu' or 'Au'.
