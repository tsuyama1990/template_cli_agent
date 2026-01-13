# UAT.md: Cycle 1 - Core Framework and Basic Pipeline

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 1 of the MLIP-AutoPipe project. The goal of this UAT is to verify that the core functionality of the Minimum Viable Product (MVP) meets the user's requirements and provides a solid foundation for future development. The tests are designed to be executed by a user or QA analyst to confirm that the command-line tool is functional, usable, and produces the expected outputs. The primary recommended method for performing this UAT is via a Jupyter Notebook, which allows for interactive execution of the tool and detailed inspection of the results.

| Scenario ID | Scenario Name                     | Priority |
| :---------- | :-------------------------------- | :------- |
| UAT-C1-001  | Successful End-to-End Pipeline Run | High     |
| UAT-C1-002  | Handling of Invalid Configuration | High     |
| UAT-C1-003  | Verification of Database Output   | High     |

### **UAT-C1-001: Successful End-to-End Pipeline Run** (Priority: High)

**Description:**
This is the primary "happy path" test scenario. It verifies that a user can successfully run the entire data generation pipeline using a valid, well-formed configuration file. The test ensures that the CLI tool executes without errors and that all four stages of the pipeline (Generation, Exploration, Sampling, Storage) are completed in the correct sequence. The user experience is central here; the tool should provide clear feedback during its execution, such as progress bars or status messages, indicating what stage is currently running. The successful creation of the output database file is the key indicator of a successful run. This test simulates the most common user interaction and validates the core value proposition of the tool. A Jupyter Notebook (`UAT_Cycle1_Validation.ipynb`) will be provided to guide the user through creating the configuration, running the tool from the notebook, and performing the initial checks.

### **UAT-C1-002: Handling of Invalid Configuration** (Priority: High)

**Description:**
A robust tool must handle user error gracefully. This scenario tests the system's resilience to invalid or malformed configuration files. The goal is to ensure that the tool provides clear, understandable error messages when it encounters problems with its input, rather than crashing with an obscure traceback. This is crucial for a positive user experience, as it helps the user to quickly diagnose and fix their configuration. This UAT will involve several sub-tests, each with a different type of invalid configuration. Examples include a YAML file with incorrect syntax, a configuration with a missing required field (e.g., `temperature_k`), a configuration with a logically invalid value (e.g., a negative number of structures), or a composition that does not sum to 1.0. The expected outcome is that the tool exits with a non-zero status code and prints a message that clearly explains which validation rule failed.

### **UAT-C1-003: Verification of Database Output** (Priority: High)

**Description:**
This scenario focuses on the integrity and correctness of the final output: the ASE database. It's not enough for the pipeline to simply run without errors; the data it produces must be correct and consistent with the input configuration. This test involves a deeper inspection of the `results.db` file generated in scenario UAT-C1-001. Using the provided Jupyter Notebook, the user will connect to the database and perform a series of checks. These include: verifying that the number of atomic structures in the database matches the `num_samples` parameter from the configuration file, checking that each structure in the database has the correct atomic species and composition, and ensuring that essential metadata such as `energy` and `forces` are present for each structure. This test provides the ultimate validation that the software is not just technically functional but also scientifically correct for the scope of Cycle 1.

## 2. Behavior Definitions

These behaviors are defined in the Gherkin style (GIVEN/WHEN/THEN) to provide clear, unambiguous descriptions of the expected system behavior for each test scenario.

---

### **Scenario: Successful End-to-End Pipeline Run**

*   **GIVEN** A clean working directory.
*   **AND** The user has a valid YAML configuration file named `config.yml` that specifies a binary Cu-Au alloy, requests 4 initial structures, an EMT-based MD exploration at 300K for 10 steps, and asks for 2 final samples using random sampling.
*   **WHEN** The user executes the command `mlip-autopipec run --config-path config.yml` in the terminal.
*   **THEN** The application should start without errors.
*   **AND** The application should display progress indicators for the Generation, Exploration, Sampling, and Storage stages.
*   **AND** The application should exit with a status code of 0, indicating success.
*   **AND** A new file named `results.db` should be present in the working directory.

---

### **Scenario: Handling of Configuration with Invalid Value**

*   **GIVEN** A clean working directory.
*   **AND** The user has a YAML configuration file named `invalid_config.yml` where the `num_structures` parameter is set to -1, which is a physically impossible value.
*   **WHEN** The user executes the command `mlip-autopipec run --config-path invalid_config.yml` in the terminal.
*   **THEN** The application should exit with a non-zero status code, indicating an error.
*   **AND** The application should print a clear error message to the console.
*   **AND** The error message should state that the `num_structures` field must be greater than 0, referencing the validation rule that failed.
*   **AND** No `results.db` file should be created.

---

### **Scenario: Handling of Configuration with Malformed YAML**

*   **GIVEN** A clean working directory.
*   **AND** The user has a file named `malformed_config.yml` that contains invalid YAML syntax (e.g., incorrect indentation).
*   **WHEN** The user executes the command `mlip-autopipec run --config-path malformed_config.yml` in the terminal.
*   **THEN** The application should exit with a non-zero status code, indicating an error.
*   **AND** The application should print a clear error message indicating that it failed to parse the YAML file.
*   **AND** No `results.db` file should be created.

---

### **Scenario: Verification of Database Content**

*   **GIVEN** A `results.db` file has been successfully generated by a pipeline run (as in UAT-C1-001) using a configuration that requested 2 final samples of a Cu-Au alloy.
*   **WHEN** The user opens the `results.db` file using a tool or script that can read ASE databases (e.g., the provided Jupyter Notebook).
*   **THEN** The database should contain exactly 2 rows (atomic configurations).
*   **AND** When inspecting a row from the database, the atoms object should contain only 'Cu' and 'Au' elements.
*   **AND** Each row in the database must have key-value pairs for `energy` and `forces`.
*   **AND** The `forces` value should be a NumPy array of shape (N, 3), where N is the number of atoms in the structure.
