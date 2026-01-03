# CYCLE 01: USER ACCEPTANCE TESTING (UAT) - Core CLI Pipeline

This document outlines the User Acceptance Testing (UAT) scenarios for the core command-line functionality delivered in Cycle 1 of the MLIP-AutoPipe project. The goal of this UAT is to ensure that a user can successfully and intuitively use the CLI to generate a valid materials database from a simple, declarative configuration file.

## 1. Test Scenarios

| Scenario ID | Title                                       | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C1-001  | Successful End-to-End Pipeline Run for an Alloy | High     |
| UAT-C1-002  | Handling of Invalid Configuration File      | High     |
| UAT-C1-003  | Pipeline Resumption (Idempotency)           | Medium   |

### Scenario UAT-C1-001: Successful End-to-End Pipeline Run for an Alloy

**(Min 300 words)**
This is the primary "happy path" scenario. It verifies that a user can take a standard, well-formed configuration file, run the main pipeline command, and receive a valid, populated ASE database as output without any errors. This test is crucial as it validates the entire integrated software stack, from the initial parsing of the user's input to the final write in the database. The user experience should be simple and predictable: the command is executed, provides some feedback that it is running, and upon completion, the specified output file exists and contains valid data.

For this test, we will provide a sample Jupyter Notebook, `UAT-C1-001.ipynb`. This notebook will serve as both the testing script and as a piece of user documentation. It will guide the user through the following steps:
1.  **Setup:** Programmatically create a temporary working directory.
2.  **Configuration:** Define a valid YAML configuration as a Python string and write it to a `config.yml` file in the temporary directory. The configuration will specify a simple binary alloy, such as Copper-Gold (CuAu), with a request for a small number of initial structures (e.g., 5) and a very short exploration phase (e.g., 10 MD steps).
3.  **Execution:** Use Python's `subprocess` module to call the CLI command (`mlip-autopipec run-pipeline --config config.yml`). The notebook will capture and display the standard output and standard error from the command, allowing the user to see the application's logging in real-time.
4.  **Verification:** After the command completes successfully, the notebook will contain cells to inspect the output. It will check for the existence of the database file (`results.db`) specified in the config. It will then use the `ase.db` library to connect to this database and perform several assertions:
    *   Verify that the database contains exactly 5 rows.
    *   Query each row and check that it has been marked as 'labeled'.
    *   Extract one of the `Atoms` objects and verify that it contains both Copper and Gold atoms.
    *   Check that the `data` field for each row contains the expected results (energy, forces, stress).
This notebook provides a seamless and interactive way for the user to confirm that the core promise of the software is met.

### Scenario UAT-C1-02: Handling of Invalid Configuration File

**(Min 300 words)**
This scenario tests the system's robustness and user-friendliness in the face of common errors, specifically an invalid configuration file. A key aspect of a good user experience is providing clear, actionable feedback when the user makes a mistake. The system should not crash with an obscure traceback; instead, it should gracefully exit and inform the user exactly what is wrong with their input. This test validates the Pydantic-based schema validation layer.

Similar to the first test, this will be demonstrated in a Jupyter Notebook, `UAT-C1-002.ipynb`. This notebook will present the user with several examples of incorrect configurations and show the expected error output for each.
1.  **Setup:** The notebook will again create a temporary directory for its operations.
2.  **Execution (Test Case A - Typo in Key):** It will create a `config.yml` with a deliberate typo in a key (e.g., `temprature_K` instead of `temperature_K`). It will then run the CLI command and assert that the program exits with a non-zero status code and prints an error message that clearly states that `temprature_K` is not a valid field and suggests the correct spelling, `temperature_K`.
3.  **Execution (Test Case B - Invalid Value Type):** It will create a `config.yml` where a value has the wrong type (e.g., `temperature_K: "hot"` instead of a number). It will run the command and assert that the error message clearly indicates that the value for `temperature_K` must be a valid number.
4.  **Execution (Test Case C - Value Out of Range):** It will create a `config.yml` with a valid type but an invalid value (e.g., `temperature_K: -100`), which violates the Pydantic model's `gt=0` constraint. It will run the command and assert that the error message clearly states that the temperature must be greater than zero.

This notebook demonstrates to the user that the system is designed to help them, validating the integrity of their inputs and providing immediate, helpful feedback, which is a cornerstone of a positive user experience.

### Scenario UAT-C1-03: Pipeline Resumption (Idempotency)

**(Min 300 words)**
This scenario tests a more advanced but critical feature for a potentially long-running scientific workflow: idempotency and the ability to resume. If a user starts a pipeline that takes several hours and it gets interrupted, they should be able to simply run the exact same command again without corrupting their results or re-doing completed work. This test ensures the `WorkflowOrchestrator` correctly checks the state of the database before executing tasks.

This test will also be presented in a Jupyter Notebook, `UAT-C1-003.ipynb`, to guide the user.
1.  **Setup:** The notebook sets up a temporary directory and a valid `config.yml` for generating 10 structures.
2.  **Initial Run:** It runs the CLI command for the first time. After it completes, the verification step connects to the database and confirms that 10 structures have been generated and labeled.
3.  **Simulating Interruption (Advanced):** To make the test more explicit, the notebook will then manually connect to the database and "un-label" half of the entries, simulating a state where the process was interrupted midway through the exploration phase.
4.  **Second Run:** The notebook will execute the *exact same* CLI command a second time.
5.  **Verification:** The core of the UAT is in this step. The notebook will:
    *   Show the user the log output from the second run. The output should indicate that the system found 10 initial structures already present and did not re-run the generation step. It should show that it only performed the "exploration" task for the 5 structures that were manually "un-labeled".
    *   Connect to the database and verify that all 10 rows are now correctly marked as 'labeled'.
    *   Check the modification times of the database entries or use other markers to confirm that the 5 already-completed rows were not modified during the second run.

This scenario demonstrates to the user that the system is reliable and efficient, saving them valuable time and computational resources by intelligently resuming from the last known good state.

## 2. Behavior Definitions

**GIVEN** a user has a clean, empty directory
**AND** the user creates a valid `config.yml` file specifying a CuAu alloy, 5 initial structures, and a database path of `my_alloy.db`
**WHEN** the user executes the command `mlip-autopipec run-pipeline --config config.yml` in the terminal
**THEN** the command should exit with a status code of 0
**AND** a file named `my_alloy.db` should exist in the directory
**AND** the `my_alloy.db` file should be a valid SQLite database
**AND** when connected to with `ase.db`, the database should contain exactly 5 rows
**AND** each row in the database should have a `labeled` key set to `True`
**AND** each row, when read as an `Atoms` object, should contain both `Cu` and `Au` elements.

---

**GIVEN** a user has a clean, empty directory
**AND** the user creates a `config.yml` file with the content `system:\n  temperature_K: -50.0`
**WHEN** the user executes the command `mlip-autopipec run-pipeline --config config.yml`
**THEN** the command should exit with a non-zero status code
**AND** the command should print an error message to the console
**AND** the error message should contain the phrase "temperature_K"
**AND** the error message should contain the phrase "must be greater than 0".

---

**GIVEN** a user has already successfully run the pipeline for a CuAu alloy, resulting in a populated `my_alloy.db` file with 5 labeled rows
**WHEN** the user executes the *exact same* command `mlip-autopipec run-pipeline --config config.yml` a second time
**THEN** the command should exit with a status code of 0
**AND** the command's log output should indicate that the generation step was skipped
**AND** the command's log output should indicate that zero new exploration tasks were run
**AND** the `my_alloy.db` file should still contain exactly 5 rows.
