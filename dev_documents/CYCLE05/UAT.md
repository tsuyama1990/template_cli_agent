# Cycle 5 User Acceptance Testing (UAT)

This document outlines the User Acceptance Testing (UAT) for the final cycle of the MLIP-AutoPipe project. The focus of this cycle is the overall user experience, including the command-line interface, logging, error handling, and documentation.

## 1. Test Scenarios

These scenarios are designed to test the application from the perspective of a new user. They cover the main user interactions, the clarity of the feedback provided by the system, and the helpfulness of the documentation.

| Scenario ID | Scenario Name                                | Priority |
|-------------|----------------------------------------------|----------|
| UAT-C5-001  | Successful End-to-End Run from CLI           | High     |
| UAT-C5-002  | Correct Logging and Progress Indication      | High     |
| UAT-C5-003  | User-Friendly Error on Invalid Input         | High     |
| UAT-C5-004  | Usefulness of Documentation                  | Medium   |

---

**Scenario ID: UAT-C5-001**
- **Description:** This test is the ultimate validation of the system. A user with a valid `input.yaml` file should be able to run the entire pipeline from the command line with a single command and receive a final, trained MLIP. This tests the main success path of the application.
- **Priority:** High
- **Preconditions:**
    - The `mlip-pipe` application is installed in a clean environment.
    - All external dependencies (like Quantum Espresso) are installed and available in the PATH.
    - A valid `input.yaml` for a simple, fast-to-compute system (e.g., a silicon dimer) is available.

---

**Scenario ID: UAT-C5-002**
- **Description:** This test verifies that the system provides clear and useful feedback to the user during a run. The user should see informative log messages about the pipeline's progress and a progress bar for long-running tasks. This ensures the user knows what the application is doing and that it hasn't frozen.
- **Priority:** High
- **Preconditions:**
    - The user is running a standard pipeline execution from the CLI.

---

**Scenario ID: UAT-C5-003**
- **Description:** This test ensures that the application handles common user errors gracefully. If the user provides a syntactically correct YAML file but with an invalid value (e.g., a non-existent element), the system should provide a clear error message and exit, rather than crashing with a confusing stack trace.
- **Priority:** High
- **Preconditions:**
    - The user has created an `input.yaml` file with a semantic error (e.g., `composition: "Si_C"` but `elements: ["Si"]`).

---

**Scenario ID: UAT-C5-004**
- **Description:** This test validates the quality and clarity of the user documentation. A new user, with a background in materials science but not in the specifics of this tool, should be able to follow the `README.md` to install the software, understand the basic workflow, and run the provided tutorial.
- **Priority:** Medium
- **Preconditions:**
    - A new user is given access to the project's repository.
    - The user has the necessary base environment (e.g., Python, `uv`).

## 2. Behavior Definitions

---

**Scenario: UAT-C5-001 - Successful End-to-End Run from CLI**

```gherkin
GIVEN a user has created a valid 'input.yaml' for a silicon system
AND the user is in a terminal with the 'mlip-pipe' command available

WHEN the user executes the command "uv run mlip-pipe input.yaml"

THEN the application should start and display an initial log message
AND the application should run through all the cycles (Generation, Labeling, Training, etc.) without crashing
AND the application should eventually finish and print a "Workflow completed successfully" message
AND a directory containing the final MLIP model file and logs should be present
AND the command should exit with a status code of 0
```

---

**Scenario: UAT-C5-002 - Correct Logging and Progress Indication**

```gherkin
GIVEN a user has started a standard pipeline run from the CLI

WHEN the system begins a long-running task, such as the surrogate MD exploration

THEN a progress bar should be displayed on the console, showing the progress of the task
AND the console should display INFO-level log messages at the start and end of each major stage (e.g., "Starting Module B: Explorer & Sampler", "Module B finished.")
AND if the user runs the same command with the "--log-file run.log" option, a file named "run.log" should be created containing the full log output
```

---

**Scenario: UAT-C5-003 - User-Friendly Error on Invalid Input**

```gherkin
GIVEN a user has created an 'input.yaml' file, but it refers to an element in the composition that is not declared in the elements list
  """
  system:
    elements: ["Fe", "O"]
    composition: "Fe2O3_X" # Contains an invalid 'X'
  """
WHEN the user executes "uv run mlip-pipe input.yaml"

THEN the application should exit within a few seconds
AND it should print a clear, user-friendly error message to the console, such as "Error: Composition 'Fe2O3_X' contains an element 'X' not declared in the elements list."
AND the command should exit with a non-zero status code
```

---

**Scenario: UAT-C5-004 - Usefulness of Documentation**

```gherkin
GIVEN a new user has cloned the project repository

WHEN the user reads the 'README.md' file

THEN the user should be able to find clear, step-by-step instructions for installing the software and its dependencies
AND the user should be able to find and run a "Quick Start" example command
AND by following the tutorial in the 'docs/' directory, the user should be able to successfully run a more complex example and understand the structure of the outputs
```
