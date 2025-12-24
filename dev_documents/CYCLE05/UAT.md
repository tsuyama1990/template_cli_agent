# CYCLE05 User Acceptance Testing (UAT)

## 1. Test Scenarios

This UAT plan focuses on verifying the final user experience, system robustness, and the quality of the documentation. It is less about the scientific output and more about the usability and polish of the tool.

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C05-001 | New User Onboarding Experience              | High     |
| UAT-C05-002 | Clarity of CLI Output and Progress          | High     |
| UAT-C05-003 | Robustness to Invalid User Configuration    | High     |
| UAT-C05-004 | Robustness to Missing External Dependencies | High     |
| UAT-C05-005 | Documentation Accuracy and Completeness     | Medium   |

---

**Scenario UAT-C05-001: New User Onboarding Experience**

*   **Description**: This test simulates the experience of a brand new user attempting to install and run the software for the first time, guided only by the documentation.
*   **Preconditions**:
    *   A clean machine or virtual environment with only Python, `git`, and `uv` installed.
    *   The user has access to the project's `README.md` and the `docs/` directory.
*   **Acceptance Criteria**:
    *   The user must be able to successfully install the MLIP-AutoPipe package and all its dependencies using the commands provided in the `docs/getting_started.md`.
    *   The user must be able to run the simple "Silicon" example described in the getting started guide.
    *   The example run must complete successfully without errors.
    *   The user should report that the instructions were clear, easy to follow, and sufficient to get the example working.

---

**Scenario UAT-C05-002: Clarity of CLI Output and Progress**

*   **Description**: This test evaluates the effectiveness of the new user interface elements, such as progress bars and structured logging.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation.
    *   A standard workflow is initiated for a system that will take at least 1-2 minutes to run.
*   **Acceptance Criteria**:
    *   During the DFT calculation phase, a progress bar must be displayed, showing the number of completed structures, the total, and an estimated time remaining.
    *   During the surrogate MD exploration, a progress bar or spinner should indicate that the system is busy.
    *   During the MLIP training phase, a spinner or a static message like "Training MLIP model..." should be displayed.
    *   The console output should be clean and well-structured, with clear headings for different stages of the workflow (e.g., "--- Starting DFT Labeling ---").
    *   A detailed `autoprope.log` file must be created, containing much more verbose, timestamped information suitable for debugging.

---

**Scenario UAT-C05-003: Robustness to Invalid User Configuration**

*   **Description**: This test verifies that the system provides clear, actionable feedback when the user provides a faulty configuration file.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation.
    *   Several invalid `input.yaml` files are prepared:
        1.  One with a typo in a key (e.g., `compostion:` instead of `composition:`).
        2.  One with an invalid value (e.g., `temperature: "high"` instead of a number).
        3.  One that is not valid YAML.
*   **Acceptance Criteria**:
    *   In all cases, the program must exit gracefully and *not* display a long, cryptic Python stack trace.
    *   For the typo, the error message should say something like: "Configuration error: Unknown field 'compostion'. Did you mean 'composition'?".
    *   For the invalid value, the Pydantic validation error should be presented clearly, e.g., "Validation error in 'temperature': expected a number, but received 'high'".
    *   For the invalid YAML, the message should indicate a parsing error, e.g., "Error: The configuration file is not a valid YAML file."

---

## 2. Behavior Definitions

**Behavior for UAT-C05-001: New User Onboarding Experience**

```gherkin
Feature: New User Experience
  As a new user,
  I want to be able to install the software and run a basic example by following the official documentation,
  So that I can quickly get started and evaluate the tool.

  Scenario: Following the Getting Started guide
    GIVEN I have a clean environment with Python and uv
    AND I have cloned the project repository
    WHEN I follow the installation and setup instructions in `docs/getting_started.md`
    AND I run the command provided for the introductory Silicon example
    THEN the installation should complete without errors
    AND the example workflow should run to completion and produce an MLIP model
    AND I should not need any information outside of the documentation to succeed.
```

---

**Behavior for UAT-C05-004: Robustness to Missing External Dependencies**

```gherkin
Feature: Graceful Failure for Missing Dependencies
  As a user,
  I want the system to fail fast with a clear message if a required external program is not installed,
  So that I can easily diagnose and fix my environment setup.

  Scenario: Attempting to run a workflow without Quantum Espresso
    GIVEN a functional MLIP-AutoPipe installation
    BUT the Quantum Espresso `pw.x` executable is not in the system PATH or the configured location
    WHEN I try to run any workflow
    THEN the program should terminate almost immediately
    AND it should display a clear, user-friendly error message, such as "Error: Quantum Espresso executable 'pw.x' not found. Please ensure it is installed and in your PATH."
    AND it must not produce a "FileNotFoundError" stack trace.
```

---

**Behavior for UAT-C05-005: Documentation Accuracy and Completeness**

```gherkin
Feature: Accurate and Helpful Documentation
  As a user,
  I want to find clear explanations of advanced features in the documentation,
  So that I can customize the workflow for my specific research needs.

  Scenario: Consulting the documentation for an advanced feature
    GIVEN a functional MLIP-AutoPipe installation
    WHEN I read the tutorial page on "Configuring the Active Learning Loop"
    THEN the page should clearly explain the meaning of the `uncertainty_threshold` parameter
    AND it should provide a good example of how to set up the `simulation` section in the config file
    AND the example configuration provided in the documentation should be valid and work correctly if I use it in a run.
```
