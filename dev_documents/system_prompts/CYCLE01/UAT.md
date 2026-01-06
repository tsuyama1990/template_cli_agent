# User Acceptance Test (UAT) Plan: CYCLE 01

This document outlines the User Acceptance Testing scenarios for the core features implemented in Cycle 1 of the MLIP-AutoPipe project. The goal of this UAT is to verify that the foundational pipeline, which handles configuration parsing, initial structure generation, and database storage, works as expected from an end-user's perspective. It is about ensuring the software is not only functionally correct but also usable, robust, and reliable.

It is highly recommended to perform this UAT using a Jupyter Notebook. This interactive format is an ideal environment for this kind of testing as it will allow the user to execute the steps one by one, inspect the outputs in real-time, and gain a clear, hands-on understanding of the system's core functionality. The notebook will serve a dual purpose: it is both a rigorous test script for validating the software and an interactive tutorial for new users who want to learn how to use the application. This approach makes the testing process itself a valuable asset for the project's documentation and user onboarding.

## 1. Test Scenarios

| ID    | Priority | Test Scenario                                   |
| :---- | :------- | :---------------------------------------------- |
| UAT-1 | High     | Verify successful pipeline run with a valid configuration. |
| UAT-2 | High     | Verify graceful pipeline failure with an invalid configuration. |
| UAT-3 | Medium   | Verify the contents and scientific integrity of the generated database. |

### Scenario Details

#### **UAT-1: Successful pipeline run with a valid configuration.** (High Priority)
This is the most fundamental "happy path" scenario. Its purpose is to test the primary, end-to-end functionality of the core pipeline in its ideal state. The user should be able to take a standard, well-formed configuration file, run the main command-line tool, and see the process complete smoothly and without any errors. A successful run of this test provides the foundational confidence that the application is working as designed. It confirms that all the basic components developed in Cycle 1—the CLI, the configuration parsing and validation, the structure generation logic, and the database storage engine—are correctly integrated and functioning in concert. This scenario is designed to deliver a moment of amazement for the user: the simplicity of going from a simple, human-readable text file that describes a chemical system to a fully populated, machine-readable database of atomic structures. This successful test is the first proof point that the application can reliably automate a previously complex and manual task. It validates the core value proposition of the software.

#### **UAT-2: Verify graceful pipeline failure with an invalid configuration.** (High Priority)
This scenario is just as crucial as the happy path, as it tests the application's robustness and user-friendliness in the face of common user errors. A high-quality application should not only work correctly with perfect input but should also fail gracefully and, most importantly, informatively when given invalid input. This test will involve attempting to run the CLI with a deliberately incorrect configuration file. Examples of invalid configurations could include a misspelled key, a chemical composition that doesn't sum to 1.0, or a negative number of requested structures. The expected outcome is not a confusing crash with a long stack trace, but a clear, easy-to-understand error message that precisely identifies the problem in the configuration file and suggests how to fix it. This test directly validates the effectiveness of the Pydantic-based schema validation. A successful outcome here demonstrates the robustness of the system and significantly improves the overall user experience by empowering users to self-correct their inputs, making the software feel helpful rather than brittle.

#### **UAT-3: Verify the contents and scientific integrity of the generated database.** (Medium Priority)
This scenario "closes the loop" and ensures that the output of the pipeline is not just created, but is scientifically correct and matches the user's intent. It is not enough for the pipeline to run without errors; it must produce the right results. After a successful run (from UAT-1), the user will use a script, ideally within the UAT Jupyter Notebook, to connect to and inspect the generated SQLite database. They will perform several specific and critical checks to validate the integrity of the output data. These checks include:
1.  Verifying that the number of atomic structures stored in the database exactly matches the number requested in the configuration file. This is a basic sanity check.
2.  Extracting one of the atomic structures and verifying that its chemical composition (i.e., the ratio of different elements) matches the composition requested in the configuration.
3.  Ensuring that the structures are physically plausible by checking that there are no unnaturally close atoms (atomic overlaps). This confirms that the internal validation logic within the generator is working correctly.
This test provides tangible, objective proof that the software is producing the scientifically correct and expected results, which is the ultimate and most important goal of this scientific tool.

## 2. Behavior Definitions

These definitions are written in the Gherkin style to clearly and unambiguously outline the expected behavior, preconditions, and post-conditions for each scenario.

---

**Scenario: UAT-1 - Successful pipeline run with a valid configuration**

*   **GIVEN** I have a valid YAML configuration file named `config_valid.yaml` in my current directory.
*   **AND** this file specifies the following configuration:
    ```yaml
    system:
      generator_type: "alloy"
      elements: ["Si", "C"]
      composition: {"Si": 0.5, "C": 0.5}
      num_structures: 10
    database_path: "uat_pipeline.db"
    ```
*   **AND** the file `uat_pipeline.db` does not currently exist in my directory.
*   **WHEN** I open my terminal and execute the command `mlip-autopipec run-pipeline --config-path config_valid.yaml`.
*   **THEN** the command should complete in a reasonable amount of time and exit with a success code of 0.
*   **AND** the terminal should display a clear, user-friendly success message, for example, "Pipeline completed successfully."
*   **AND** I should see a new file named `uat_pipeline.db` in my current directory.

---

**Scenario: UAT-2 - Verify graceful pipeline failure with an invalid configuration**

*   **GIVEN** I have an invalid YAML configuration file named `config_invalid.yaml` in my current directory.
*   **AND** this file contains a composition where the values sum to 0.9 instead of 1.0, which is physically incorrect:
    ```yaml
    system:
      generator_type: "alloy"
      elements: ["Si", "C"]
      composition: {"Si": 0.5, "C": 0.4} # Invalid: 0.5 + 0.4 = 0.9
      num_structures: 10
    database_path: "uat_invalid.db"
    ```
*   **WHEN** I open my terminal and execute the command `mlip-autopipec run-pipeline --config-path config_invalid.yaml`.
*   **THEN** the command should fail quickly with a non-zero exit code.
*   **AND** the terminal should display a clear, specific error message that explicitly states that the composition values must sum to 1.0.
*   **AND** the application should not create an empty or partial database file; no file named `uat_invalid.db` should be created.

---

**Scenario: UAT-3 - Verify the contents and integrity of the generated database**

*   **GIVEN** I have successfully completed the pipeline run from scenario UAT-1.
*   **AND** the database file `uat_pipeline.db` exists in my current directory.
*   **WHEN** I execute a Python script (or a Jupyter Notebook cell) that uses the `ase.db` library to connect to the `uat_pipeline.db` file.
*   **THEN** the script should be able to successfully connect to the database.
*   **AND** a query for the number of rows in the database should return exactly `10`.
*   **AND** when I select and read the first row from the database, it should be a valid `ase.Atoms` object.
*   **AND** the chemical formula of this `Atoms` object should correctly represent a 50/50 composition of Silicon and Carbon (e.g., for a 64-atom supercell, it should contain 32 Si and 32 C atoms).
*   **AND** a check of the minimum interatomic distances in the structure should confirm that all atoms are separated by a physically reasonable distance (e.g., greater than 1.0 Angstrom), proving that the internal physical validation was applied correctly.
