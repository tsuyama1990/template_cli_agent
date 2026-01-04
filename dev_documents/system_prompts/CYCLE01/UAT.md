# CYCLE01/UAT.md

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 1 of the MLIP-AutoPipe project. The primary focus of this cycle is to validate the core command-line functionality for generating and storing initial atomic structures. The tests are designed to be executed by a user to confirm that the software meets its basic requirements and provides a satisfactory and robust user experience. To facilitate this, the primary UAT artifact will be a single, well-documented Jupyter Notebook (`UAT_Cycle01.ipynb`). This notebook will provide a narrative walkthrough, combining explanations with executable code cells that allow the user to interactively run the tests and verify the outcomes themselves. This approach is superior to a simple text document as it provides a live, hands-on validation experience.

| Scenario ID | Test Scenario | Priority |
| :--- | :--- | :--- |
| UAT-001 | Generate a set of simple binary alloy structures via the CLI | High |
| UAT-002 | Handle invalid and malformed configuration files gracefully | High |
| UAT-003 | Verify the contents and structural integrity of the generated database | Medium |

### UAT-001: Generate a set of simple binary alloy structures via the CLI

This test scenario represents the primary "happy path" for the application in Cycle 1. Its purpose is to confirm that a user can successfully execute the end-to-end workflow with a valid configuration and receive the expected output. A successful test here demonstrates that the core logic of configuration parsing, structure generation, and database storage is functioning correctly. The Jupyter Notebook will contain a dedicated section for this scenario. It will begin by programmatically creating a valid YAML configuration file. This file will specify a simple, well-understood binary alloy system, for instance, a request for 10 structures of a Copper-Gold (CuAu) alloy in the Face-Centered Cubic (fcc) crystal structure. The notebook will then provide a code cell that uses the `!` shell command to execute the `mlip-autopipec run-pipeline --config ...` command. The user can run this cell and see the real-time output from the CLI. The final part of the test will involve using the ASE library within the notebook to connect to the newly created database file and perform basic queries to confirm that it contains exactly 10 structures. This provides a seamless and interactive way for the user to verify the core functionality from start to finish.

### UAT-002: Handle invalid and malformed configuration files gracefully

This scenario is crucial for testing the robustness and user-friendliness of the application's input validation. A tool that crashes with a cryptic traceback on invalid input provides a poor user experience. This test will verify that the application handles common configuration errors in a controlled and informative way. The Jupyter Notebook will have a series of subsections for this scenario, each demonstrating a different type of error. For example, one subsection will attempt to run the pipeline with a path to a non-existent configuration file, and the user will be asked to verify that the application prints a clear "File not found" error. Another subsection will use a configuration file with a missing required field (e.g., `elements`), and the user will verify that the Pydantic validation error is caught and presented in a user-friendly format, such as "Error: Missing required field 'system.elements'". A third example will use a file with the wrong data type (e.g., `n_structures: "ten"` instead of `n_structures: 10`), again verifying a clear error message. This demonstrates to the user that the system is robust against common mistakes and provides helpful feedback.

### UAT-003: Verify the contents and structural integrity of the generated database

This scenario goes a step deeper than UAT-001 by performing more rigorous checks on the output data itself. Its purpose is to confirm that the generated structures are not just present, but are also physically plausible and consistent with the input configuration. After running a successful generation (as in UAT-001), the Jupyter Notebook will provide code cells to perform more detailed analysis of the output database. The first check will be to iterate through every structure in the database and verify its chemical composition, ensuring that the ratio of Cu to Au atoms is correct. The second, more advanced check will be to calculate the minimum interatomic distance for each structure. The test will verify that this distance is above a reasonable physical threshold (e.g., 1.5 Angstroms), confirming that the generation process is not creating unrealistic structures with overlapping atoms. The notebook might also include a simple 3D visualization of one of the generated structures using a library like `nglview`, providing the user with direct visual confirmation of the output's quality. This gives the user a high degree of confidence in the physical validity of the generated data.

## 2. Behavior Definitions

The following Gherkin-style definitions provide a formal, unambiguous description of the expected system behavior for each test scenario. This level of detail is essential for clear communication between developers, testers, and stakeholders.

### Scenario: Successful Generation of Alloy Structures

**GIVEN** a valid YAML configuration file named `config.yml` is present in the current directory.
**AND** the configuration file specifies a `system` with `elements: ['Cu', 'Au']`, a `composition` of `{'Cu': 0.5, 'Au': 0.5}`, a `crystal_structure` of `'fcc'`, and `n_structures: 10`.
**AND** the configuration file specifies a `db` path of `'./cu_au.db'`.

**WHEN** the user executes the command `mlip-autopipec run-pipeline --config config.yml` in the terminal.

**THEN** the application should execute without raising any exceptions and display a success message in the console.
**AND** the application should exit with a status code of 0.
**AND** a new SQLite database file named `cu_au.db` should be created in the current directory.
**AND** when connecting to the `cu_au.db` database, it should report that it contains exactly 10 rows (structures).
**AND** upon inspecting each structure in the database, the chemical formula for each should be 'AuXCuY' where the ratio of X to Y is consistent with the 50/50 composition.
**AND** upon inspecting each structure, the crystal lattice should be consistent with the fcc structure.

### Scenario: Handling a Missing Configuration File

**GIVEN** the current directory does not contain a file named `non_existent_config.yml`.

**WHEN** the user executes the command `mlip-autopipec run-pipeline --config non_existent_config.yml`.

**THEN** the application should not crash or produce a Python traceback.
**AND** the application should print a clear, user-friendly error message to the console, such as "Error: Configuration file not found at path 'non_existent_config.yml'".
**AND** the application should exit with a non-zero status code (e.g., 1), indicating that an error occurred.
**AND** no new database file should be created in the current directory.

### Scenario: Handling a Configuration with a Semantically Invalid Field

**GIVEN** a YAML configuration file named `invalid_config.yml` is present.
**AND** the configuration file is syntactically valid YAML, but contains a semantic error, such as `system.composition` values that do not sum to 1.0 (e.g., `{'Cu': 0.6, 'Au': 0.6}`).

**WHEN** the user executes the command `mlip-autopipec run-pipeline --config invalid_config.yml`.

**THEN** the application should not crash or produce a Python traceback.
**AND** the application should print a clear, user-friendly error message resulting from the Pydantic model validation. The message should precisely identify the source of the error, for example: "Error: Validation error in configuration: The sum of composition values must be equal to 1.0.".
**AND** the application should exit with a non-zero status code.
**AND** no new database file should be created.
