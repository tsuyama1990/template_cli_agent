# User Acceptance Testing (UAT): Cycle 1

This document outlines the User Acceptance Testing (UAT) scenarios for the first cycle of the MLIP-AutoPipe project. The focus of this cycle is the core command-line functionality for generating a database of initial atomic structures. These tests are designed from the perspective of a user whose primary goal is to quickly and reliably generate a set of starting configurations for a materials science project. The successful execution of these scenarios will confirm that the application is not only functional but also robust, user-friendly, and provides a tangible benefit for its target audience.

## 1. Test Scenarios

The following scenarios are designed to be executed by a user from the command line to verify that the core requirements of Cycle 1 have been met. They cover the primary success path, common user errors, and important quality-of-life features like idempotency.

### Scenario ID: UAT-C1-001
-   **Priority**: High
-   **Title**: Successful Generation of a Binary Alloy Database
-   **Description**: This is the primary success-path scenario and represents the most common use case for the functionality delivered in Cycle 1. The user will define a simple Iron-Platinum (Fe-Pt) binary alloy in a YAML configuration file, specifying the composition, the number of structures to generate, and the desired size of the simulation cell (supercell). They will then execute the main command-line interface command, pointing to this configuration file. The successful completion of this test is a critical milestone, as it validates the entire end-to-end pipeline for the core workflow. It provides the user with immediate confidence that the fundamental components of the system—configuration parsing by Hydra, Pydantic model validation, the `AlloyGenerator`'s structure creation logic, the physical validation checks (e.g., for atomic overlap), and the final storage of data into the ASE database—are all working together correctly. A successful run will produce a tangible and useful output: a `FePt_structures.db` file. This file is not a black box; it is a standard ASE database that can be immediately used with other tools in the materials science ecosystem for visualisation or further analysis, demonstrating the immediate practical value of the application even at this early stage. The user should be able to inspect this database and confirm that the number of structures and their chemical compositions match the input request, providing a clear and unambiguous definition of success.

### Scenario ID: UAT-C1-002
-   **Priority**: High
-   **Title**: Handling of Invalid Configuration - Composition Error
-   **Description**: This scenario tests the system's robustness and its ability to provide a positive user experience in the face of incorrect input, which is a critical aspect of usable software. The user will intentionally create a syntactically valid YAML configuration file but one that is semantically incorrect according to the physical rules of the system. Specifically, the chemical compositions for the alloy components will not sum to 1.0 (e.g., Fe: 0.6, Pt: 0.5, for a total of 1.1). When the user attempts to run the CLI with this invalid configuration, the system must *not* crash with an unhandled exception or a cryptic traceback. Instead, it should fail gracefully and provide a clear, informative error message directly to the terminal. The error message should explicitly state what the problem is (e.g., "Validation Error: Compositions must sum to 1.0") and ideally should point to the specific field in the configuration file that is incorrect. This test is crucial for verifying that the Pydantic-based validation layer is functioning as designed. It demonstrates that the application is actively guiding the user towards correcting their input, which saves time and frustration and is a hallmark of a well-designed tool. Success in this scenario means the user is empowered to fix their own mistakes without needing to consult documentation or, worse, read the source code.

### Scenario ID: UAT-C1-003
-   **Priority**: Medium
-   **Title**: Idempotent Execution (Pipeline Resume)
-   **Description**: This scenario verifies the pipeline's ability to avoid re-doing completed work, a feature often referred to as idempotency or checkpointing. The user will first execute the successful generation scenario (UAT-C1-001). This will result in the creation of an intermediate `initial_structures.xyz` file and the final `FePt_structures.db` database. The user will then immediately run the exact same command a second time, without cleaning the directory. The expected and correct behavior is that the tool will detect the existence of the intermediate file, recognise that the 'Generation' stage has already been successfully completed for this configuration, and skip it. The application should print a clear message to the log, such as "Found existing initial structures file, skipping Generation stage." This demonstrates that the state isolation and checkpointing mechanisms are effective. This is a crucial quality-of-life feature for a pipeline that will eventually perform computationally expensive, long-running simulations in Cycle 2. It gives the user confidence that they can safely re-run scripts without fear of corrupting existing data or unnecessarily wasting valuable computational time and resources. It also provides a simple mechanism for resuming a partially failed workflow, making the overall system more robust and practical for real-world use.

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behavior for each UAT scenario in a structured, unambiguous format.

---

**Scenario**: `UAT-C1-001: Successful Generation of a Binary Alloy Database`

*   **GIVEN** I am in a clean directory that contains no existing database files or `.xyz` files.
*   **AND** I have a working installation of the `mlip-autopipec` application.
*   **AND** I have used a text editor to create a configuration file named `config.yaml` with the following valid and physically plausible content:
    ```yaml
    system:
      elements: [Fe, Pt]
      composition: {Fe: 0.5, Pt: 0.5}
      lattice_constant: 3.8
      lattice: 'fcc' # Face-centered cubic
    generation:
      num_structures: 10
      supercell_size: [3, 3, 3]
    database:
      path: FePt_structures.db
      overwrite: false
    ```
*   **WHEN** I open my terminal, navigate to the directory containing `config.yaml`, and execute the command `mlip-autopipec run --config-path config.yaml`.
*   **THEN** the process should start executing and display informative log messages, indicating that it is proceeding through the "Generation", "Exploration", "Sampling", and "Storage" stages.
*   **AND** the process should complete successfully with an exit code of 0, indicating no errors occurred.
*   **AND** upon completion, a new file named `FePt_structures.db` and an intermediate file `initial_structures.xyz` should be present in the directory.
*   **AND** when I inspect the `FePt_structures.db` file using an external tool like ASE's database browser (`ase db FePt_structures.db`), it must report that it contains exactly 10 rows (atomic structures).
*   **AND** when I inspect any single structure from the database, it must contain a total of 108 atoms (a 3x3x3 supercell of a 4-atom FCC unit cell), comprising exactly 54 Iron (Fe) atoms and 54 Platinum (Pt) atoms, confirming the requested composition was honored.

---

**Scenario**: `UAT-C1-002: Handling of Invalid Configuration - Composition Error`

*   **GIVEN** I am in a clean directory.
*   **AND** I have created a configuration file named `invalid_config.yaml` where the composition sum is intentionally incorrect (1.1 instead of 1.0):
    ```yaml
    system:
      elements: [Fe, Pt]
      composition: {Fe: 0.6, Pt: 0.5} # Sum is 1.1
      lattice_constant: 3.8
      lattice: 'fcc'
    generation:
      num_structures: 10
    database:
      path: test.db
    ```
*   **WHEN** I execute the command `mlip-autopipec run --config-path invalid_config.yaml` in my terminal.
*   **THEN** the process should fail quickly and return a non-zero exit code, signaling an error to any calling scripts.
*   **AND** the terminal must display a user-friendly error message, clearly indicating that a "Pydantic Validation Error" has occurred.
*   **AND** the error message must contain text that explicitly identifies the problem, for example, "The sum of composition values must be 1.0".
*   **AND** no files named `test.db` or `initial_structures.xyz` should be created in the directory, confirming that the pipeline halted before performing any actions.

---

**Scenario**: `UAT-C1-003: Idempotent Execution (Pipeline Resume)`

*   **GIVEN** I have already successfully completed the entire scenario `UAT-C1-001`.
*   **AND** the files `FePt_structures.db` and `initial_structures.xyz` already exist in my current directory.
*   **AND** I have recorded the last modification timestamp of the `FePt_structures.db` file.
*   **WHEN** I execute the exact same command `mlip-autopipec run --config-path config.yaml` for a second time in the same directory.
*   **THEN** the process should complete much more quickly than the first run, ideally in less than a second.
*   **AND** the process should exit with a success code of 0.
*   **AND** the terminal output must contain a clear log message indicating that a stage was skipped, such as "Found existing initial structures file 'initial_structures.xyz', skipping Generation stage."
*   **AND** the last modification timestamp of the `FePt_structures.db` file should be unchanged from what I recorded previously, proving that no new data was written to it.
