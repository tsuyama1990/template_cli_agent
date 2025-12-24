# CYCLE01 User Acceptance Testing (UAT)

## 1. Test Scenarios

This UAT plan focuses on verifying the core functionality of the CYCLE01 deliverable: the ability to process a set of structures, label them using Quantum Espresso, and train a basic MLIP model from the resulting data.

| Scenario ID | Scenario Name                                | Priority |
| :---------- | :------------------------------------------- | :------- |
| UAT-C01-001 | Successful End-to-End Workflow (Happy Path) | High     |
| UAT-C01-002 | Handling of DFT Calculation Failure          | High     |
| UAT-C01-003 | Configuration File Validation                | Medium   |
| UAT-C01-004 | Data Provenance and Database Integrity       | High     |
| UAT-C01-005 | Delta Learning Functionality Verification    | Medium   |

---

**Scenario UAT-C01-001: Successful End-to-End Workflow (Happy Path)**

*   **Description**: This test verifies that the system can successfully execute the entire workflow on a valid, simple input without any errors. It is the primary "happy path" test case that ensures all components are correctly integrated. The test will use a small set of silicon (Si) crystal structures.
*   **Preconditions**:
    *   A working Quantum Espresso (`pw.x`) executable is available in the system's PATH or at a location specified in the config.
    *   `uv` and a valid Python environment are set up.
    *   A directory containing 2-3 simple silicon atomic structures (e.g., different lattice constants) in POSCAR format is prepared.
    *   A valid `input.yaml` is created, pointing to the structures and specifying valid DFT and training parameters.
*   **Acceptance Criteria**:
    *   The CLI command `uv run mlip-pipe run ...` must execute and exit with a status code of 0.
    *   A database file (e.g., `mlip_autoprope.db`) must be created.
    *   The database must contain entries corresponding to each of the input silicon structures.
    *   Each entry in the database must have valid, non-null values for energy, forces, and stress.
    *   A trained MLIP model file (e.g., `model.yace`) must be created in the specified output directory.
    *   The command's log output should indicate successful completion of both the DFT labeling and the model training stages.

---

**Scenario UAT-C01-002: Handling of DFT Calculation Failure**

*   **Description**: This test ensures the system is resilient to errors during the DFT calculation stage. It verifies that the failure of one calculation does not crash the entire workflow and that the system correctly logs the error. The test will include one valid structure and one intentionally invalid structure designed to make QE fail (e.g., atoms too close together).
*   **Preconditions**:
    *   Same as UAT-C01-001.
    *   The input directory of structures contains one valid Si structure and one invalid structure (e.g., two Si atoms placed 0.1 Angstroms apart).
*   **Acceptance Criteria**:
    *   The CLI command must still execute and exit with a status code of 0 (or a specific non-zero code indicating partial success, to be decided).
    *   The system must log a clear error message indicating which specific structure failed the DFT calculation and why (e.g., "SCF did not converge").
    *   The database must contain an entry for the *successful* calculation.
    *   The database must *not* contain a completed entry for the failed calculation (or it may contain a record with a "failed" status).
    *   The `TrainingEngine` must proceed to train a model using only the data from the successfully calculated structure.
    *   A trained MLIP model file must still be created.

---

## 2. Behavior Definitions

**Behavior for UAT-C01-001: Successful End-to-End Workflow**

```gherkin
Feature: Core Workflow Execution
  As a materials scientist,
  I want to run the full data generation and training pipeline on a set of structures,
  So that I can obtain a basic MLIP model.

  Scenario: Running a standard, valid calculation
    GIVEN a directory containing valid silicon crystal structures
    AND a valid `input.yaml` configuration file specifying DFT and ACE training parameters
    AND a functional Quantum Espresso installation
    WHEN I execute the command `uv run mlip-pipe run --structures-dir <path> input.yaml`
    THEN the process should complete successfully with an exit code of 0
    AND a database file named `mlip_autoprope.db` should be created
    AND the database should contain a record for each input structure
    AND each record should include calculated energy, forces, and stress values
    AND a trained model file named `model.yace` should be saved to the output directory.
```

---

**Behavior for UAT-C01-002: Handling of DFT Calculation Failure**

```gherkin
Feature: System Resilience to DFT Failures
  As a user,
  I want the system to gracefully handle failures in individual DFT calculations,
  So that a single bad structure does not halt the entire data generation process.

  Scenario: Processing a dataset with one invalid structure
    GIVEN a directory containing one valid and one intentionally invalid atomic structure
    AND a valid `input.yaml` configuration file
    AND a functional Quantum Espresso installation
    WHEN I execute the command `uv run mlip-pipe run --structures-dir <path> input.yaml`
    THEN the process should complete without crashing
    AND the system's log output should clearly indicate that one structure failed to be calculated
    AND the database should contain a record for the valid structure
    AND the database should not contain a completed record for the invalid structure
    AND the training engine should proceed using only the data from the valid structure
    AND a trained model file should still be created.
```

---

**Behavior for UAT-C01-004: Data Provenance and Database Integrity**

```gherkin
Feature: Data Provenance Tracking
  As a researcher,
  I need to trust that all calculations are tracked and reproducible,
  So that I can verify the integrity of my trained models.

  Scenario: Verifying database records after a successful run
    GIVEN a successful run of the workflow has been completed
    WHEN I inspect the generated database file
    THEN I should find a 'calculations' table
    AND each row in the table should correspond to one successfully processed input structure
    AND each row must contain a unique ID, the input structure details, the full DFT results, and a timestamp.
    AND I should also find a 'models' table
    AND this table should contain a record for the trained MLIP
    AND the model record should contain a link to the calculation IDs used for its training, ensuring traceability.
```
