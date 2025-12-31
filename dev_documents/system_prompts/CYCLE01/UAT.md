# UAT.md: Cycle 01 - The Core Engine

## 1. Test Scenarios

User Acceptance Testing (UAT) for the foundational Cycle 01 is designed to be executed from the perspective of a technical user, such as a developer, a system administrator, or a power-user scientist who is comfortable with the command line and basic scripting. Since there is no formal user interface in this cycle, the UAT process is centered on running a prescribed Python script that drives the core engine and then meticulously verifying that the resulting outputs and data artifacts are correct, complete, and stored as expected. The overarching goal is to build deep confidence in the fundamental reliability of the data processing chain. Success in this UAT means that the system can be trusted to take a single, well-defined atomic structure, process it through the DFT and MLIP training engines without error, and correctly record the entire transaction in the database. This validation is critical, as all future, more complex functionalities will depend entirely on the robustness of this core workflow.

| Scenario ID | Scenario Description                                       | Priority |
| :---------- | :--------------------------------------------------------- | :------- |
| UAT-C01-01  | **Successful End-to-End Workflow Verification ("Happy Path")** | High     |
|             | This is the most critical scenario for the cycle, testing the ideal "happy path." A technical user provides a simple, physically valid, and computationally inexpensive atomic structure (like a H2 molecule). The expectation is that the system will flawlessly execute the entire end-to-end workflow without any manual intervention. This includes successfully launching and completing a Quantum Espresso calculation, correctly parsing the results, storing these results in the database, and then successfully training a basic MLIP model from this single data point. A successful outcome for this scenario provides the definitive demonstration that the core components (Labelling Engine, Training Engine, Database) are not only functional in isolation but are also correctly integrated and communicating with each other as designed. It is the primary proof of concept for the entire project's architecture. |          |
| UAT-C01-02  | **Graceful Handling of a Predictable DFT Calculation Failure**           | High     |
|             | This scenario is designed to test the system's robustness and resilience in the face of a common real-world problem. The user provides a deliberately malformed or physically unstable atomic structure that is known to cause convergence issues in Quantum Espresso (e.g., two atoms placed unphysically close). The absolute requirement is that the system must not crash, hang, or enter an infinite loop. Instead, the user expects the Labelling Engine to correctly identify the failure from the DFT output, log the error in a detailed and informative way, cleanly terminate the processing for that structure, and correctly mark the structure as 'failed' in the database. This test validates the error-handling and state-management logic, which is crucial for building a system that can run long, unsupervised workflows without being derailed by a single problematic data point. |          |
| UAT-C01-03  | **Verification of Data Provenance and Database Integrity**                 | Medium   |
|             | This scenario focuses on the "single source of truth": the database. The user's goal is to verify that the system is meticulously and accurately recording its activities, ensuring data provenance and traceability. After running a successful workflow (as in UAT-C01-01), the user will use a standard SQLite client to manually inspect the contents of the generated `mlip.db` file. They will perform checks to ensure that the original atomic structure, the complete set of calculated DFT results (energy, forces, stress), and the final status ('labelled') are all correctly stored and logically linked together. This test is vital for ensuring that the data produced by the pipeline is reliable, auditable, and can be trusted for scientific analysis. It confirms that the database is not just a data dump, but a structured and trustworthy record of the entire computational process. |          |

## 2. Behavior Definitions

The behaviors for Cycle 01 UAT will be verified by a technical user who will execute a Python test script from the command line. The verification process will involve observing the script's console output for success or failure messages, checking for the existence and content of output files, and using an external tool (a SQLite database browser) to inspect the state of the database after the script has run.

### Scenario: UAT-C01-01 - Successful End-to-End Workflow Verification ("Happy Path")

*   **GIVEN** a completely clean and empty project environment, specifically ensuring that no `mlip.db` file or trained model output files exist.
*   **AND** the user has access to a Python script (`run_cycle01_test.py`) that imports the orchestrator, defines a simple, physically valid `ase.Atoms` object (e.g., an H2 molecule with a bond length near its equilibrium of ~0.74 Angstroms), and passes this object to the main workflow function.
*   **AND** a working installation of Quantum Espresso, specifically the `pw.x` executable, is correctly configured and accessible in the system's PATH.
*   **WHEN** the user executes the test script from their terminal (e.g., `python run_cycle01_test.py`).
*   **THEN** the script must execute to completion without raising any unhandled Python exceptions and must terminate with a success status code of 0.
*   **AND** the user should observe clear, informative messages being printed to the console log, indicating the successful completion of each major stage: "DFT calculation successful," "Results written to database," and "Model training complete."
*   **AND** upon completion, a new file representing the trained MLIP model (e.g., `model_cycle01.pt`) must be created in the designated output directory.
*   **AND** a final, detailed inspection of the SQLite database must confirm that it contains a single record for the H2 molecule and that this record's `state` column is correctly and unambiguously marked as 'labelled'. This confirms that the state transition logic is working correctly for the success case.

### Scenario: UAT-C01-02 - Graceful Handling of a Predictable DFT Calculation Failure

*   **GIVEN** a completely clean and empty project environment.
*   **AND** the user has access to a Python script (`run_cycle01_failure_test.py`) that defines an `ase.Atoms` object that is deliberately and predictably problematic for SCF convergence (e.g., two Argon atoms placed at an extremely short distance of 0.5 Angstroms, which would cause a massive, unphysical repulsive force).
*   **AND** a working installation of Quantum Espresso is available.
*   **WHEN** the user executes this specific failure-testing script from their terminal.
*   **THEN** the script must execute to completion without crashing, hanging, or entering any kind of infinite retry loop. It should terminate cleanly.
*   **AND** the user should observe clear and explicit error messages in the console log. These messages should unambiguously state that the DFT calculation failed and should, critically, include the key error messages captured from the Quantum Espresso output (e.g., "Error: SCF not converged in 100 steps").
*   **AND** the user must verify that no trained MLIP model file has been created, as the training step should have been skipped due to the lack of valid input data.
*   **AND** a final, detailed inspection of the SQLite database must show that it contains a single record for the problematic Argon dimer structure, but its `state` column must be correctly marked as 'failed' or 'error'. This confirms that the system can correctly identify, flag, and isolate problematic data points without halting the entire workflow.

### Scenario: UAT-C01-03 - Verification of Data Provenance and Database Integrity

*   **GIVEN** that the successful workflow from scenario UAT-C01-01 has been executed and has completed without errors.
*   **WHEN** the user opens the generated `mlip.db` file using a standard, external SQLite database client (such as `sqlite3` on the command line or a graphical tool like DB Browser for SQLite).
*   **THEN** the user must be able to execute a SQL query (e.g., `SELECT * FROM structures;`) and verify that the table contains exactly one entry, corresponding to the H2 molecule.
*   **AND** this entry must have its 'state' field explicitly set to the string 'labelled'.
*   **AND** the database must also contain the complete, associated DFT results. The user will verify the presence and correctness of: a single floating-point value for the energy, a NumPy array of forces with the correct dimensions (2 atoms x 3 dimensions), and a NumPy array for the stress tensor.
*   **AND** the user will perform a "sanity check" on the physical values. The stored energy value must be physically reasonable for a stable H2 molecule (i.e., a significant negative value, indicating a bound state relative to two isolated H atoms).
*   **AND** the stored forces on the two atoms should be very close to zero (e.g., < 1e-4 eV/Angstrom), as the input structure was intentionally created near its equilibrium bond length. This confirms that the DFT calculation and the data parsing were not just syntactically correct, but also physically meaningful.
