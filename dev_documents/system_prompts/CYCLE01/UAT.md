# User Acceptance Testing (UAT): CYCLE01 - Core Engine & Workflow Foundation

## 1. Test Scenarios

This document outlines the User Acceptance Testing scenarios for CYCLE01. The goal of this cycle is to establish the fundamental workflow of labelling a single atomic structure and training a model from it. The "user" for this UAT is a developer or a power-user who is comfortable with the command line and verifying file outputs.

| Scenario ID | Test Scenario                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C01-01  | Successful execution of the labelling step  | High     |
| UAT-C01-02  | Successful execution of the training step   | High     |
| UAT-C01-03  | Verification of persisted data in database  | Medium   |

---

### **Scenario UAT-C01-01: Successful execution of the labelling step**

**(Min 300 words)**
This scenario is the most critical for CYCLE01 as it validates the core functionality of the `LabelingEngine`. The user, acting as a developer, will test the system's ability to take a predefined atomic structure, correctly interface with the Quantum Espresso (QE) executable, parse the output, and store the results. The acceptance criteria are focused on the successful orchestration of this process. The user will first prepare a simple input: a database file containing a single 'unlabeled' water molecule. They will then execute the command-line interface's `label` command.

The primary goal is to confirm that the system generates a valid QE input file based on the atomic structure and the provided configuration. The user will be expected to inspect this generated input file (`qe_input.in`) to ensure it has the correct atomic positions, cell parameters, and DFT settings (like `ecutwfc`). They will then confirm that the system correctly invokes the `pw.x` command via the command line. For this test, a mock script might be used in place of the actual `pw.x` to verify the command-line arguments passed to it. Finally, the user must verify that the system can parse the resulting output file (`qe_output.out`). They will check the system's logs or screen output to see that the energy, forces, and stress have been successfully extracted and reported. Success in this scenario demonstrates that the fundamental link between the pipeline and the external DFT code is functioning correctly.

---

### **Scenario UAT-C01-02: Successful execution of the training step**

**(Min 300 words)**
This scenario tests the second half of the core workflow: the `TrainingEngine`. It validates that the system can take the data generated from the labelling step and successfully train a Machine Learning Interatomic Potential (MLIP). The user will begin this test with the database state from the successful completion of UAT-C01-01, where at least one structure has been labelled. They will then execute the `train` command from the CLI.

The acceptance criteria for this scenario focus on the data pipeline and the model creation process. The user will first verify, through log outputs, that the `TrainingEngine` correctly queries the database and retrieves the labelled data. The system should report how many structures it has found and is using for training. The main goal is to confirm that an MLIP model file is actually created. The user will navigate to the project's output directory and confirm the existence of a new file, for example, `ace_model.json` or a similar format produced by the `pacemaker` library. While the scientific accuracy of this model is not under scrutiny in this early cycle (as it's trained on very little data), its successful creation is paramount. The user might also inspect the training logs to ensure the training process ran without any numerical or library-specific errors. A successful test here proves that the data flow from the database to the training library is correct and that the system can produce the primary artifact—the MLIP—that is the goal of the entire pipeline.

---

### **Scenario UAT-C01-03: Verification of persisted data in database**

**(Min 300 words)**
This scenario is designed to ensure data integrity and traceability, which are core principles of the system. It tests the `AseDBWrapper`'s ability to correctly manage the state of each atomic structure as it moves through the pipeline. The user will interact directly with the SQLite database file to verify its contents at different stages of the workflow.

The test will proceed in two stages. First, after adding an initial structure but before running the labelling step, the user will open the `asedb.db` file with a tool like `sqlite3` or `ase db`. They will query the database and verify that the row for the new structure exists and its `state` key is correctly set to `'unlabeled'`. After successfully running the `label` command (as in UAT-C01-01), the user will inspect the database again. This time, the acceptance criteria are stricter. They must verify that the `state` key for that row has been updated to `'labeled'`. Furthermore, they must check for the existence of the `dft_result` key. They will be expected to query this key, retrieve its value (which will be a JSON string), and manually inspect it to confirm that it contains the correct, non-zero values for energy, forces, and stress that were parsed from the QE output. This detailed verification ensures that the database is not just a black box but a transparent and reliable record of the entire computational history.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C01-01 - Successful execution of the labelling step**

*   **GIVEN** a clean project environment set up with `uv pip install -e .`.
*   **AND** a valid `asedb.db` file is present in the working directory.
*   **AND** the database contains a single atomic structure (e.g., H2O) with its `state` marked as `'unlabeled'`.
*   **AND** a configuration file is present which specifies the command to run Quantum Espresso.
*   **AND** the specified Quantum Espresso command points to a valid (or mock) executable.
*   **WHEN** the user executes the command `cdd label --id 1` from the command line.
*   **THEN** the system should log a message indicating that it is starting the labelling process for structure ID 1.
*   **AND** the `LabelingEngine` should create a valid Quantum Espresso input file named `qe_input.in` in a temporary directory.
*   **AND** this `qe_input.in` file should contain the correct atomic symbols and positions for the H2O molecule.
*   **AND** the system should execute the Quantum Espresso command as a subprocess.
*   **AND** after the subprocess completes, the system should parse the `qe_output.out` file.
*   **AND** the system should log the extracted values for energy, forces, and stress.
*   **AND** the `state` of the structure with ID 1 in the `asedb.db` should be updated to `'labeled'`.
*   **AND** the `dft_result` key-value pair for structure ID 1 should be populated with a JSON string containing the parsed energy, forces, and stress.

---

**Scenario: UAT-C01-02 - Successful execution of the training step**

*   **GIVEN** the state of the system after the successful completion of scenario UAT-C01-01.
*   **AND** the `asedb.db` contains at least one structure with its `state` marked as `'labeled'`.
*   **AND** the `'labeled'` structure has a valid `dft_result` entry.
*   **AND** a configuration file is present which specifies the MLIP model type as 'ACE'.
*   **WHEN** the user executes the command `cdd train` from the command line.
*   **THEN** the system should log a message indicating that it is starting the training process.
*   **AND** the `TrainingEngine` should query the database and find all structures marked as `'labeled'`.
*   **AND** the system should log the total number of structures it is using for the training set.
*   **AND** the `TrainingEngine` should initiate the ACE model training process using the `pacemaker` library.
*   **AND** the training process should complete without raising any errors.
*   **AND** a new file named `ace_model.json` (or similar) should be created in the output directory.
*   **AND** the system should log a final message indicating that training is complete and the model has been saved.
