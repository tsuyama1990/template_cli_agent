# Cycle 01 User Acceptance Test (UAT): Core Engine

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 01. The focus is on verifying the fundamental capabilities of the core engine from a user's perspective. This means ensuring its ability to correctly process atomic structures, execute DFT calculations, and train a basic MLIP model in a way that is both correct and transparent. The scenarios are designed to be run from the command line by a user or an automated script, validating the end-to-end workflow of this initial cycle. Each scenario tests a critical aspect of the core functionality, from the successful processing of a single data point to the handling of inevitable calculation failures and the verification of the final trained model's fidelity. Passing these tests will provide high confidence that the foundational components of the pipeline are robust and reliable, ready for the addition of more complex features in subsequent cycles.

| Scenario ID | Description                                                              | Priority |
|-------------|--------------------------------------------------------------------------|----------|
| UAT-01-001  | **Successful Single-Point Calculation:** Verify that the system can take a single, valid atomic structure, perform a DFT calculation using automated parameters, and correctly store the resulting energy, forces, and stress in the database. This test is the most fundamental unit of work for the entire system and validates the core data generation capability. | High     |
| UAT-01-002  | **Successful Training from Labelled Data:** Verify that the system can use a pre-populated database of labelled structures to successfully train an MLIP and save the resulting model file. This test confirms that the training engine can correctly consume the data produced by the labelling engine and that the underlying machine learning framework is correctly integrated. | High     |
| UAT-01-003  | **Full End-to-End Workflow for a Simple System:** Verify the complete linear workflow for a simple, well-known system (e.g., an Argon dimer). This involves starting with raw structure files, running the labelling process, and then automatically proceeding to the training process, resulting in a final model. This test demonstrates that the orchestrator can correctly chain the modules together. | High     |
| UAT-01-004  | **Handling of Failed DFT Calculation:** Verify that if a Quantum Espresso calculation fails (e.g., due to SCF non-convergence), the system correctly identifies the failure, logs an appropriate and informative error message, and does not incorrectly pollute the database with erroneous or incomplete data. This tests the system's robustness and error-handling capabilities. | Medium   |
| UAT-01-005  | **Verification of Delta Learning:** Verify that the training process correctly applies the Delta Learning methodology. This will be checked by writing a post-processing script to confirm that the trained model, when its predictions are added to the baseline potential's contribution, can accurately reproduce the original DFT energies, forces, and stresses, demonstrating the fidelity of the learned model. | Medium   |

## 2. Behaviour Definitions

The following Gherkin-style definitions describe the expected behaviour for each test scenario in detail.

---

### **Scenario: UAT-01-001 - Successful Single-Point Calculation**

This scenario tests the most fundamental piece of the workflow: the ability to correctly label a single atomic structure with its DFT properties. It validates the input generation, process execution, and output parsing capabilities of the Labelling Engine. A standard, simple crystalline material like bulk silicon is chosen because its properties are well-known and calculations are expected to converge easily. The success of this test is a prerequisite for any further functionality, as it confirms that the pipeline can generate valid training data.

**GIVEN** a clean and empty ASE database file has been initialised at a known location.
**AND** a user has a single CIF file representing a standard 8-atom conventional cell of bulk silicon.
**AND** a script is used to load this structure into the database, resulting in a single entry with the state "unlabelled".
**WHEN** the user executes the `LabellingEngine`'s main script, providing it with the path to the database and the relevant configuration.
**THEN** the system should log an informational message indicating it has found one unlabelled structure and is beginning processing.
**AND** the system should generate a valid Quantum Espresso input file for a static (SCF) single-point calculation in a temporary directory. The input file should contain physically sensible parameters for silicon, automatically determined by the system's SSSP-based heuristics (e.g., correct pseudopotentials, appropriate `ecutwfc`).
**AND** the system should execute the `pw.x` command as a subprocess, capturing its standard output.
**AND** upon completion of the subprocess, the system should parse the captured output, specifically looking for the final energy, forces, and stress tensor blocks.
**AND** the corresponding entry for the silicon structure in the ASE database should be updated to have the state "labelled".
**AND** the database entry's key-value pairs should be populated with the calculated floating-point value for energy, and NumPy arrays for forces and the stress tensor, consistent with the values present in the Quantum Espresso output file.

---

### **Scenario: UAT-01-002 - Successful Training from Labelled Data**

This scenario tests the functionality of the Training Engine in isolation. It verifies that the engine can correctly query the database for completed calculations, process the data, and successfully execute a training routine to produce a valid MLIP file. This test is crucial for ensuring the machine learning component of the pipeline is correctly integrated and that the data format produced by the Labelling Engine is compatible with the Training Engine's requirements.

**GIVEN** an ASE database file that has been pre-populated with at least 10 "labelled" structures of bulk silicon. These structures could represent the crystal with different amounts of strain applied.
**AND** each of these 10 entries contains valid, non-null data for energy, forces, and stress.
**WHEN** the user executes the `TrainingEngine`'s main script, providing it with the path to the database and the relevant configuration.
**THEN** the system should log a message indicating that it has found 10 labelled structures and is preparing the dataset for training.
**AND** the system should perform the Delta Learning calculations internally, subtracting the contribution of the baseline potential from the DFT values.
**AND** the system should successfully execute the ACE model training process without raising any errors or exceptions.
**AND** a new file representing the trained MLIP (e.g., `model.ace`) should be created in the specified output directory.
**AND** this file should be non-empty and have a file size consistent with a trained ACE model.

---

### **Scenario: UAT-01-003 - Full End-to-End Workflow for a Simple System**

This scenario tests the integration of the two core modules by the orchestrator script. It validates the entire linear workflow for Cycle 01, from raw input files to a final trained model, for a system where the physical behaviour is trivial and well-understood. An Argon dimer is an ideal candidate, as its interaction is dominated by van der Waals forces, which can be captured by even a simple MLIP.

**GIVEN** an empty output directory and a clean ASE database.
**AND** a directory containing two POSCAR files for an Argon dimer: one at its approximate equilibrium bond distance (e.g., 3.8 Å) and one significantly compressed (e.g., 2.5 Å).
**WHEN** the user executes the main orchestrator script for Cycle 01, pointing it to the directory of POSCAR files.
**THEN** the orchestrator should first log that it is populating the database and create two entries with the state "unlabelled".
**AND** the `LabellingEngine` should then be automatically invoked. Its log messages should indicate the successful processing of both structures, which are then updated to the "labelled" state in the database.
**AND** after the labelling is complete, the `TrainingEngine` should be automatically invoked.
**AND** the training process should complete successfully.
**AND** a final MLIP file (e.g., `model.ace`) should be created in the output directory.
**AND** the main console output should clearly indicate the successful completion of both the labelling and training stages, and the script should exit with a status code of 0.

---

### **Scenario: UAT-01-004 - Handling of Failed DFT Calculation**

This scenario tests the robustness and error-handling capability of the Labelling Engine. Real-world DFT calculations can often fail for a variety of numerical reasons, and the system must be able to handle these failures gracefully without crashing or corrupting the dataset.

**GIVEN** a clean ASE database file.
**AND** an `ase.Atoms` object for a structure that is known to be difficult for SCF calculations to converge (e.g., a metal cluster with a challenging electronic structure) is added to the database as "unlabelled".
**WHEN** the user executes the `LabellingEngine`'s main script.
**THEN** the system should attempt to execute the `pw.x` command as normal.
**AND** the system's output parser should detect from the captured standard output that the calculation did not converge successfully (e.g., by finding the "convergence NOT achieved" message from Quantum Espresso).
**AND** the system should write a clear and specific error message to the log file or console (e.g., "ERROR: DFT calculation for database ID 1 failed to converge after 100 iterations.").
**AND** the state of the entry in the ASE database should either remain "unlabelled" or, preferably, be moved to a distinct "error" state to prevent repeated attempts.
**AND** the energy, forces, and stress fields for that entry in the database must not be populated with any partial or incorrect results from the failed run.

---

### **Scenario: UAT-01-005 - Verification of Delta Learning**

This scenario provides a quantitative check on the fidelity of the trained model, specifically validating that the Delta Learning approach has been implemented correctly. It ensures that the final model, when combined with the baseline potential it was trained against, accurately reproduces the ground-truth DFT data.

**GIVEN** a trained MLIP model that was produced by the `TrainingEngine` for a simple system (e.g., the Argon dimer from UAT-01-003).
**AND** the original DFT-labelled data for that same system, which is still present in the ASE database.
**WHEN** a user writes and executes a separate analysis script (or a dedicated subcommand of the main application) that loads the trained MLIP and an implementation of the baseline potential used during training.
**AND** the script iterates through the original dataset, and for each structure, it calculates the predicted total energy as `E_predicted = E_MLIP + E_baseline`.
**THEN** the array of `E_predicted` values should closely match the original `E_DFT` values stored in the database.
**AND** the Root Mean Square Error (RMSE) between the predicted and DFT energies should be below a predefined, small threshold (e.g., < 0.01 eV/atom), which quantitatively demonstrates that the model has learned the training data to a high degree of accuracy.
**AND** a similar comparison for the forces should also show a low RMSE (e.g., < 0.05 eV/Å).