# User Acceptance Test (UAT): Cycle 1

**Version:** 1.0.0
**Status:** Ready for Testing

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 1. The focus is on verifying that the core automation pipeline for DFT labeling and MLIP training is functional, robust, and produces scientifically valid results. Each scenario is designed to build confidence in a specific aspect of the Cycle 1 functionality, from the high-level workflow to the detailed correctness of the underlying calculations.

| Scenario ID | Test Scenario Description                                                                | Priority |
|-------------|------------------------------------------------------------------------------------------|----------|
| UAT-C1-001  | **Verify Successful End-to-End Run on a Standard System:** The user provides a set of pre-defined structures for crystalline silicon and a minimal `input.yaml`. The system must run from start to finish without errors and produce a valid MLIP file. This is the primary "happy path" test, ensuring all components are correctly integrated and functional. | High     |
| UAT-C1-002  | **Verify Correctness of DFT Labels:** The user provides a single, simple structure (e.g., a water molecule). The system runs the labeling engine, and the user verifies that the calculated energy and atomic forces in the database are consistent with a manually run, reference Quantum Espresso calculation. This test is crucial for ensuring the scientific validity of the generated data. | High     |
| UAT-C1-003  | **Verify "Two-Tier Config" Expansion:** The user provides a minimal `input.yaml` for a magnetic system (e.g., bcc Iron). The user inspects the generated `exec_config_dump.yaml` to confirm that the system correctly enabled spin-polarized DFT settings (`nspin=2`) and set an initial magnetic moment. This validates the intelligence of the heuristic engine. | Medium   |
| UAT-C1-004  | **Verify Handling of Failed DFT Calculation:** The user provides a deliberately malformed or physically impossible input structure that is known to cause SCF convergence failure. The system should gracefully handle the error, log it appropriately, and continue to process other valid structures without crashing. This tests the robustness and fault tolerance of the labeling engine. | Medium   |
| UAT-C1-005  | **Verify Delta Learning Implementation:** The user runs the full pipeline and inspects the trained model. The user verifies that the potential's behavior at very short interatomic distances is repulsive, consistent with the ZBL baseline, confirming that the delta learning is active. This ensures the physical plausibility of the trained potential. | High     |


## 2. Behavior Definitions

This section provides detailed Gherkin-style behavior definitions for the test scenarios, outlining the exact steps and expected outcomes for a user performing the acceptance test.

---

### **Scenario: UAT-C1-001 - Successful End-to-End Run on a Standard System**

*   **GIVEN** a clean working directory containing two files:
    *   A minimal configuration file named `input.yaml`, with content specifying the elements involved, for example: `composition: "Si"` and `elements: ["Si"]`.
    *   An ASE trajectory file named `structures.traj`, which contains a representative set of 10 to 20 varied atomic configurations of crystalline silicon. This set should include the equilibrium structure, structures with applied volumetric and shear strains, and structures with random atomic displacements ("rattled" structures).
*   **WHEN** the user executes the pipeline from the command line, for instance, by running the command `mlip-pipe run-cycle1 input.yaml structures.traj`.
*   **THEN** the system should begin execution, providing clear log messages indicating the start of the `ConfigExpander`, `LabelingEngine`, and `TrainingEngine` stages.
*   **AND** the entire workflow must complete without any crashes, Python exceptions, or fatal error messages being displayed to the user.
*   **AND** upon completion, a new file named `exec_config_dump.yaml` must exist in the directory. A user inspection of this file should reveal a comprehensive set of DFT and training parameters, demonstrating that the minimal input was successfully expanded.
*   **AND** a database file, typically `results.db`, must be created. When queried (e.g., using `ase db results.db -c`), this database must show a count of entries that is exactly equal to the number of structures that were in the input `structures.traj` file.
*   **AND** a trained model file, for example `final_model.ace`, must be present in the output directory. This file should be non-empty and have a file size consistent with a trained potential.
*   **AND** the total time taken for the execution should be reasonable and commensurate with the workload (e.g., for 20 silicon structures, it might take 10-30 minutes, depending on the machine).

---

### **Scenario: UAT-C1-002 - Verify Correctness of DFT Labels**

*   **GIVEN** a directory containing:
    *   A minimal `input.yaml`.
    *   An `structures.traj` file that contains only a single, precisely defined structure of a water molecule in a vacuum.
*   **AND** the user has, independently of the pipeline, run a manual Quantum Espresso calculation on an identical water molecule structure, using the exact same pseudopotentials and cutoff energies that the pipeline is configured to use. The user has noted the resulting total energy (`E_ref`) and the forces on each of the three atoms (`F_ref`).
*   **WHEN** the user runs only the `LabelingEngine` part of the pipeline on the provided single-structure file.
*   **THEN** the system should create or update the results database file.
*   **AND** when the user queries the database for the results of the calculation on the water molecule, the retrieved energy value (`E_db`) must be in very close agreement with the reference energy `E_ref`, for example, `abs(E_db - E_ref) < 1e-5` eV.
*   **AND** the retrieved forces for each atom (`F_db`) must be in very close agreement with the reference forces `F_ref`, with the root-mean-square difference between the force vectors being very small, for example, `sqrt(mean((F_db - F_ref)^2)) < 1e-4` eV/Angstrom.

---

### **Scenario: UAT-C1-003 - Verify "Two-Tier Config" Expansion**

*   **GIVEN** a directory containing a minimal `input.yaml` with the content `composition: "Fe"` and `elements: ["Fe"]`.
*   **WHEN** the user runs the `ConfigExpander` component of the pipeline, which is the first step of any run.
*   **THEN** a file named `exec_config_dump.yaml` must be generated.
*   **AND** when the user opens this file for inspection, under the main `dft_compute` section, there must be a parameter `magnetism` set to a physically correct value such as `"ferromagnetic"`.
*   **AND** within the subsection for the Quantum Espresso parameters, the configuration must explicitly include `nspin = 2` to enable spin polarization.
*   **AND** the configuration for the Iron atom type must include a non-zero starting magnetization, for example, `starting_magnetization(1) = 0.5`. This confirms the heuristic engine correctly identified a magnetic element and applied the appropriate, more complex DFT settings.

---

### **Scenario: UAT-C1-004 - Verify Handling of Failed DFT Calculation**

*   **GIVEN** a directory containing an `structures.traj` file with five atomic structures. Four of these are physically reasonable, but one is deliberately malformed, containing two atoms placed at an impossibly close distance (e.g., 0.1 Angstroms apart), which is guaranteed to make the DFT calculation fail.
*   **WHEN** the user runs the `LabelingEngine` part of the pipeline on this mixed-quality dataset.
*   **THEN** the system's log output must contain clear warning or error messages indicating that the calculation for one of the structures failed to complete successfully. The message should ideally mention the likely cause, such as "SCF failed to converge."
*   **AND** the pipeline must not crash or terminate prematurely. It should continue processing the remaining, valid structures.
*   **AND** upon completion, the final results database must contain exactly four entries, corresponding to the four valid input structures.
*   **AND** there must be no entry in the database associated with the malformed structure. This demonstrates the system's ability to isolate failures and maintain data integrity.

---

### **Scenario: UAT-C1-005 - Verify Delta Learning Implementation**

*   **GIVEN** a valid trained MLIP file (e.g., `final_model.ace`) that has been successfully produced by the full Cycle 1 pipeline.
*   **WHEN** the user executes a small, purpose-written script (which could be provided as part of the examples) that loads this potential file.
*   **AND** the script then uses the loaded potential to calculate the interaction energy of a dimer of the relevant atoms (e.g., two Silicon atoms) as a function of their separation distance, particularly focusing on very short distances (e.g., from 1.0 down to 0.5 Angstroms).
*   **THEN** the calculated energy must be a large, positive value, indicating a strong repulsive force, as expected from physical principles.
*   **AND** as the script decreases the separation distance, the energy must increase steeply and monotonically, following the repulsive wall of the baseline pair potential (like ZBL) used in the delta learning. This behavior confirms that the model has not learned some unphysical attractive potential at short range and that the delta learning framework is correctly ensuring physically sound behavior in this regime.
