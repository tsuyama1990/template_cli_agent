# Cycle 1 User Acceptance Tests: Core Engine and Automation Foundation

**Version:** 1.0.0
**Status:** Final

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 1 of the MLIP-AutoPipe project. The purpose of UAT is to validate the system's functionality from the perspective of an end-user, in this case, a computational materials scientist. The tests are designed not just to check for bugs, but to ensure that the software is genuinely useful, intuitive, and reliable for its intended purpose. For this foundational cycle, the scenarios focus on verifying the core promise of the `LabellingEngine` and `TrainingEngine`: to automate the tedious and expert-driven tasks of DFT data generation and subsequent MLIP training. We need to build confidence that the system can be trusted to perform these critical tasks correctly and robustly. Success in these UAT scenarios is a prerequisite for building more complex features in later cycles, as all future developments will depend on the correctness and reliability of this core engine.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C1-001** | **Successful End-to-End Run on a Standard System:** This is the most fundamental "happy path" test. It verifies that a user can take a simple, well-understood material, provide a few common structures, and have the pipeline execute from start to finish without error, producing the expected output files. This test confirms that all the basic components are correctly integrated and that the system is fundamentally functional. Its successful execution is the primary indicator of the cycle's success. | **High** |
| **UAT-C1-002** | **Correctness of DFT Parameter Automation:** This scenario tests the "expert in a box" feature of the `LabellingEngine`. A key value proposition of the software is that it removes the need for the user to manually look up and specify dozens of complex DFT parameters. This test ensures that the automation is not just present, but correct. It verifies that the system embeds the best practices of the scientific community (specifically the SSSP protocol) for selecting high-precision, reliable DFT settings, thereby ensuring the quality and reproducibility of the generated data. | **High** |
| **UAT-C1-003** | **Robustness to DFT Calculation Failure:** Real-world DFT calculations often fail for a variety of reasons, with electronic convergence failure (SCF failure) being one of the most common. A naive pipeline would simply crash, forcing the user to manually debug the input. This test ensures that our system is more robust. It validates the automated error recovery mechanism, confirming that the system can diagnose a common failure mode and intelligently attempt to fix it, for example, by adjusting the electronic mixing parameters and retrying the calculation. This demonstrates the system's resilience and its ability to handle non-ideal inputs gracefully. | **Medium** |
| **UAT-C1-004** | **Verification of Delta Learning:** The delta learning approach is a critical design choice for ensuring the physical realism and stability of the final MLIP. This test is designed to confirm that this feature is implemented correctly. It's not enough that the training process completes; the user must be confident that the resulting model is a true representation of the underlying DFT data. This test verifies that the combination of the simple baseline potential and the machine-learned delta accurately reproduces the original, high-fidelity DFT calculations, confirming the mathematical correctness of the training engine. | **High** |
| **UAT-C1-005** | **Data Provenance and Persistence:** Scientific reproducibility is paramount. A user must be able to trace back every single piece of data to the exact conditions under which it was generated. This scenario tests the system's data provenance capabilities. It ensures that the ASE database, which stores the results, also stores a rich set of metadata alongside each calculation. This metadata provides an immutable audit trail, allowing the user (or a reviewer) to later inspect a data point and know precisely which version of the DFT code, which pseudopotentials, and which convergence parameters were used to create it. This feature is essential for debugging, publication, and building trust in the automated system. | **Medium** |

---

## 2. Behavior Definitions

These behaviors are described in the Gherkin style (GIVEN/WHEN/THEN) to provide clear, unambiguous, and detailed test cases from the user's point of view.

### **UAT-C1-001: Successful End-to-End Run on a Standard System**

**Scenario:** A user wants to generate a basic MLIP for bulk Silicon, a standard benchmark material in computational science. They have a couple of structure files representing the equilibrium and a slightly strained crystal, and they want the pipeline to handle everything else.
> **GIVEN** a freshly created and activated Python virtual environment.
> **AND** the MLIP-AutoPipe package has been successfully installed via pip into this environment.
> **AND** a valid Quantum Espresso `pw.x` executable is available on the system's PATH, and the user has confirmed it is runnable.
> **AND** the user has created a simple configuration file that, at a minimum, specifies the command to run `pw.x` (e.g., `mpirun -np 4 pw.x`) and the desired name for the output model.
> **AND** the user has created an input directory named `si_structures` containing two valid POSCAR files: `si_eq.poscar` (the equilibrium diamond-structure silicon cell) and `si_strained.poscar` (the same cell with a 1% hydrostatic strain applied).
>
> **WHEN** the user executes the main pipeline command from their terminal, pointing it to the `si_structures` directory and their configuration file.
>
> **THEN** the system should start execution immediately without raising any configuration or input errors.
> **AND** the console output should clearly log the start of the process and indicate that it has found and is processing the two provided structure files.
> **AND** the system should log messages indicating that it is generating Quantum Espresso inputs and launching the `pw.x` subprocess for each structure.
> **AND** after the DFT calculations are complete, the log should indicate that the data is being saved to a database and that the training process is starting.
> **AND** upon completion, the system must exit with a status code of 0, indicating success.
> **AND** an ASE database file, likely named `results.db`, must exist in the specified output directory and its size should be non-zero.
> **AND** a trained model file, for example `silicon.model`, must also exist in the output directory and have a non-zero file size.

### **UAT-C1-002: Correctness of DFT Parameter Automation**

**Scenario:** A user is working with a more complex, multi-element material, Gallium Arsenide (GaAs), and wants to trust the system to select the correct, high-precision DFT settings automatically, as recommended by the materials science community.
> **GIVEN** the MLIP-AutoPipe system is installed and configured.
> **AND** the user creates an input directory containing a single structure file, `GaAs.cif`, representing the zincblende crystal structure of Gallium Arsenide.
>
> **WHEN** the user initiates the labelling pipeline for this single structure.
>
> **THEN** the system must create a working directory for the calculation, and within it, a Quantum Espresso input file (e.g., `gaas.in`).
> **AND** upon inspection of this `gaas.in` file, the `ATOMIC_SPECIES` card must specify the exact, correct SSSP Precision pseudopotential filenames for Gallium (`Ga.us.z_5.ld1.psl.v1.0.0-high.upf`) and Arsenic (`As.us.z_5.ld1.psl.v1.0.0-high.upf`), assuming version 1.0.0 of the SSSP library is being used.
> **AND** the `SYSTEM` namelist within the input file must specify an energy cutoff for the wavefunctions (`ecutwfc`) of at least 60 Ry and an energy cutoff for the charge density (`ecutrho`) of at least 413 Ry, as these are the maximum recommended values for the chosen Ga and As pseudopotentials in the SSSP Precision protocol. This demonstrates the system correctly applies the "use the max cutoff of all elements" rule.

### **UAT-C1-003: Robustness to DFT Calculation Failure**

**Scenario:** A user is studying a metallic system, which is notoriously prone to electronic convergence issues. They have a structure that, with default QE settings, is likely to fail. They want to see the system handle this problem automatically.
> **GIVEN** a configuration file where the error recovery mechanism is explicitly enabled, with `max_retries` set to 2 and a retry strategy that involves lowering the `mixing_beta` parameter from its default of 0.7 to 0.3 on the first retry.
> **AND** an input directory contains a structure for a slab of aluminum, which is known to be somewhat difficult to converge with aggressive electronic mixing settings.
> **AND** for the purpose of the test, the initial Quantum Espresso parameters are deliberately set to be very aggressive to guarantee a first-run SCF failure.
>
> **WHEN** the user runs the labelling pipeline on this structure.
>
> **THEN** the system's log output should clearly indicate that the first DFT run has been launched.
> **AND** the log must then show a message explicitly stating that the run failed, and the reason must be identified as "SCF not converged".
> **AND** immediately following the failure message, the log should show a message indicating that an error recovery attempt is being made, such as "SCF failed. Retrying calculation with modified mixing parameters... (Attempt 1 of 2)".
> **AND** the system should then launch a second DFT calculation. The input file for this second calculation must show that the `mixing_beta` parameter in the `ELECTRONS` namelist has been changed to 0.3.
> **AND** assuming this second run with the more conservative parameters succeeds, the final data should be correctly parsed and saved to the database, and the log should indicate a successful calculation on the second attempt.

### **UAT-C1-004: Verification of Delta Learning**

**Scenario:** A user, who understands the importance of correct physical behavior at short interatomic distances, wants to rigorously verify that the delta learning approach is correctly implemented and that the final model accurately reflects the training data.
> **GIVEN** a complete set of labelled DFT data for a simple diatomic molecule, like Carbon Monoxide (CO), has been generated and stored in the database. This dataset should include structures at, above, and below the equilibrium bond length.
> **AND** the training engine is configured to use a baseline ZBL potential (to handle short-range repulsion) and then train an ACE model on the residual between the ZBL predictions and the DFT data.
>
> **WHEN** the user runs the training pipeline, which completes successfully and produces a final model file, `co_model.ace`.
>
> **THEN** the user should be able to run a provided post-processing validation script that loads this `co_model.ace` and the original DFT dataset from the database.
> **AND** this script will iterate through each structure in the dataset, calculate the energy and atomic forces using the final model, and compare them to the original DFT values.
> **AND** the script's output must show that for all training structures, the root-mean-square error (RMSE) for the energies is less than 1 meV/atom, and the RMSE for the forces is less than 0.01 eV/Angstrom. This confirms that the final combined potential accurately reproduces the training data to a high degree of fidelity.

### **UAT-C1-005: Data Provenance and Persistence**

**Scenario:** Six months after a successful pipeline run, a collaborator asks for the precise computational details used to generate a specific data point in the training set to ensure a fair comparison with their own results.
> **GIVEN** a pipeline run has successfully completed and generated a `results.db` file several months ago.
>
> **WHEN** the user locates this `results.db` file and connects to it using the ASE command-line tool (`ase db results.db`).
>
> **THEN** running `ase db results.db -L` (or similar listing command) should show a list of all the structures that were successfully calculated during that run.
> **AND** when the user inspects a single entry corresponding to a specific structure (e.g., by running `ase db results.db id=17`), the detailed output must contain a `key_value_pairs` or `data` section.
> **AND** this metadata section must contain key-value pairs that explicitly state the DFT parameters used for that specific calculation. This must include, at a minimum: the `ecutwfc` used (e.g., `ecutwfc_ry: 90.0`), the exact names of the pseudopotential files used for each element, and the version of the Quantum Espresso code that was called (e.g., `qe_version: 7.1`). This provides a complete and unambiguous audit trail for the data.
