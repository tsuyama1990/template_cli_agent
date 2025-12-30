# Cycle 04 User Acceptance Test (UAT): Active Learning and Advanced Simulations

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 04. The primary focus is on verifying the newly implemented active learning loop and the system's ability to run long-timescale simulations autonomously from a user's perspective. These tests are designed to confirm that the pipeline can correctly identify its own weaknesses by quantifying uncertainty, gather new data to address those weaknesses, and iteratively improve its own predictive accuracy in a fully automated closed loop. The successful completion of these tests will demonstrate that the system has achieved its core design goal of becoming a truly autonomous and intelligent agent for developing machine learning potentials.

| Scenario ID | Description                                                              | Priority |
|-------------|--------------------------------------------------------------------------|----------|
| UAT-04-001  | **Successful Uncertainty-Triggered Retraining:** Verify that during a simulation, the system can detect a high-uncertainty structure, pause the simulation, correctly label the new structure via DFT, and retrain the MLIP. This is the most critical test, validating the fundamental mechanism of the active learning feedback loop. | High     |
| UAT-04-002  | **Correct Dynamic Threshold Adjustment:** Verify that the uncertainty threshold used to trigger retraining becomes stricter and more refined as the model is trained on more data across multiple generations. This test ensures the system's "curiosity" adapts appropriately as its knowledge grows. | High     |
| UAT-04-003  | **Full Active Learning Loop Convergence:** Verify that the entire active learning loop can run for multiple generations and eventually terminate gracefully when no more high-uncertainty structures are found during a full simulation run. This test confirms that the process is stable and has a natural, successful end-state. | High     |
| UAT-04-004  | **Physically Meaningful Structure Extraction:** Verify that the boundary treatment process for structures extracted from periodic simulations results in a chemically sensible, non-periodic cluster that is suitable for a subsequent DFT calculation. This is crucial for ensuring that the new data being added is physically valid and not corrupted by simulation artefacts. | Medium   |

## 2. Behaviour Definitions

The following Gherkin-style definitions describe the expected behaviour for each test scenario in detail.

---

### **Scenario: UAT-04-001 - Successful Uncertainty-Triggered Retraining**

This scenario tests the single most important feature of this cycle: the correct functioning of the feedback loop. It verifies that all the components (simulation, uncertainty calculation, structure extraction, labelling, and retraining) are correctly wired together. The test uses a high temperature simulation to deliberately induce atomic configurations that are likely to be outside the initial training set's experience, providing a predictable way to trigger the uncertainty mechanism.

**GIVEN** an MLIP has been initially trained on a small set of equilibrium and slightly strained structures for bulk silicon (the "generation 0" model).
**WHEN** the user starts the `SimulationEngine` to run a molecular dynamics simulation at a high temperature (e.g., 1200K, which is near the melting point of silicon) where significant atomic displacements are expected.
**THEN** the system should log that it is running an On-the-Fly (OTF) simulation and is actively monitoring the uncertainty at each step.
**AND** after some number of simulation steps, the system should log a clear message indicating that a high-uncertainty structure has been detected because its calculated uncertainty score exceeded the current threshold.
**AND** a new structure, corresponding to the atomic configuration that triggered the event, should be added to the ASE database with the state "unlabelled".
**AND** the main orchestrator should detect the presence of this new unlabelled structure and automatically trigger the `LabellingEngine`.
**AND** the `LabellingEngine` should then compute the DFT properties of this new structure and update its state in the database to "labelled".
**AND** the orchestrator should then automatically trigger the `TrainingEngine`.
**AND** the `TrainingEngine` should retrain the MLIP using the original dataset plus this new, high-uncertainty data point.
**AND** a new, updated MLIP model file (e.g., `model_gen1.ace`) should be saved to disk, representing the improved, "generation 1" model.

---

### **Scenario: UAT-04-02 - Correct Dynamic Threshold Adjustment**

This scenario tests the intelligence of the uncertainty mechanism. A fixed, user-defined uncertainty threshold is brittle. A good active learning system should adapt its threshold as it improves. Initially, many structures may seem novel, so the threshold should be high. After retraining, the model's domain of applicability has expanded, so the threshold should become stricter to find the next, more subtle, areas of weakness.

**GIVEN** a system that is configured to run an active learning workflow for at least two generations.
**WHEN** the first simulation generation ("generation 0") starts.
**THEN** the system should log the initial uncertainty threshold it is using, which has been calculated from the initial training set. The log message should be specific, e.g., "Dynamic uncertainty threshold for generation 0 set to 0.53 eV/A".
**AND** after the system finds a new structure, labels it, and retrains the model to create the "generation 1" potential.
**WHEN** the second simulation generation ("generation 1") starts, using this new potential.
**THEN** the system should re-calculate the threshold based on the now larger training set and log the new value it is using, e.g., "Dynamic uncertainty threshold for generation 1 set to 0.45 eV/A".
**AND** the new threshold Y should be different from, and typically stricter (lower) than, the initial threshold X. This demonstrates that the system is correctly adapting its criterion for novelty as its knowledge base expands.

---

### **Scenario: UAT-04-03 - Full Active Learning Loop Convergence**

This scenario tests the stability and end-state of the active learning process. A successful active learning run should not continue forever; it should reach a point where the MLIP is sufficiently robust for the simulated conditions that it no longer encounters configurations it considers novel. This test verifies that the system can reach this converged state and terminate gracefully.

**GIVEN** a user provides a minimal `input.yaml` for a simple, well-behaved system (e.g., bulk Aluminium).
**AND** the active learning is configured with a maximum of 10 generations and a total simulation time.
**WHEN** the user executes the full end-to-end workflow.
**THEN** the system should proceed through multiple generations of the active learning loop, with the log file showing the generation number incrementing: "Starting active learning generation 1", "Starting active learning generation 2", etc.
**AND** in the early generations, the log should show that high-uncertainty structures are being found.
**AND** eventually, after several retraining cycles, a simulation run should complete its entire allotted time without finding any atomic configurations that exceed its (now much stricter) uncertainty threshold.
**AND** upon completing this simulation, the system should log a clear and explicit message indicating that the active learning process has converged, e.g., "No new high-uncertainty structures found in the last simulation. Active learning has converged successfully."
**AND** the entire workflow should then terminate gracefully, well before reaching the configured maximum of 10 generations.

---

### **Scenario: UAT-04-04 - Physically Meaningful Structure Extraction**

This scenario tests a crucial but subtle part of the active learning process: the correct handling of periodic boundary conditions. When a high-uncertainty event is detected in a bulk simulation, simply saving the entire periodic cell is not ideal for retraining. The uncertainty is often localised to a specific atom, and a smaller, non-periodic cluster around that atom is a much more efficient structure for a DFT calculation. This test ensures this extraction is done correctly.

**GIVEN** a simulation of a covalent material (e.g., amorphous silica) is running in a periodic cell.
**AND** the system detects a high-uncertainty event and identifies that the atom with the highest uncertainty is near a periodic boundary of the simulation cell.
**WHEN** the system extracts the atomic neighbourhood around this atom to create a new training structure.
**THEN** the resulting `ase.Atoms` object that is saved to the database should be non-periodic (its `pbc` attribute should be `[False, False, False]`).
**AND** a visual inspection of the saved structure (e.g., by exporting it to an XYZ file) should show a sensible, roughly spherical cluster of atoms.
**AND** crucially, any Si-O bonds that were cut by the periodic boundary during the extraction should be correctly terminated with passivating atoms (e.g., hydrogen atoms).
**AND** the resulting cluster should be chemically sensible and should not contain any unphysical, dangling bonds, making it a valid input for a high-quality quantum chemical calculation.