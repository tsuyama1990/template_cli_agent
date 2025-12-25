# User Acceptance Test (UAT): Cycle 4

**Version:** 1.0.0
**Status:** Ready for Testing

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 4. The focus is on `Module E: Simulation Engine` and the correct implementation of the "on-the-fly" (OTF) active learning feedback loop. The tests are designed to verify that the system can detect its own uncertainty during a simulation, pause, trigger re-training, and resume with an improved model.

| Scenario ID | Test Scenario Description                                                                | Priority |
|-------------|------------------------------------------------------------------------------------------|----------|
| UAT-C4-001  | **Verify Uncertainty Trigger:** The user runs a simulation that is designed to enter a state not represented in the initial training data. The system must detect high uncertainty, pause the simulation, and log the event. | High     |
| UAT-C4-002  | **Verify Re-training and Resumption:** Following a triggered uncertainty event, the system must correctly label the new structure, re-train the MLIP, and resume the simulation from the point where it was paused, now using the new, improved potential. | High     |
| UAT-C4-003  | **Verify Correct Structure Extraction:** When an uncertainty event is triggered, the user inspects the structure that is extracted for re-training. It must be a valid, periodic atomic configuration that correctly represents the state that caused the trigger. | High     |
| UAT-C4-004  | **Verify Successful Completion of a Full OTF Run:** The user runs a complete active learning workflow on a simple system. The simulation should start, trigger one or more re-training cycles, and eventually complete its full simulation time without crashing, resulting in a "simulation-hardened" final MLIP. | High     |
| UAT-C4-005  | **Verify Stability Improvement:** The user compares two simulations: one run with the initial "generation 0" potential, and one run with the final potential after several active learning cycles. The simulation with the final potential should be noticeably more stable and produce more physical trajectories. | Medium   |


## 2. Behavior Definitions

This section provides detailed Gherkin-style behavior definitions for the test scenarios, outlining the exact steps and expected outcomes.

---

### **Scenario: UAT-C4-001 - Verify Uncertainty Trigger**

*   **GIVEN** an initial training set and a trained "generation 0" MLIP for Silicon that contains only bulk-like structures at normal density.
*   **AND** a simulation is configured to run MD in an NPT ensemble where the pressure is set very high, forcing the system to compress.
*   **AND** the `uncertainty_threshold` is set to a reasonable value.
*   **WHEN** the user launches the OTF simulation.
*   **THEN** the simulation should run for some number of steps.
*   **AND** as the simulation box volume decreases, the system must eventually pause.
*   **AND** a log message must be generated, stating "Uncertainty threshold exceeded. Pausing simulation for re-training."
*   **AND** the structure at the moment of pausing must have a significantly smaller volume than the starting structure.
*   **AND** the system should report which atom(s) triggered the uncertainty, and the reason for the uncertainty (e.g., "novel local environment").

---

### **Scenario: UAT-C4-02 - Verify Re-training and Resumption**

*   **GIVEN** the system has just completed the actions from scenario UAT-C4-001 and is paused.
*   **WHEN** the user allows the active learning pipeline to proceed.
*   **THEN** the system must automatically execute the DFT labeling for the new compressed structure.
*   **AND** the system must automatically execute the training engine to produce a new "generation 1" MLIP file.
*   **AND** the simulation must then resume automatically from the exact frame where it was paused.
*   **AND** the logs must clearly indicate that it is now running with the "generation 1" potential.
*   **AND** the simulation should now be able to continue running in the high-pressure compressed state without immediately triggering another uncertainty event.
*   **AND** the user should be able to see the new data point being added to the training set in the log output.

---

### **Scenario: UAT-C4-03 - Verify Correct Structure Extraction**

*   **GIVEN** the system has just been paused due to an uncertainty trigger.
*   **WHEN** the user inspects the structure file (e.g., a CIF or POSCAR) that was saved for the DFT labeling engine.
*   **THEN** the file must contain a valid atomic structure with the correct number of atoms and species.
*   **AND** the lattice parameters and atomic positions in the file must match those of the simulation at the moment it was paused.
*   **AND** the structure must have periodic boundary conditions, confirming the "periodic embedding" extraction was successful.
*   **AND** the user should be able to visualize the extracted structure and confirm that it is a physically reasonable configuration.

---

### **Scenario: UAT-C4-04 - Verify Successful Completion of a Full OTF Run**

*   **GIVEN** a minimal `input.yaml` for a simple system.
*   **AND** the simulation is configured to run for a total of 10 picoseconds.
*   **AND** the active learning loop is configured with `max_active_learning_cycles = 5`.
*   **WHEN** the user runs the full OTF workflow.
*   **THEN** the system must run to completion without crashing.
*   **AND** the final log output must state "Simulation completed successfully."
*   **AND** the output directory should contain the final trajectory file covering the full 10 picoseconds.
*   **AND** the output directory should also contain several generations of the MLIP (e.g., `model_gen0.ace`, `model_gen1.ace`, `model_gen2.ace`), indicating that the re-training loop was successfully triggered.
*   **AND** the user should be able to see a summary of the active learning process, including the number of cycles, the number of new data points added, and the final uncertainty level.

---

### **Scenario: UAT-C4-05 - Verify Stability Improvement**

*   **GIVEN** the "generation 0" MLIP and the final "generation N" MLIP from a completed OTF run.
*   **WHEN** the user runs two separate, identical, long MD simulations, one with each potential.
*   **THEN** the simulation run with the "generation 0" potential might show signs of instability, such as large energy drifts or unphysical atomic configurations (e.g., "atom fusion").
*   **AND** the simulation run with the final "generation N" potential must show stable energy conservation and maintain a physically realistic trajectory for the entire duration.
*   **AND** the user should be able to plot the potential energy over time for both simulations, and the plot for the final potential should show much smaller fluctuations and no unphysical energy drifts.
