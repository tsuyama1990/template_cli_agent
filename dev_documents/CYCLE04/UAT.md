# CYCLE04 User Acceptance Testing (UAT)

## 1. Test Scenarios

This UAT plan focuses on verifying the new autonomous, self-improving capabilities of the system, primarily the On-the-fly (OTF) active learning loop and the integration of the LAMMPS simulation engine.

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C04-001 | Successful End-to-End Active Learning Loop  | High     |
| UAT-C04-002 | Simulation Pause and Resume Integrity       | High     |
| UAT-C04-003 | Correct Structure Extraction on Uncertainty | Medium   |
| UAT-C04-004 | Incremental Model Retraining Verification   | High     |
| UAT-C04-005 | Adaptive kMC for Rare Event Discovery       | Medium   |

---

**Scenario UAT-C04-001: Successful End-to-End Active Learning Loop**

*   **Description**: This is the most critical test for this cycle. It verifies that the system can autonomously improve its own potential. The test will start with a deliberately incomplete MLIP (e.g., trained only on perfect crystal structures) and run a simulation at high temperature, forcing the system to encounter new, unseen configurations.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation from CYCLE03.
    *   A compatible LAMMPS executable is available in the system's PATH.
    *   An initial MLIP has been trained on a very limited dataset (e.g., only a single Si crystal structure).
    *   The configuration is set for an MD simulation at a high temperature (e.g., 2000K) with the active learning loop enabled.
*   **Acceptance Criteria**:
    *   The command to run the simulation must execute and run without crashing.
    *   The system log must show that the simulation started with the initial model (`model_v1`).
    *   The log must show that high uncertainty was detected, triggering the pause of the simulation.
    *   The log must confirm that a new structure was extracted and sent for DFT calculation.
    *   The database must show a new entry corresponding to this extracted structure.
    *   The log must show that the training engine was invoked to create an improved model (`model_v2`).
    *   The log must show that the simulation was successfully resumed using `model_v2`.
    *   The system should be able to run for a longer time with `model_v2` before the next uncertainty trigger, demonstrating improvement.

---

**Scenario UAT-C04-002: Simulation Pause and Resume Integrity**

*   **Description**: This test focuses on the technical robustness of the control mechanism over the external LAMMPS process. It ensures that pausing and resuming the simulation does not lead to corrupted data or an inconsistent simulation state.
*   **Preconditions**:
    *   Same as UAT-C04-001.
    *   The uncertainty threshold is set very low to guarantee an early trigger.
*   **Acceptance Criteria**:
    *   The LAMMPS subprocess must be verifiably paused (e.g., `SIGSTOP`) and resumed (`SIGCONT`).
    *   After resuming, the physical properties of the simulation (e.g., total energy, temperature) must be consistent with the values just before the pause, showing no unphysical jumps.
    *   The final trajectory file should be continuous and show no temporal or spatial gaps corresponding to the pause.

---

## 2. Behavior Definitions

**Behavior for UAT-C04-001: Successful End-to-End Active Learning Loop**

```gherkin
Feature: Autonomous Potential Refinement via Active Learning
  As a scientist,
  I want the system to automatically detect its own weaknesses during a simulation and improve itself,
  So that I can generate a robust potential without manual trial-and-error.

  Scenario: Running a high-temperature simulation with a weak initial potential
    GIVEN a system with a trained, but incomplete, MLIP for Silicon
    AND the active learning loop is enabled with a reasonable uncertainty threshold
    WHEN I launch a high-temperature Molecular Dynamics simulation
    THEN the system should start the simulation using the initial MLIP
    AND the system log must show the detection of high uncertainty, triggering a pause
    AND a new atomic configuration should be extracted and submitted for a DFT calculation
    AND a new, improved MLIP should be trained incorporating the new data point
    AND the simulation must automatically resume from where it left off, using the improved MLIP.
```

---

**Behavior for UAT-C04-004: Incremental Model Retraining Verification**

```gherkin
Feature: Efficient Model Retraining
  As a user,
  I want the retraining process during active learning to be as fast as possible,
  So that the simulation can resume quickly.

  Scenario: The system adds a new data point to the training set
    GIVEN a simulation is running and the active learning loop is triggered
    WHEN the system invokes the Training Engine to incorporate the new data point
    THEN the system log should indicate that it is performing an "incremental update" or "fine-tuning" of the existing model
    AND the time taken for this retraining step should be significantly shorter (e.g., >5x faster) than the time it took to train the original model from scratch.
```

---

**Behavior for UAT-C04-005: Adaptive kMC for Rare Event Discovery**

```gherkin
Feature: Long-Timescale Simulation with Kinetic Monte Carlo
  As a researcher studying material defects,
  I want to simulate rare atomic events like diffusion,
  So that I can understand long-term material evolution.

  Scenario: Simulating vacancy diffusion in a crystal
    GIVEN a trained MLIP for a metal (e.g., Aluminum)
    AND an initial structure of an Aluminum crystal containing a single vacancy
    AND the configuration is set to run an "Adaptive kMC" simulation
    WHEN I launch the simulation
    THEN the system should successfully identify potential vacancy hop events and their corresponding saddle points
    AND the kMC simulation should proceed, correctly simulating the vacancy moving through the crystal lattice over an extended period of simulation time (e.g., microseconds)
    AND the log should report the events being discovered and the advancement of simulation time.
```
