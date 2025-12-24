# Cycle 4 User Acceptance Testing (UAT)

This document provides the User Acceptance Testing (UAT) plan for Cycle 4. This cycle's primary feature is the fully autonomous, on-the-fly (OTF) active learning loop, which enables the system to improve itself. Tests are designed to verify this self-refinement capability.

## 1. Test Scenarios

These scenarios test the core of the active learning loop and the advanced simulation capabilities. The user's role is to initiate the process and then verify that the system correctly and autonomously improves the MLIP.

| Scenario ID | Scenario Name                                      | Priority |
|-------------|----------------------------------------------------|----------|
| UAT-C4-001  | Successful Detection of an Uncertain Structure     | High     |
| UAT-C4-002  | Completion of a Full Active Learning Cycle         | High     |
| UAT-C4-003  | Successful Run of an Advanced kMC Simulation       | Medium   |

---

**Scenario ID: UAT-C4-001**
- **Description:** This test verifies that the `SimulationEngine` can correctly identify an atomic configuration where the MLIP is uncertain. The system will be baited with an MLIP that is known to be poor for a specific simulation condition. The engine should pause the simulation and correctly extract the structure that caused the high uncertainty.
- **Priority:** High
- **Preconditions:**
    - An MLIP is provided that has been trained only on low-temperature, stable crystal structures.
    - The simulation is configured to run at a high temperature, ensuring the system will explore configurations the MLIP has never seen.
    - The uncertainty threshold is set to a reasonable value.

---

**Scenario ID: UAT-C4-002**
- **Description:** This is the most critical test for this cycle. It verifies that the entire feedback loop works: an uncertain structure is detected, it is sent back for DFT labeling, the MLIP is retrained with the new data, and the simulation resumes with the improved MLIP. The key outcome is to see a measurable improvement in the MLIP's quality.
- **Priority:** High
- **Preconditions:**
    - The full pipeline (Modules C, D, and E) is integrated.
    - An initial, deliberately incomplete MLIP is provided.
    - The system is configured to run for at least two active learning generations.

---

**Scenario ID: UAT-C4-003**
- **Description:** This test verifies that the system can run an advanced simulation, specifically an adaptive Kinetic Monte Carlo (kMC) simulation, to identify rare events. The focus is on ensuring the kMC engine runs, correctly identifies transition states, and evolves the system over a much longer timescale than standard MD.
- **Priority:** Medium
- **Preconditions:**
    - A robust MLIP for a system with known diffusion mechanisms (e.g., a vacancy in a crystal) is provided.
    - The simulation is configured to run in kMC mode.

## 2. Behavior Definitions

---

**Scenario: UAT-C4-001 - Successful Detection of an Uncertain Structure**

```gherkin
GIVEN the system is provided with an MLIP trained only on silicon crystal at 10K
AND the 'SimulationEngine' is configured to run an MD simulation at 2000K (above melting point)
AND the uncertainty threshold is set to 0.1

WHEN the 'SimulationEngine' starts the MD simulation

THEN the simulation should pause automatically before reaching the total number of configured steps
AND the engine's log file should report a message like "Uncertainty threshold exceeded. Extracting structure for retraining."
AND the engine should output a new 'ase.Atoms' object for DFT calculation
AND the atoms in the new structure should be in a disordered, liquid-like configuration, visibly different from the initial crystal
```

---

**Scenario: UAT-C4-002 - Completion of a Full Active Learning Cycle**

```gherkin
GIVEN the system starts with an incomplete MLIP as in UAT-C4-001
AND the full pipeline is configured to run for 2 active learning generations

WHEN the user initiates the active learning workflow

THEN the system should first detect an uncertain structure and trigger Generation 1 of retraining
AND a new DFT calculation should be automatically performed for the uncertain structure
AND a new MLIP, version 2, should be trained using the original data plus the new structure
AND the simulation should then resume using MLIP version 2
AND the log should show that the simulation with MLIP version 2 runs for a longer time before another uncertainty event occurs, or completes without any
AND a final, refined MLIP (version 2 or higher) should be produced at the end
```

---

**Scenario: UAT-C4-003 - Successful Run of an Advanced kMC Simulation**

```gherkin
GIVEN the system is provided with a high-quality MLIP for copper (Cu)
AND the initial structure is a copper crystal containing a single atomic vacancy
AND the 'SimulationEngine' is configured to run in adaptive kMC mode for a simulated time of 1 nanosecond

WHEN the 'SimulationEngine' starts the kMC simulation

THEN the simulation should complete successfully
AND the final atomic structure should show that the vacancy has moved from its original position, indicating a successful diffusion event
AND the log file should contain information about the transition states and event rates found by the kMC algorithm
```
