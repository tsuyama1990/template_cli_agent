# Cycle 3 User Acceptance Testing (UAT)

This document outlines the User Acceptance Testing (UAT) for Cycle 3 of the MLIP-AutoPipe project. The focus of this cycle is on the system's ability to efficiently explore the conformational space and intelligently select data points, which is a critical step in reducing the overall computational cost.

## 1. Test Scenarios

These scenarios test the new exploration and sampling capabilities. The user is an observer in this case, as this module is fully autonomous. The goal is to verify that the system is making intelligent and efficient choices.

| Scenario ID | Scenario Name                                | Priority |
|-------------|----------------------------------------------|----------|
| UAT-C3-001  | Successful Surrogate-Based Exploration       | High     |
| UAT-C3-002  | Correct and Diverse Structure Sampling       | High     |
| UAT-C3-003  | Performance Benchmark of Sampling Process    | Medium   |

---

**Scenario ID: UAT-C3-001**
- **Description:** This test verifies that the system can successfully use a pre-trained foundation model to run a molecular dynamics simulation and generate a trajectory. This is the "Explorer" part of the module. The system should be able to load the surrogate model, run the simulation for a specified number of steps, and save the resulting atomic positions to a trajectory file.
- **Priority:** High
- **Preconditions:**
    - A pre-trained foundation model (e.g., MACE-MP) is downloaded and available.
    - An initial `ase.Atoms` structure is provided.
    - The configuration specifies the MD simulation parameters (temperature, duration).

---

**Scenario ID: UAT-C3-002**
- **Description:** This test verifies that the DIRECT sampling algorithm correctly processes a large trajectory and selects a small, diverse subset of structures. The key outcome is to ensure the selected structures are not all from one part of the simulation (e.g., not all low-energy, equilibrium structures) but represent a good spread of the explored space.
- **Priority:** High
- **Preconditions:**
    - A large trajectory file (e.g., >10,000 frames) from a surrogate MD simulation is available.
    - The configuration for the sampling algorithm (number of clusters, final sample count) is set.

---

**Scenario ID: UAT-C3-003**
- **Description:** This test validates the performance optimisations made using Numba. The processing of a large trajectory file, especially the descriptor calculation, should be reasonably fast. This test ensures that the system is scalable and can handle the large amounts of data generated during the exploration phase without becoming a major bottleneck.
- **Priority:** Medium
- **Preconditions:**
    - A standard benchmark trajectory file (e.g., 50,000 frames) is available.
    - The hardware (CPU cores) for the test is specified.

## 2. Behavior Definitions

---

**Scenario: UAT-C3-001 - Successful Surrogate-Based Exploration**

```gherkin
GIVEN the system has a set of initial structures from Module A
AND the configuration specifies using the 'MACE-MP' foundation model for exploration
AND the configuration sets the MD simulation to run for 1,000 steps at 500K

WHEN the 'Explorer' component of Module B is executed

THEN the process should complete successfully
AND a trajectory file (e.g., 'surrogate_md.traj') should be created on disk
AND the trajectory file should contain exactly 1,000 atomic configurations (frames)
```

---

**Scenario: UAT-C3-002 - Correct and Diverse Structure Sampling**

```gherkin
GIVEN the system has a trajectory file containing 20,000 frames from a melt-quench simulation, which includes solid, liquid, and intermediate states
AND the configuration for the 'Sampler' is set to select 100 final structures

WHEN the 'Sampler' component of Module B is executed

THEN it should return a list of exactly 100 'ase.Atoms' objects
AND a visual inspection of the selected structures (e.g., by plotting their potential energies) should show a wide distribution of energy values, not just values from the low-energy solid state
AND the log files should indicate that the system identified multiple clusters (e.g., more than 10) in the data
```

---

**Scenario: UAT-C3-003 - Performance Benchmark of Sampling Process**

```gherkin
GIVEN a benchmark trajectory file containing 50,000 frames of a 256-atom silicon cell is available
AND the system is running on a standard 8-core CPU

WHEN the 'Sampler' component of Module B is executed on this trajectory

THEN the entire sampling process, from loading the file to selecting the final structures, should complete in under a predefined time limit (e.g., 5 minutes)
AND the log output should confirm that a Numba-accelerated version of the descriptor calculation was used
```
