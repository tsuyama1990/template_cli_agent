# CYCLE03 User Acceptance Testing (UAT)

## 1. Test Scenarios

This UAT plan focuses on verifying the new high-throughput exploration and intelligent sampling capabilities of Module B, as well as the performance enhancements from using `Numba`.

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C03-001 | End-to-End Workflow with Surrogate Exploration | High     |
| UAT-C03-002 | Verification of Data Reduction via Sampling | High     |
| UAT-C03-003 | Performance Gain from Numba JIT             | High     |
| UAT-C03-004 | GPU Acceleration for Surrogate MD           | Medium   |
| UAT-C03-005 | Correct Handling of Missing Surrogate Model | Medium   |

---

**Scenario UAT-C03-001: End-to-End Workflow with Surrogate Exploration**

*   **Description**: This is the primary "happy path" test for the cycle. It verifies that the entire, newly modified workflow runs from start to finish, correctly incorporating the surrogate MD exploration and sampling step.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation from CYCLE02.
    *   The pre-trained surrogate model (e.g., MACE-MP) is downloaded and accessible.
    *   A minimal `input.yaml` file for a system like silicon (`composition: "Si"`).
*   **Acceptance Criteria**:
    *   The command `uv run mlip-pipe run input.yaml` must execute and exit with a status code of 0.
    *   The system log must clearly indicate that the `Explorer & Sampler (Module B)` has started.
    *   The log must show that a surrogate MD simulation is being run for a specified number of steps.
    *   The log must report the total number of frames in the trajectory and the final number of structures selected for DFT (e.g., "Sampled 150 structures from a 100,000 frame trajectory").
    *   The number of new DFT calculations logged to the database must exactly match the number of sampled structures.
    *   The workflow must successfully train and save a final MLIP model based on the intelligently sampled data.

---

**Scenario UAT-C03-002: Verification of Data Reduction via Sampling**

*   **Description**: This test explicitly verifies that the sampling process is working and is effective at reducing the computational load. It checks that the number of expensive DFT calculations is far less than the number of initial candidates.
*   **Preconditions**:
    *   Same as UAT-C03-001.
    *   The `exec_config_dump.yaml` (or `input.yaml`) is configured to run a 10,000-step MD and select 50 final structures.
*   **Acceptance Criteria**:
    *   After the run completes, the log file or terminal output must confirm that the surrogate MD trajectory contained 10,000 frames (or a number close to it).
    *   The log must explicitly state that 50 structures were selected for DFT labeling.
    *   By querying the ASE database, the number of new calculations added during this run must be exactly 50. This confirms the filtering was effective.

---

## 2. Behavior Definitions

**Behavior for UAT-C03-001: End-to-End Workflow with Surrogate Exploration**

```gherkin
Feature: Intelligent Data Selection with a Surrogate Model
  As a researcher,
  I want the system to efficiently explore the material's energy surface using a fast model,
  So that I can reduce the number of expensive DFT calculations and save time and resources.

  Scenario: Running the full workflow with the exploration module enabled
    GIVEN a minimal `input.yaml` for Silicon
    AND a pre-trained MACE surrogate model is available
    WHEN I execute the command `uv run mlip-pipe run input.yaml`
    THEN the process should complete successfully with an exit code of 0
    AND the system log should show the start and end of the surrogate MD simulation
    AND the log should state the number of candidates generated and the smaller number of structures selected
    AND the number of DFT calculations performed should match the number of selected structures
    AND a final MLIP model file should be created.
```

---

**Behavior for UAT-C03-003: Performance Gain from Numba JIT**

```gherkin
Feature: High-Performance Descriptor Calculation
  As a developer,
  I need to ensure that the trajectory analysis is fast enough to handle large datasets,
  So that the sampling step does not become a new bottleneck in the workflow.

  Scenario: Benchmarking the Numba-optimised descriptor calculation
    GIVEN a pre-generated MD trajectory file with at least 5,000 frames
    AND a benchmarking script that can run the descriptor calculation function with and without Numba's JIT compiler
    WHEN I run the benchmark script
    THEN the execution time for the Numba-jitted version must be at least 10 times faster than the pure Python version
    AND the numerical output of both versions must be identical to within a tolerance of 1e-6.
```

---

**Behavior for UAT-C03-004: GPU Acceleration for Surrogate MD**

```gherkin
Feature: GPU Utilization for Accelerated Exploration
  As a user with access to GPU hardware,
  I want the system to automatically use the GPU for surrogate model calculations,
  So that the exploration phase is completed as quickly as possible.

  Scenario: Running the workflow on a machine with a CUDA-enabled GPU
    GIVEN a machine with a compatible NVIDIA GPU and the CUDA toolkit installed
    AND a minimal `input.yaml` for any material
    WHEN I execute the command `uv run mlip-pipe run input.yaml`
    THEN the system log should contain a message indicating that a CUDA device was detected and is being used for the surrogate model
    AND the surrogate MD simulation should complete noticeably faster than on a CPU-only machine.
```
