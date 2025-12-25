# User Acceptance Test (UAT): Cycle 3

**Version:** 1.0.0
**Status:** Ready for Testing

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 3. The primary focus is on `Module B: Explorer & Sampler`, verifying its ability to efficiently explore the conformational space and intelligently select a small, diverse set of structures for expensive DFT calculations. The key outcome is a significant reduction in computational cost compared to the brute-force approach of the previous cycles.

| Scenario ID | Test Scenario Description                                                                | Priority |
|-------------|------------------------------------------------------------------------------------------|----------|
| UAT-C3-001  | **Verify Reduction in DFT Calculations:** The user runs the pipeline on a system. The number of structures selected by Module B for DFT calculation must be exactly the number requested in the configuration and significantly smaller than the total number of frames in the exploratory MD. | High     |
| UAT-C3-002  | **Verify Diversity of Sampled Structures:** The user runs the exploration and sampling module and visualizes the results. The sampled structures, when projected into a low-dimensional space (e.g., using PCA on their descriptors), should be spread out and cover the main regions visited during the exploratory MD, not be clustered in one small area. | High     |
| UAT-C3-003  | **Verify Use of Universal Potential:** The user runs the pipeline and monitors the process. It should be clear from the logs and resource usage (e.g., GPU utilization) that the initial, large-scale MD simulation is being run with a fast, universal MLIP, not with slow DFT calculations. | High     |
| UAT-C3-004  | **Verify Numba Optimization:** The user runs the sampling part of the module on a large trajectory. The descriptor calculation and clustering steps should complete in a reasonably short amount of time, demonstrating the effectiveness of the Numba JIT compilation. | Medium   |
| UAT-C3-005  | **Verify Full Pipeline Integration:** The user runs the complete pipeline (A -> B -> C -> D). The workflow must complete successfully, with the data flowing from the sampler to the labeler and trainer, resulting in a valid final MLIP. | High     |


## 2. Behavior Definitions

This section provides detailed Gherkin-style behavior definitions for the test scenarios, outlining the exact steps and expected outcomes.

---

### **Scenario: UAT-C3-001 - Verify Reduction in DFT Calculations**

*   **GIVEN** an `input.yaml` configured to run an exploratory MD simulation for 10,000 steps.
*   **AND** the configuration specifies that `num_samples_to_select` is `100`.
*   **WHEN** the user executes the full pipeline up to the end of `Module B`.
*   **THEN** the logs should indicate that the exploratory MD produced 10,000 frames.
*   **AND** the logs should indicate that the DIRECT sampler processed these 10,000 frames.
*   **AND** the final list of structures passed to `Module C` for DFT labeling must contain exactly 100 structures.
*   **AND** the system should report that it has reduced the number of structures to be calculated from 10,000 to 100, a reduction of 99%.

---

### **Scenario: UAT-C3-002 - Verify Diversity of Sampled Structures**

*   **GIVEN** a completed run of `Module B` that has produced a set of sampled structures and saved the intermediate data (e.g., the calculated descriptors and the dimensionality reduction results).
*   **WHEN** the user runs a provided plotting script to visualize the results. The script will generate a 2D scatter plot of all MD frames in the PCA/UMAP space, highlighting the 100 selected samples.
*   **THEN** the plot should show a large "cloud" of points representing all the explored configurations.
*   **AND** the 100 highlighted points should be distributed across this cloud, with points appearing in the dense central regions, the sparser outer regions, and any distinct clusters that may have formed.
*   **AND** there should not be a large, obvious region of the explored space from which no points were selected.
*   **AND** the sampled points should not be concentrated in a single small area of the plot, which would indicate a failure of the diversity sampling.

---

### **Scenario: UAT-C3-003 - Verify Use of Universal Potential**

*   **GIVEN** an `input.yaml` configured to use a specific universal potential (e.g., `MACE-MP`).
*   **WHEN** the user runs `Module B` and monitors the system's processes (e.g., using `top` or `nvidia-smi`).
*   **THEN** the logs must explicitly state "Starting exploratory MD with potential: MACE-MP".
*   **AND** a Python process should consume significant CPU or GPU resources.
*   **AND** there should be no `pw.x` or other DFT-related processes running during the exploration phase.
*   **AND** the MD simulation should run significantly faster (e.g., orders of magnitude) than an equivalent AIMD simulation would.
*   **AND** the user should be able to see the MD simulation progressing in real-time through the log output, with regular updates on the simulation time and temperature.

---

### **Scenario: UAT-C3-004 - Verify Numba Optimization**

*   **GIVEN** a large MD trajectory file containing 50,000+ frames.
*   **WHEN** the user runs the sampling part of the code (`_run_direct_sampling`) on this trajectory.
*   **THEN** the process of calculating descriptors for all 50,000 frames should complete in a matter of minutes, not hours.
*   **AND** the logs during this step might show information from Numba about JIT compilation of specific functions.
*   **AND** the user should be able to see a progress bar for the descriptor calculation, indicating that the process is running and not stalled.

---

### **Scenario: UAT-C3-005 - Verify Full Pipeline Integration**

*   **GIVEN** a directory containing only a minimal `input.yaml` for a system like "SiGe".
*   **AND** the pipeline is configured to run all modules (A, B, C, D).
*   **WHEN** the user executes the full pipeline command (`mlip-pipe run input.yaml`).
*   **THEN** the system should execute without crashing.
*   **AND** the logs should clearly show the sequential execution of the modules:
    1.  `Module A` generates initial seeds.
    2.  `Module B` takes the seeds, runs exploratory MD, and samples a small subset.
    3.  `Module C` takes the sampled subset and runs DFT on them.
    4.  `Module D` takes the labeled data and trains the final model.
*   **AND** the final `results.db` should contain the exact number of entries specified in the sampling configuration, not the number of initial seeds or MD frames.
*   **AND** a valid MLIP file must be created in the output directory.
*   **AND** the user should be able to load the final MLIP and use it to run a short MD simulation, confirming that it is a valid and functional potential.
