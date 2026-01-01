# CYCLE 03: UAT.md - Surrogate-Based Exploration and Sampling

## 1. Test Scenarios

This User Acceptance Test (UAT) for Cycle 03 focuses on validating the newly implemented `Explorer & Sampler` module. The core objective is to demonstrate the significant gain in computational efficiency achieved by using a surrogate model for exploration before committing to expensive DFT calculations. The user experience will highlight the "smart" data selection process, contrasting the large number of configurations explored with the small, information-rich subset that is ultimately chosen for labeling.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C03-01** | **Verification of Intelligent Sampling for Cost Reduction**| **High** |

**Scenario UAT-C03-01: Verification of Intelligent Sampling for Cost Reduction**

This scenario will use amorphous Silicon (a-Si) as a test case, a system where structural diversity is key to creating a good potential. The user will start with the same minimal input as in previous cycles. They will observe the system generate a few initial structures, but then, critically, they will see the new `Explorer & Sampler` module in action. This module will run a short molecular dynamics simulation using the fast MACE surrogate model, exploring thousands of different atomic configurations. The user will then be able to verify that from this vast pool of candidates, the system intelligently selects only a small, representative number (e.g., 50) for the expensive DFT labeling step.

The `C03_UAT.ipynb` Jupyter Notebook will provide a clear, step-by-step demonstration of this process. It will allow the user to:
1.  Define the minimal `input.yaml` for Silicon.
2.  Run the pipeline via the CLI.
3.  Inspect the database and log files to explicitly see the number of initial structures, the number of total frames explored by the MD simulation, and the final, much smaller number of structures that were selected for labeling.
This provides tangible proof of the system's efficiency, showing how it avoids redundant calculations and focuses resources where they are most needed, a cornerstone of the project's design philosophy.

## 2. Behavior Definitions

**Feature:** Surrogate-based exploration and intelligent data sampling.

**Scenario:** Successful and efficient data selection for an amorphous Silicon system.

*   **GIVEN** a clean working directory.
*   **GIVEN** a minimal YAML configuration file named `input.yaml` with the following content:
    ```yaml
    system:
      elements: ["Si"]
    simulation:
      temperature: [300, 1000]
    ```
*   **GIVEN** an expanded configuration that specifies running a 100-step MD simulation for exploration and selecting `n_samples: 50` final structures.
*   **GIVEN** an empty ASE database file (`uat_c03.db`).
*   **GIVEN** a mock of the Quantum Espresso executable.

*   **WHEN** the user executes the command `uv run mlip-pipe --config-file input.yaml --database-file uat_c03.db`.

*   **THEN** the command should execute and exit with a success code (0).
*   **AND** the log output should clearly indicate the start and end of the `Explorer & Sampler` stage.
*   **AND** the log output should contain a message similar to: "Explorer generated 100 trajectory frames."
*   **AND** the log output should contain a message similar to: "DIRECT sampling selected 50 representative structures for labeling."
*   **AND** after the `Explorer & Sampler` stage, the `uat_c03.db` database should contain some structures marked as 'generated' or 'archived'.
*   **AND** the `uat_c03.db` database should contain exactly 50 structures marked as 'selected_for_labeling'.
*   **AND** the subsequent log messages for the `LabelingEngine` should indicate that it is processing exactly 50 structures.
*   **AND** a final trained model file named `model.ace` should be created.
