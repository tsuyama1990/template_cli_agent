# CYCLE 04: UAT.md - Active Learning and On-the-Fly Simulation

## 1. Test Scenarios

The User Acceptance Test for Cycle 04 is designed to showcase the most innovative and powerful feature of the MLIP-AutoPipe system: the fully autonomous, closed-loop active learning process. This UAT will provide a clear and compelling demonstration of the system's ability to "learn on the fly" by intelligently identifying its own weaknesses and automatically seeking new data to improve itself. The user experience will be centered on observing this self-improvement in action.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C04-01** | **Demonstration of a Full Active Learning Loop Iteration** | **High** |

**Scenario UAT-C04-01: Demonstration of a Full Active Learning Loop Iteration**

This scenario will simulate a realistic research problem: modeling the behavior of point defects in a crystal, a classic case where a potential trained only on bulk crystal data will fail. The UAT will begin with an MLIP that has only been trained on a perfect silicon crystal. The user will then instruct the system to run a molecular dynamics simulation on a silicon crystal containing a vacancy (a missing atom).

The user will observe the following autonomous behavior:
1.  **Simulation & Detection:** The system starts the simulation. As atoms near the vacancy move into configurations the model has never seen, the model's uncertainty will rise. The system will automatically detect this and pause.
2.  **Extraction & Labeling:** The system will extract the local atomic environment around the vacancy, send it to the (mocked) DFT engine for labeling, and add the new, labeled data point to the database.
3.  **Retraining:** The system will automatically retrain the MLIP, incorporating the new information about the vacancy defect.
4.  **Convergence:** The system will then restart the simulation. This time, because it has learned about the vacancy, it will be able to continue the simulation without triggering further uncertainty, demonstrating that it has successfully improved itself and "converged."

The `C04_UAT.ipynb` Jupyter Notebook will provide a narrative for this process, allowing the user to trigger the workflow and then use pre-written cells to inspect the database and log files at each stage, verifying that the retraining loop was successfully triggered and completed.

## 2. Behavior Definitions

**Feature:** On-the-fly active learning loop for autonomous model improvement.

**Scenario:** The system successfully identifies an unknown defect structure during a simulation, labels it, and retrains the model.

*   **GIVEN** an ASE database (`uat_c04.db`) that contains labeled data for a **perfect** 64-atom silicon crystal.
*   **GIVEN** a trained MLIP model file (`model_initial.ace`) that was trained only on this perfect crystal data.
*   **GIVEN** a starting `Atoms` object for the simulation which is the same 64-atom silicon crystal but with **one atom removed** (a vacancy).
*   **GIVEN** a configuration file that specifies a maximum of 3 active learning generations.
*   **GIVEN** a mock MLIP calculator that is programmed to return a high uncertainty value whenever it encounters an atom with fewer than 4 neighbors (i.e., the atoms around the vacancy).
*   **GIVEN** a mock of the Quantum Espresso executable.

*   **WHEN** the user executes the active learning workflow with the vacancy structure as input.

*   **THEN** the command should execute and exit with a success code (0).
*   **AND** the log output should clearly show the start of "Active Learning Generation 1".
*   **AND** during Generation 1, the log should contain a message like "High uncertainty detected. Extracting structure for re-labeling."
*   **AND** a new structure should be added to the `uat_c04.db` database with the status 'unlabeled_active'.
*   **AND** the `LabelingEngine` should run and update this new structure's status to 'labeled'.
*   **AND** the `TrainingEngine` should run, creating an updated model file (`model_gen1.ace`).
*   **AND** the log output should then show the start of "Active Learning Generation 2".
*   **AND** during Generation 2, the simulation should run to completion **without** detecting high uncertainty.
*   **AND** the final log message should indicate convergence, e.g., "No new uncertain structures found. Active learning converged."
