# CYCLE04/UAT.md

## 1. Test Scenarios

User Acceptance Testing (UAT) for Cycle 4 is critically important as it validates the transition from a linear pipeline to a dynamic, self-improving system. The core user story for this cycle is: "As a computational materials scientist, I want the system to not only build an initial potential but also to iteratively refine it by running simulations, detecting where the model is uncertain, and automatically acquiring new DFT data to patch these knowledge gaps." This UAT will be conducted via a Jupyter Notebook (`uat_cycle04.ipynb`) that allows the user to monitor the active learning process across several generations, providing tangible proof of the system's ability to learn and improve autonomously.

| Scenario ID | Scenario Description                                                                | Priority |
| :---------- | :---------------------------------------------------------------------------------- | :------- |
| UAT-C04-01  | **Verify Active Learning Loop Execution**                                           | High     |
|             | The user will run the pipeline with the active learning mode enabled for a small number of generations (e.g., 3). The primary goal is to verify that the system correctly cycles through the `Simulate -> Detect Uncertainty -> Re-label -> Re-train` loop for the specified number of iterations without crashing. |          |
| UAT-C04-02  | **Confirm Uncertainty-Driven Data Acquisition**                                     | High     |
|             | The user will inspect the ASE database after the active learning run. They will verify that new atomic structures were added to the database during the simulation phase and that these new structures were subsequently labeled with DFT results. This provides direct evidence that the uncertainty mechanism is successfully identifying and acquiring new data. |          |
| UAT-C04-03  | **Demonstrate Model Improvement Across Generations**                                | Medium   |
|             | The user will load the MLIP models generated at the end of each active learning generation. They will test the models' performance on a small, independent "holdout" set of structures. The expectation is to see a measurable improvement (i.e., a decrease in prediction error) in the models from later generations compared to the initial model, demonstrating that the system is truly learning. |          |

## 2. Behavior Definitions

The behavior of the system will be defined using the Gherkin-style (GIVEN/WHEN/THEN) syntax. These definitions will be implemented as interactive steps within the `uat_cycle04.ipynb` notebook.

### Scenario: UAT-C04-01 - Verify Active Learning Loop Execution

*   **GIVEN** a minimal `input.yaml` for a simple material.
*   **AND** a configuration file that enables active learning for 3 generations (e.g., `active_learning: {max_generations: 3}`).
*   **WHEN** the user executes the main CLI command `uv run mlip-pipe run --input input.yaml`.
*   **THEN** the command should execute and finish without unhandled exceptions.
*   **AND** the system's log output should clearly show the progression through the learning loop. The user should see distinct log blocks for "Running Active Learning Generation 1", "Running Active Learning Generation 2", and "Running Active Learning Generation 3".
*   **AND** within each generation block, the logs should show the sequence of "Running simulation engine", "New uncertain structures found", "Running labeling engine", and "Running training engine".
*   **AND** three distinct model files should be created, e.g., `trained_model_gen0.pt`, `trained_model_gen1.pt`, and `trained_model_gen2.pt`.

### Scenario: UAT-C04-02 - Confirm Uncertainty-Driven Data Acquisition

*   **GIVEN** the successful completion of the workflow from scenario UAT-C04-01.
*   **AND** the initial training set (after the exploration phase) contained N structures.
*   **WHEN** the user connects to the output ASE database (`mlip_autopipec.db`) and queries for all structures that have a DFT label.
*   **THEN** the total number of labeled structures in the database should be greater than N.
*   **AND** the user will query the database for structures added during "generation 1".
*   **AND** these structures should have a non-null value for their DFT energy and forces, confirming they were successfully passed to the labeling engine and processed.

### Scenario: UAT-C04-03 - Demonstrate Model Improvement Across Generations

*   **GIVEN** the successful completion of the workflow from scenario UAT-C04-01.
*   **AND** the user has prepared a small, independent set of "holdout" atomic structures that were not included in any training.
*   **AND** the user has loaded the three trained models: `model_gen0`, `model_gen1`, and `model_gen2`.
*   **WHEN** the user calculates the Root Mean Square Error (RMSE) of the force predictions for the holdout set using `model_gen0`.
*   **AND** the user calculates the force RMSE for the holdout set using `model_gen1`.
*   **AND** the user calculates the force RMSE for the holdout set using `model_gen2`.
*   **THEN** the RMSE of `model_gen1` should ideally be lower than the RMSE of `model_gen0`.
*   **AND** the RMSE of `model_gen2` should ideally be lower than the RMSE of `model_gen1`.
*   **AND** the notebook will generate a bar chart showing the RMSE for each model generation, providing a clear visual representation of the model's improvement over time. (Note: Due to the stochastic nature of training, a slight increase is possible, but a general downward trend is expected).
