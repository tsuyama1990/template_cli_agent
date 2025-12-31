# CYCLE04/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 4 of the MLIP-AutoPipe project. With the preceding cycles, we have built a powerful linear pipeline capable of generating, exploring, and training on a static dataset. This cycle marks a pivotal evolution, transforming that linear process into a dynamic, self-improving, closed-loop system. The primary objective is to implement a full active learning workflow. This is the ultimate expression of the system's core philosophy of "removing the human expert from the loop," as the system will now be able to identify and learn from its own areas of ignorance autonomously.

The centerpiece of this cycle is the implementation of **Module E: Simulation Engine**. This module will take the MLIP trained in Module D and use it to run large-scale molecular dynamics simulations. The key innovation is that this will not be a passive simulation. Module E will continuously perform "on-the-fly" (OTF) uncertainty quantification. It will monitor the predictions of the MLIP, and if it encounters an atomic configuration where the model's uncertainty exceeds a dynamic threshold, it will automatically flag this structure as a point of high information value. This "interesting" structure will then be captured and sent back to the DFT Labeling Engine (Module C) for calculation. The `Orchestrator` will be significantly upgraded to manage this cyclical workflow, feeding the newly labeled data back into the Training Engine (Module D) to retrain and improve the MLIP. A crucial technical challenge to be addressed is the proper handling of periodic boundary conditions, ensuring that structures extracted from simulations are physically meaningful. By the end of this cycle, MLIP-AutoPipe will be a true learning machine, capable of iteratively and intelligently refining its own potential to become progressively more accurate and robust.

## 2. System Architecture

The architecture for Cycle 4 undergoes a significant change, moving from a linear pipeline to a cyclical graph. The new `SimulationEngine` module becomes the driver of the learning loop, and the `Orchestrator`'s logic becomes more complex to manage the iterative nature of the workflow.

**File Structure:**

The following ASCII tree highlights the new or significantly modified files for this cycle. The main additions are the new Module E and the corresponding updates to the orchestrator and configuration.

```
.
├── src/
│   └── mlip_autopipec/
│       ├── data/
│       │   └── models.py     # Modified to include active learning config
│       ├── modules/
│       │   ├── ...
│       │   └── **e_simulation_engine.py**
│       ├── orchestrator.py # Heavily modified to manage the active learning loop
│       └── utils/
│           ├── **boundary_utils.py** # Utilities for handling periodic boundaries
│           └── ...
├── tests/
│   ├── unit/
│   │   ├── modules/
│   │   │   └── **test_e_simulation_engine.py**
│   │   └── utils/
│   │       └── **test_boundary_utils.py**
│   └── e2e/
│       └── **test_cycle04_workflow.py**
└── ...
```

**Component Blueprint:**

*   **`modules/e_simulation_engine.py`**: This new file will house the `SimulationEngine` class. This class is responsible for running production simulations and identifying uncertain structures. Its `execute` method will:
    1.  Load the latest trained MLIP model.
    2.  Set up and run a large-scale MD simulation.
    3.  During the simulation, at regular intervals, calculate the uncertainty of the MLIP's predictions for the current atomic configuration. A common method for MACE models is to use the disagreement between the multiple feature readouts as a proxy for uncertainty.
    4.  Compare the uncertainty against a dynamic threshold.
    5.  If the threshold is exceeded, it will pause the simulation, extract the uncertain structure, and use helper functions from `boundary_utils.py` to prepare it for a DFT calculation (e.g., creating a suitable cluster with passivation).
    6.  Save the prepared structure to the `AseDB` with a status indicating it's a new candidate for active learning.
    7.  The method will return a status indicating whether new structures were found, which will signal the orchestrator to continue the learning loop.
*   **`utils/boundary_utils.py`**: This new utility module will provide functions to correctly handle the extraction of atomic structures from simulations with periodic boundary conditions. A key function, `extract_and_passivate_cluster`, will take a full periodic cell and an atom index, cut out a spherical cluster around that atom, and apply appropriate boundary treatment (e.g., adding hydrogen atoms to terminate any broken covalent bonds) to create a physically realistic, non-periodic structure suitable for a DFT calculation.
*   **`orchestrator.py` (Heavily Modified)**: The `Orchestrator` will be refactored to manage the active learning loop. Instead of a single `run_full_pipeline` method, its logic will be structured around a loop:
    ```
    run_initial_pipeline() # A -> B -> C -> D
    for i in range(max_active_learning_cycles):
        found_new = SimulationEngine.execute()
        if found_new:
            LabelingEngine.execute()
            TrainingEngine.execute()
        else:
            break # Converged
    ```
    It will be responsible for managing the state across iterations, such as loading the correct model for each new simulation run.
*   **`data/models.py` (Modified)**: The Pydantic configuration will be updated to include parameters that control the active learning loop.

## 3. Design Architecture

The design for Cycle 4 focuses on creating a robust and configurable active learning loop. The Pydantic schemas will be extended to give the user control over the simulation and the stopping criteria for the iterative process.

**Pydantic Schema Design:**

*   **`ActiveLearningConfig` (Model)**: A new top-level model to group active learning parameters.
    *   `max_generations`: An integer setting the maximum number of active learning cycles (e.g., 10).
    *   `uncertainty_threshold`: A string or float. Can be a fixed value or a dynamic strategy like `"dynamic_95percentile"`, which sets the threshold based on the uncertainty distribution of the existing training set.
    *   `md_steps_per_generation`: The number of MD steps to run in Module E before checking for new structures.

*   **`SimulationEngineConfig` (Model, part of FullConfig)**: Encapsulates settings for the simulation itself.
    *   `simulation_temperature`: The temperature for the main MD simulation.
    *   `pressure`: The pressure for the simulation.

*   **`FullConfig` (Top-level Model, modified)**:
    *   ... (previous sections)
    *   **`simulation`**: An instance of `SimulationEngineConfig`.
    *   **`active_learning`**: An instance of `ActiveLearningConfig`.

**Data Flow and Consumers/Producers:**

The data flow in Cycle 4 is cyclical, orchestrated by the main controller.

1.  **Initial Model**: The `SimulationEngine` (Module E) **consumes** the MLIP model artifact produced by the `TrainingEngine` (Module D).
2.  **Simulation Trajectory**: During its run, Module E **produces** a simulation trajectory.
3.  **Uncertainty Calculation**: The engine internally **consumes** its own trajectory frames and the MLIP model to **produce** an uncertainty score for each frame.
4.  **Uncertain Structures**: When the uncertainty score exceeds the threshold, the engine **produces** one or more new `Atoms` objects. These are extracted and processed using `boundary_utils.py`.
5.  **Database Update**: These new, uncertain structures are saved to the `AseDB`, where they are now available for **consumption** by the `LabelingEngine` (Module C).
6.  **Loop Trigger**: The `Orchestrator` detects that new structures have been added to the database. It then re-triggers Module C, which **consumes** the new structures and **produces** DFT labels. These labels are then **consumed** by Module D, which **produces** a new, improved MLIP model. This new model is then consumed by Module E in the next iteration of the loop.

This circular data flow is the essence of the active learning process, allowing the system to continuously refine its own training data and improve its predictive accuracy.

## 4. Implementation Approach

The implementation will focus on building Module E, correctly handling boundary conditions, and then refactoring the orchestrator to manage the loop.

1.  **Configuration Update**:
    *   Add the `SimulationEngineConfig` and `ActiveLearningConfig` Pydantic models to `data/models.py`.
    *   Update the `ConfigExpander` to populate these new sections with sensible defaults in the `exec_config_dump.yaml`.

2.  **Boundary Utilities Implementation**:
    *   Create `src/mlip_autopipec/utils/boundary_utils.py`.
    *   Implement the `extract_and_passivate_cluster` function. This is a non-trivial task that requires careful geometric and chemical logic. It needs to correctly identify which bonds are "broken" by the cluster extraction and add terminating atoms (e.g., H) at those locations.
    *   Write extensive unit tests in `tests/unit/utils/test_boundary_utils.py` with various crystal structures to ensure the cluster extraction and passivation are working correctly and not creating physical artifacts.

3.  **Simulation Engine Implementation**:
    *   Create the `SimulationEngine` class in `src/mlip_autopipec/modules/e_simulation_engine.py`.
    *   The `execute` method will first load the latest MLIP `torch` model.
    *   It will set up an ASE MD simulation using the loaded model as the calculator.
    *   Implement the uncertainty quantification logic. This will involve inspecting the internals of the MACE model to get the disagreement between its different outputs, which serves as a proxy for uncertainty.
    *   Implement the main simulation loop. Inside the loop, periodically check the uncertainty. If it's high, call the `boundary_utils` to extract the structure and save it to the database.
    *   The method should return `True` if new structures were added, and `False` otherwise.
    *   Unit tests for this module will be complex. They will need to mock the MLIP model to return controllable uncertainty values, allowing the tests to trigger the structure extraction logic deterministically.

4.  **Orchestrator Refactoring**:
    *   Overhaul the `Orchestrator` class.
    *   The main `run` method will now implement the master loop as described in the architecture section. It will first call the initial A->B->C->D pipeline. Then, it will enter a `for` loop for the active learning generations, calling Module E, then C, then D in each iteration, based on the return value of Module E.
    *   This logic needs to be carefully managed to ensure the correct model is loaded and the correct data is processed in each iteration.

5.  **End-to-End Test Update**:
    *   Create a new E2E test, `test_cycle04_workflow.py`. This test will validate the active learning loop itself.
    *   It will require significant mocking. The MLIP model will need to be mocked to report high uncertainty for one specific frame in a simulation. The test will run the entire pipeline for one or two active learning iterations.
    *   The assertions will be critical:
        *   Verify that after the first simulation run, a new structure is added to the database for labeling.
        *   Verify that the `LabelingEngine` is re-run on this new structure.
        *   Verify that the `TrainingEngine` is re-run to produce a new model.
        *   Verify that the `SimulationEngine` is called again in the second iteration. This provides a complete validation of the closed-loop process.

## 5. Test Strategy

Testing in Cycle 4 is focused on validating the dynamic, iterative behavior of the active learning loop and the correctness of the complex structure extraction logic.

**Unit Testing Approach (Min 300 words):**
*   **`boundary_utils.py`**: The unit tests for this module must be rigorous to prevent the introduction of physical artifacts into the training data. In `test_boundary_utils.py`, we will create test cases for different kinds of materials (e.g., a covalent solid like silicon, an ionic crystal like NaCl). For each material, we will create a periodic ASE `Atoms` object. We will then call our `extract_and_passivate_cluster` function to extract a cluster around a central atom. The assertions will be multi-faceted:
    1.  Geometric: Verify that the extracted cluster contains the correct number of atoms within the specified cutoff radius.
    2.  Chemical: For covalent systems, verify that the correct number of passivating hydrogen atoms have been added and that their positions are chemically sensible (i.e., they have reasonable bond lengths).
    3.  Charge Neutrality: For ionic systems, future tests would need to verify that the extracted cluster is charge-neutral.
*   **`SimulationEngine`**: Testing the `SimulationEngine` requires controlling its behavior, which depends on the MLIP's uncertainty. In `test_e_simulation_engine.py`, we will create a mock MACE model. This mock model will be programmed to return a low uncertainty value for all inputs *except* for one specific, pre-defined atomic configuration, for which it will return a very high uncertainty. We will run the engine's `execute` method with a short MD simulation that is guaranteed to produce this specific configuration. We will then assert that the engine correctly detects the high uncertainty, calls the (mocked) `boundary_utils` function, and saves the new structure to the (mocked) database. This test deterministically validates the entire uncertainty detection and structure extraction workflow.

**Integration Testing Approach (Min 300 words):**
*   **Orchestrator Loop Logic**: We need to test the new cyclical logic of the `Orchestrator`. An integration test will be created where all modules (A, B, C, D, E) are replaced with mocks. These mocks will be simple objects that record whether their `execute` method was called. The mock for Module E will be configured to return `True` (found new structures) on its first call, and `False` on its second. The test will then run the orchestrator. The assertions will check the sequence of calls: A, B, C, D should be called once initially. Then, E, C, D should be called for the first active learning iteration. Finally, E should be called a second time, and because it returns `False`, C and D should *not* be called again. This test validates the entire control flow of the active learning loop without needing any real calculations.
*   **End-to-End Workflow Test**: `test_cycle04_workflow.py` will be the definitive test of this cycle's functionality. It will be a complex E2E test of the full loop, using the `CliRunner`. It will involve extensive mocking similar to the integration test described above but triggered from the top-level CLI command.
    1.  The initial pipeline (A->B->C->D) will run with mocks.
    2.  The `SimulationEngine` (E) will be called. Its internal MD simulation will be mocked, but its uncertainty logic will be real. The mock MLIP will be configured to flag one structure as uncertain.
    3.  The test will assert that a new structure has been added to the database.
    4.  The test will then verify that the `Orchestrator` correctly loops and re-runs the `LabelingEngine` (C) and `TrainingEngine` (D) on this new data point.
This test provides the highest level of confidence that the entire active learning machinery is correctly wired and functional from the user's perspective.
