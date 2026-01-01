# CYCLE 04: SPEC.md - Active Learning and On-the-Fly Simulation

## 1. Summary

Cycle 04 represents the culmination of the MLIP-AutoPipe project's core vision, transforming the linear pipeline into a dynamic, self-improving, autonomous system. The central objective of this cycle is to implement the active learning loop via `Module E: Simulation Engine`. This module will deploy the MLIP trained in the earlier stages to run large-scale simulations (e.g., Molecular Dynamics). The key innovation is the "on-the-fly" (OTF) uncertainty quantification. During the simulation, the engine will continuously monitor the reliability of the MLIP's predictions. When the system encounters an atomic configuration that is significantly different from its training data—an "out-of-distribution" or extrapolation scenario—the model's uncertainty will spike.

The `SimulationEngine` will automatically detect these high-uncertainty events, pause the simulation, and extract the local atomic environment that caused the trigger. To ensure this extracted cluster is physically valid for a high-precision DFT calculation, the system will apply sophisticated boundary treatments, such as passivating dangling bonds. This new, information-rich structure is then fed back into the `LabelingEngine` (Cycle 01) to be labeled with DFT. Finally, the `TrainingEngine` retrains the MLIP with this new data point, improving its accuracy and expanding its domain of applicability. This entire process—simulate, detect, extract, label, retrain—forms a closed loop that runs autonomously until the MLIP is robust enough to handle the desired simulation conditions without encountering new, unknown configurations. This cycle elevates the system from a simple automation tool to an intelligent agent for scientific discovery.

## 2. System Architecture

This cycle introduces the final core module, `e_simulation_engine.py`. The most significant architectural change is the refactoring of the `WorkflowOrchestrator` from a simple sequential runner into a stateful loop manager that can repeatedly cycle through the simulation, labeling, and training stages.

**File Structure (Cycle 04 Focus):**

```
.
├── src/
│   └── mlip_autopipec/
│       ├── modules/
│       │   ├── ...
│       │   └── **e_simulation_engine.py** # NEW module implementation
│       ├── **utils/**
│       │   └── **cluster_utils.py**   # NEW helper for cluster extraction/passivation
│       └── **workflow.py**            # HEAVILY MODIFIED for active learning loop
└── pyproject.toml
```

The files marked in **bold** are the primary deliverables. `e_simulation_engine.py` is the new module responsible for running simulations and detecting uncertainty. A new utility, `cluster_utils.py`, is introduced to handle the complex and critical logic of extracting and treating the boundary conditions of atomic clusters. The `workflow.py` file will undergo a major redesign to manage the new, non-linear, iterative workflow.

## 3. Design Architecture

The design of this cycle revolves around creating a robust feedback loop, which requires careful state management and a clear definition of the loop's stopping conditions.

**Pydantic Schema (`configs/models.py`):**

The configuration will be updated to control the active learning loop.

*   `ActiveLearningConfig` (to be added to `MainConfig`):
    *   `max_generations: int` (The maximum number of retraining iterations)
    *   `uncertainty_threshold: Union[float, str]` (The trigger for retraining. A float for a fixed value, or a string like "dynamic_95percentile" for an adaptive threshold).
    *   `md_steps_per_generation: int` (The number of MD steps to run in each simulation phase).

**Class and Module Design:**

*   **`SimulationEngine` (`modules/e_simulation_engine.py`):**
    *   Constructor accepts `ActiveLearningConfig`, `AseDBWrapper`, and the path to the current MLIP model file.
    *   Main public method, `run()`, executes one generation of the simulation.
    *   **Internal Logic:**
        1.  **Load MLIP:** Load the trained potential and attach it as an ASE calculator.
        2.  **Run MD:** Run an MD simulation for `md_steps_per_generation`.
        3.  **Monitor Uncertainty:** At each MD step, it will query the calculator for an uncertainty metric. (This assumes the MLIP model/calculator, like MACE, exposes this.)
        4.  **Trap & Extract:** If uncertainty > `uncertainty_threshold`, pause the simulation. It will then call a helper function in `cluster_utils.py` to extract the local atomic environment around the atom(s) with the highest uncertainty.
        5.  **Save for Labeling:** The newly extracted cluster is passed to the `AseDBWrapper` to be saved to the database with a unique status, e.g., 'unlabeled_active'. The engine can be configured to find one or multiple uncertain structures per generation.

*   **`cluster_utils.py` (`utils/cluster_utils.py`):**
    *   `extract_cluster(full_atoms, center_indices)`: A function that takes a large `Atoms` object and the indices of the high-uncertainty atoms, and returns a new `Atoms` object representing the local cluster.
    *   `passivate_cluster(cluster_atoms)`: A function that takes a raw cluster, identifies dangling bonds at its surface, and adds passivating atoms (e.g., Hydrogen) to create a chemically saturated, physically realistic structure suitable for DFT.

*   **`WorkflowOrchestrator` (`workflow.py`):**
    *   The orchestrator will be refactored to manage the main loop.
    *   The `run()` method will now look something like this:
        ```python
        def run(self):
            # Initial steps (Config, Structure Gen, Initial Sampling)
            ...
            # Initial Labeling and Training
            self.labeling_engine.run()
            self.training_engine.run()

            for generation in range(self.config.active_learning.max_generations):
                self.logger.info(f"Starting active learning generation {generation + 1}")

                # 1. Simulate and find new structures
                self.simulation_engine.run() # This populates the DB

                # 2. Check if we need to continue
                new_structures = self.db.get_rows(status='unlabeled_active')
                if not new_structures:
                    self.logger.info("No new uncertain structures found. Converged.")
                    break

                # 3. Label the new structures
                self.labeling_engine.run(label_status='unlabeled_active')

                # 4. Retrain the model with all labeled data
                self.training_engine.run()
        ```

## 4. Implementation Approach

The implementation will focus on building the simulation engine and then refactoring the orchestrator to manage the loop.

1.  **Update Configuration:** Add the `ActiveLearningConfig` to the Pydantic models and update the `ConfigExpander`.
2.  **Implement Cluster Utilities:** Begin with the scientifically challenging part in `cluster_utils.py`. Implement robust logic for extracting a cluster and, crucially, for passivating its surface.
3.  **Implement SimulationEngine:** Create the `SimulationEngine` class.
    *   Integrate the MLIP model as an ASE calculator.
    *   Write the MD simulation loop.
    *   Implement the core uncertainty-trapping logic. This will require careful study of the chosen MLIP's API to get the uncertainty values.
    *   Integrate the functions from `cluster_utils.py` to process the trapped structures before saving them to the database.
4.  **Refactor LabelingEngine:** Modify the `LabelingEngine`'s `run` method to optionally accept a status parameter, so it can be directed to label either the initial structures or the new 'unlabeled_active' structures.
5.  **Refactor WorkflowOrchestrator:** This is a major step. Rewrite the `run` method from a linear script into the stateful loop described in the design section. This involves careful management of the state between generations.
6.  **Integrate LAMMPS (Optional):** If ASE's built-in MD is not performant enough, this step would involve creating a wrapper to use LAMMPS as the MD engine, which is more complex but offers higher performance. For this cycle, starting with ASE's engine is sufficient.

## 5. Test Strategy

Testing the active learning loop requires careful mocking to create a reproducible test environment.

**Unit Testing Approach:**
(Located in `tests/unit/`)
*   **`utils/cluster_utils.py`:** Write tests for the passivation logic. Create an artificial, unpassivated cluster of a covalent material (like a chunk of silicon) and assert that the `passivate_cluster` function correctly adds hydrogen atoms to the surface atoms with dangling bonds.
*   **`modules/e_simulation_engine.py`:** This is the most critical unit test of the cycle.
    *   **Setup:** Create a mock MLIP calculator object. This mock will be programmed to return a high uncertainty value only when it sees a very specific atomic arrangement.
    - **Execution:** Run the `simulation_engine.run()` method with an MD trajectory that is guaranteed to contain this specific arrangement.
    *   **Verification:**
        1.  Assert that the engine correctly identifies the high-uncertainty frame.
        2.  Assert that the `extract_cluster` and `passivate_cluster` utilities are called.
        3.  Assert that the final, passivated `Atoms` object is passed to the mocked database wrapper with the 'unlabeled_active' status.

**Integration Testing Approach:**
(Located in `tests/e2e/`)
*   **Single Loop E2E Test:** This test will verify that the entire loop is connected and functions for one full iteration.
    *   **Setup:** Prepare a database that has already been through an initial training run and contains a trained MLIP model file.
    *   **Execution:** Run the main `WorkflowOrchestrator`. Use the same mock MLIP calculator from the unit test to ensure that the `SimulationEngine` will find an uncertain structure. Mock the DFT `subprocess.run` call and the final `TrainingEngine` call.
    *   **Verification:** The test will assert the following sequence of events:
        1.  The `SimulationEngine` runs and adds a new 'unlabeled_active' structure to the database.
        2.  The `LabelingEngine` runs and changes the status of this new structure to 'labeled'.
        3.  The `TrainingEngine` runs again, signifying that a retraining iteration was triggered.
        This confirms that the orchestrator's loop logic is correct and the modules communicate correctly through the database.
