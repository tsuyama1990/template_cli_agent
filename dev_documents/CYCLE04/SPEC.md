# Cycle 04 Specification: Active Learning and Advanced Simulations

## 1. Summary

This document provides the technical specifications for Cycle 04 of the MLIP-AutoPipe project. This is arguably the most critical cycle in the project, as it transforms the pipeline from a linear, one-shot process into a dynamic, self-correcting, and truly autonomous system. The central goal of this cycle is to implement a closed-loop **Active Learning** workflow. This will be achieved by developing the **Simulation Engine** module, a sophisticated component that uses a trained MLIP to run large-scale simulations while simultaneously monitoring the model's own uncertainty about the atomic configurations it encounters. This cycle embodies the core philosophy of "removing the human expert from the loop" by automating the iterative process of model refinement, which is the most complex and intuition-driven part of creating a state-of-the-art potential.

The core of this cycle is the implementation of an "On-the-Fly" (OTF) inference loop. The `SimulationEngine` will take an MLIP trained in the previous stages and use it to run molecular dynamics (MD) or kinetic Monte Carlo (kMC) simulations. During these simulations, for every new atomic configuration generated, the engine will calculate an uncertainty score. This score quantifies the model's confidence in its own prediction. If this score exceeds a dynamically adjusted threshold, it serves as a signal that the model is operating in an "extrapolation domain"—a region of the chemical space that was not well represented in its original training data—where its predictions may be unreliable. When this triggering event occurs, the simulation is automatically paused, the high-uncertainty structure is carefully extracted, and it is sent back to the `LabellingEngine` for a high-fidelity DFT calculation. The newly labelled data point is then used to retrain and improve the MLIP, after which the simulation can be resumed with the newly refined, more knowledgeable model.

This cycle also specifies the implementation of advanced techniques that are crucial for the robustness and efficiency of the active learning loop. A **Dynamic Uncertainty Threshold** will be implemented, which allows the system's "curiosity" to adapt as the model matures, starting with a high tolerance for novelty and becoming progressively stricter. Furthermore, a sophisticated **Periodic Boundary Treatment** module will be developed to ensure that atomic cluster structures extracted from periodic simulations are physically meaningful and free from the artefacts that can poison the training data. The successful completion of this cycle will result in a truly autonomous pipeline that intelligently and efficiently identifies and patches the weaknesses in its own MLIP, ensuring the final potential is robust, accurate, and reliable across all the relevant regions of the material's phase space.

## 2. System Architecture

The architecture for Cycle 04 is a significant evolution of the pipeline. It "closes the loop" by connecting the output of the new `SimulationEngine` back to the input of the `LabellingEngine`, transforming the previously linear data flow into a cycle of continuous improvement.

1.  **Initial Training:** The workflow proceeds as in Cycle 03, starting from a minimal `input.yaml`. The `ConfigExpander`, `StructureGenerator`, `Explorer & Sampler`, `LabellingEngine`, and `TrainingEngine` all run in sequence. This results in an initial, "generation 0" version of the MLIP being trained. This initial model is reasonably accurate but may have significant gaps in its knowledge of the potential energy surface.
2.  **Module E: `SimulationEngine`:** This new module is invoked by the orchestrator after the initial training is complete. It is the core of the active learning loop.
    *   It loads the most recently trained MLIP from the file system.
    *   It initiates a long-running simulation (e.g., MD or kMC) using this potential. The type of simulation and its parameters (temperature, pressure, duration) are specified in the `active_learning` section of the `exec_config_dump.yaml` configuration file.
    *   **OTF Uncertainty Quantification:** At each step of the simulation, it calculates an uncertainty metric for the current atomic configuration. The specific method used will depend on the MLIP architecture but could include the variance in predictions from a committee of models, the distance of a configuration from the training data in a latent space, or other model-specific metrics.
    *   **Dynamic Threshold Check:** The calculated uncertainty score is compared against a dynamic threshold. This threshold is not a fixed, arbitrary value. Instead, at the beginning of each active learning generation, it is recalculated based on the distribution of uncertainties observed on the *existing* training set (e.g., as the 99th percentile of these values). This means as the model sees more data and becomes more confident, the threshold for what it considers "surprising" automatically becomes stricter.
    *   **Structure Extraction:** If the uncertainty score exceeds the current threshold, the simulation is paused. The current atomic structure is flagged. A sophisticated boundary treatment is applied here to handle periodicity correctly, for example, by extracting a non-periodic cluster around the atom with the highest local uncertainty, and passivating any dangling bonds that were created by cutting the cluster out of the bulk.
    *   **Feedback to Database:** The extracted, physically meaningful high-uncertainty structure is saved to the ASE database with the state "unlabelled", marking it as a new candidate for DFT calculation.
3.  **Active Learning Loop (Orchestrator's Role):** The main orchestrator is upgraded to manage the control flow for this new, cyclical workflow.
    *   After the `SimulationEngine` finishes a run (either by completing its allotted simulation time or by finding a new structure), the orchestrator checks the database for any new "unlabelled" structures.
    *   If new structures are found, the orchestrator re-invokes the **`LabellingEngine`** to perform DFT calculations on them.
    *   After they are labelled, the orchestrator re-invokes the **`TrainingEngine`** to retrain the MLIP with the newly augmented dataset. This creates a "generation 1" model.
    *   The orchestrator then calls the `SimulationEngine` again, which now loads the newly improved MLIP and can either restart the simulation or resume from where it left off.
    *   This loop continues—training, simulating, labelling, retraining—until a termination criterion specified in the configuration is met. This could be a maximum number of retraining iterations (e.g., 10 generations) or the simulation completing a full run without finding any more high-uncertainty configurations, at which point the model is considered to have converged.

## 3. Design Architecture

This cycle introduces the `SimulationEngine` class, a new utility module for boundary treatment, and significantly enhances the main orchestrator to manage the new cyclical workflow state.

-   **`src/mlip_autopipec/modules/simulation_engine.py`:**
    -   **`SimulationEngine` class:** This class encapsulates all the logic for running simulations and detecting uncertainty.
        -   `__init__(self, db: AseDB, config: dict)`: Initialises with the database and configuration dependencies.
        -   `run()`: The main public method. It loads the latest MLIP, sets up the simulation, and starts the OTF loop.
        -   `_load_latest_mlip(self)`: A private method to locate and load the most recent version of the trained potential file. It will need to be able to distinguish between different generations of the model.
        -   `_setup_simulation(self, atoms: ase.Atoms, mlip)`: Sets up the appropriate ASE dynamics object (e.g., `Langevin`, `NPT`) with the loaded MLIP as the calculator.
        -   `_run_otf_loop(self, dynamics)`: The core loop where the simulation runs step-by-step. Inside this loop, it attaches a callback function to the ASE dynamics object that is executed at every step. This callback function will trigger the uncertainty calculation and check it against the threshold.
        -   `_calculate_uncertainty(self, atoms: ase.Atoms) -> float`: Calculates the uncertainty score for the current structure. The exact implementation is a key scientific choice and will be designed to be pluggable, but an initial implementation might use force variance from a committee of models.
        -   `_handle_high_uncertainty(self, atoms: ase.Atoms)`: This method is called when the threshold is triggered. It will call the structure extraction and boundary treatment logic from the new utility module and then save the new candidate structure to the database.
-   **`src/mlip_autopipec/utils/boundary_treatment.py`:**
    -   This new module will contain pure functions for processing structures extracted from periodic simulations to make them suitable for retraining.
    -   `extract_and_passivate(atoms: ase.Atoms, center_atom_index: int) -> ase.Atoms`: An example function that will be implemented. It extracts a spherical cluster of a given radius around a central atom, identifies any bonds that were cut by the periodic boundary, and passivates the resulting dangling bonds (e.g., with hydrogen atoms). This is crucial for creating physically and chemically meaningful cluster models for the DFT retraining calculation.
-   **`src/mlip_autopipec/orchestrator.py`:**
    -   The main orchestrator logic will be significantly refactored from a linear script into a stateful class or a loop-based function.
    -   It will now manage a `generation` counter and the main active learning loop.
    -   The main loop will be clearly structured:
        ```python
        for generation in range(max_generations):
            training_engine.run()  # Train or retrain the model
            simulation_engine.run() # Run simulation until a new structure is found or time is up

            # Check for convergence
            new_structures = db.get_atoms_by_state('unlabelled')
            if not new_structures:
                print("Convergence reached.")
                break

            labelling_engine.run() # Label the new structures
        ```
-   **`src/mlip_autopipec/config/expander.py`:**
    -   The `ConfigExpander` will be updated to populate a new `active_learning` section in the `exec_config_dump.yaml`. This will include parameters that the user can optionally control, such as `max_generations`, the initial uncertainty threshold, the choice of uncertainty metric, and the type of simulation to run in the `SimulationEngine`.

## 4. Implementation Approach

The implementation of the active learning loop will be staged to manage its complexity.

1.  **Basic Simulation Setup:** The first step is to implement the basic simulation setup in the `SimulationEngine`. This involves loading the ACE model trained in Cycles 01-03 and using it to run a standard, non-OTF ASE MD simulation. The initial goal is simply to verify that the trained potential is stable enough to run dynamics for an extended period.
2.  **Uncertainty Metric:** A first, simple version of the `_calculate_uncertainty` method will be implemented. The choice of metric is key. If the MLIP model natively supports an uncertainty metric (like MACE's latent space), that will be used. If not, a simpler heuristic (e.g., based on the magnitude of atomic forces or local energy contributions) can be used as a starting point, with the design allowing for more complex metrics to be added later.
3.  **Orchestrator Refactoring:** The main orchestrator script will be refactored into a loop that can manage the `generation` state. This is a significant structural change. Initially, the loop can be tested by having the `SimulationEngine` always produce a "fake" high-uncertainty structure after a few steps. This will allow for the verification that the orchestrator correctly calls the labelling and training engines in response, validating the control flow before the real uncertainty metric is integrated.
4.  **Dynamic Threshold Logic:** The logic for the dynamic uncertainty threshold will be implemented. This requires a new utility function that can take the current training set from the database, calculate the uncertainty for all structures within it using the chosen metric, and then use NumPy to determine a high percentile value (e.g., the 99th percentile). The `SimulationEngine` will call this function at the start of each new generation to get an updated threshold.
5.  **Boundary Treatment:** The `boundary_treatment` utility module will be implemented. This is a complex but crucial step. The implementation will start with a simple "buffer" approach, where a larger-than-needed cluster is extracted, and only the forces on the central atoms are used for retraining. The more advanced passivation logic, which is chemically complex, will be added and tested subsequently.
6.  **Full Loop Integration:** Finally, all the pieces will be connected. The `SimulationEngine` will be updated to call the real uncertainty metric, compare it to the real dynamic threshold, and on triggering, call the new boundary treatment utility before saving the structure to the database. The orchestrator will be set to run the full, real loop. This final integration will require careful end-to-end debugging of the entire autonomous workflow.

## 5. Test Strategy

Testing Cycle 04 is challenging due to the stochastic and complex nature of the active learning loop. The strategy will rely on a combination of unit tests for the pure-logic components and carefully designed integration tests that verify the control flow of the loop with mocked, deterministic behaviours.

**Unit Testing Approach (Min 300 words):**
-   **`SimulationEngine`:** The `SimulationEngine` itself will be tested by mocking the ASE dynamics object. Instead of running a real, time-consuming MD simulation, the mock object will be configured to `yield` a pre-defined sequence of `ase.Atoms` objects, simulating the steps of a trajectory. A test will verify that the engine correctly calls the (mocked) `_calculate_uncertainty` method for each of these structures. Another, more complex test will involve configuring the mock uncertainty calculator to return a value that exceeds the threshold on, for example, the 5th step. The test will then assert that the `_handle_high_uncertainty` method is called correctly at this point, and that the database's `add_atoms` method is subsequently called with the correct structure.
-   **Dynamic Threshold Calculator:** The function that calculates the dynamic threshold will be unit-tested as a pure function. It will be given a list of mock `Atoms` objects with pre-defined uncertainty scores attached, and the test will assert that it correctly calculates the specified percentile using `numpy.percentile`, verifying the statistical logic.
-   **`boundary_treatment`:** The functions in this module are also pure functions and are well-suited for unit testing. For example, the `extract_and_passivate` function will be given a periodic `ase.Atoms` object with a known, simple structure (e.g., bulk silicon). The test will assert that the returned cluster has the correct number of atoms, that its `pbc` flags are all `False`, and that new hydrogen atoms have been added at the correct locations to passivate any dangling bonds that were created by the virtual "cut". This validates the complex geometric and chemical logic of this crucial utility.

**Integration Testing Approach (Min 300 words):**
-   **Simulation -> DB -> Labelling:** A key integration test will verify the feedback loop. It will start with a trained model and run the `SimulationEngine`. The uncertainty model will be mocked to trigger deterministically on the 5th step of a short MD run. The test will then assert that a new "unlabelled" structure appears in the (real, temporary) `AseDB`. The test will then instantiate and run the `LabellingEngine` (with a mocked DFT call that returns immediately) and will verify that this new structure's state is correctly updated to "labelled". This test is critical as it validates that the data flows correctly from the end of the pipeline (simulation) back to the beginning (labelling).
-   **Full Active Learning Loop "Dry Run":** A comprehensive end-to-end test will be created to validate the entire orchestration of the active learning loop. It will use a simple, fast-to-train model. All external processes (DFT) and long-running simulations (MD) will be heavily mocked or configured to run for a trivial number of steps (e.g., MD runs for 10 steps, the mock uncertainty calculator triggers on step 5). The test will configure the `max_generations` parameter to 2. It will then run the main orchestrator and assert that the sequence of high-level operations is correct: `TRAIN-gen-0` -> `SIMULATE-gen-0` (finds a new structure) -> `LABEL-gen-1` -> `TRAIN-gen-1` -> `SIMULATE-gen-1` (finds no new structures) -> `TERMINATE`. This test does not check for physical correctness but provides high confidence that the complex control flow of the entire autonomous loop is implemented correctly.