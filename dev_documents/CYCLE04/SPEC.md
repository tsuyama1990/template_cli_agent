# Cycle 4 Specification: Active Learning & Advanced Simulation

## 1. Summary

Cycle 4 is the capstone of the MLIP-AutoPipe's automation strategy, transforming the system from a linear pipeline into a fully autonomous, self-improving "active learning" loop. The central component of this cycle is the **Simulation Engine (Module E)**. This module's purpose is to use the MLIP generated in the previous cycles to run large-scale molecular dynamics (MD) or kinetic Monte Carlo (kMC) simulations. Its most critical function, however, is not just to run simulations but to actively seek out the weaknesses of its own MLIP.

The core innovation of this cycle is the implementation of an "on-the-fly" (OTF) uncertainty quantification mechanism. As the simulation runs, Module E continuously monitors the MLIP's predictions. If it encounters an atomic configuration where the model is uncertain (i.e., it is extrapolating into unknown territory), it pauses the simulation, extracts the uncertain structure, and sends it back to the **Labeling Engine (Module C)** for DFT calculation. This new, high-value data point is then used to retrain and improve the MLIP. This feedback loop allows the system to iteratively refine its potential, patching its own deficiencies until it becomes robust and reliable across a wide range of configurations. This cycle also incorporates advanced simulation techniques like adaptive kMC to efficiently discover rare but important events, such as diffusion or chemical reactions, and focuses on optimising these simulation engines with Numba for maximum performance.

## 2. System Architecture

The architecture of Cycle 4 closes the loop, connecting the output of the pipeline (the MLIP) back to the input (the data generation stage).

**Module E: Simulation Engine (OTF & kMC)**
This module is the final piece of the automated workflow. It takes the latest trained MLIP from Module D and uses it to explore the material's behaviour over long timescales.

1.  **On-the-Fly (OTF) Active Learning Loop:**
    *   **Input:** A trained MLIP and a set of simulation parameters (temperature, pressure, duration).
    *   **Logic:** The engine initiates an MD simulation using a high-performance engine like LAMMPS. At each simulation step (or every few steps), a callback function is triggered.
    *   **Uncertainty Quantification:** This callback function calculates an "uncertainty score" for the current atomic configuration. This score quantifies how much the MLIP is extrapolating. The specific metric will depend on the MLIP model (e.g., the variance among predictions from a committee of models).
    *   **Thresholding:** The uncertainty score is compared against a pre-defined threshold from the configuration file.
    *   **Extraction:** If the score exceeds the threshold, the simulation is paused. The system then extracts the atomic configuration that caused the high uncertainty. To avoid artificial surface effects from cutting a structure out of a large simulation cell, it uses a "periodic embedding" technique, saving a small, fully periodic cell centered around the uncertain region.
    *   **Feedback:** This newly extracted structure is then sent back to the start of the pipeline to be labeled by Module C. The workflow then continues, with Module D retraining an improved MLIP, which is then passed back to Module E to resume the simulation. This cycle repeats until the simulation completes without triggering uncertainty, or a maximum number of iterations is reached.

2.  **Advanced Simulation (Adaptive kMC):**
    *   Beyond simple MD, Module E will include engines for exploring rare events. Adaptive kMC is a powerful technique for this.
    *   The core kMC loop (searching for saddle points, calculating transition rates, selecting an event) can be slow in pure Python. Key parts of this loop, especially the rate calculation and event selection, will be heavily optimised with Numba for performance.
    *   Uncertainty monitoring will also be active during these simulations, allowing the system to improve the accuracy of energy barriers for reactions and diffusion.

## 3. Design Architecture

The design focuses on creating a robust simulation engine that can be controlled by the main workflow orchestrator and can manage the active learning loop.

**`modules/e_simulation_engine.py`**:
-   **`SimulationEngine` Class**:
    -   **`__init__(self, config: FullConfig)`**: The constructor takes the full configuration, which includes settings for the OTF loop (uncertainty threshold, max iterations) and the simulation parameters.
    -   **`run_active_learning_loop(self, initial_mlip: Any) -> Any`**: The main public method. It takes the first trained MLIP as input. It contains the primary loop (`for generation in range(max_generations)`). In each iteration, it will:
        1.  Call `_run_simulation` with the current MLIP.
        2.  If the simulation returns a list of uncertain structures, it sends them to the labeling/training modules.
        3.  Receive a new, improved MLIP.
        4.  Repeat.
        The method returns the final, fully refined MLIP.
    -   **`_run_simulation(self, mlip: Any) -> list[ase.Atoms]`**: This method is responsible for a single "generation" of simulation. It will interface with an external code like LAMMPS. It prepares the input files for LAMMPS, including the MLIP, and starts the simulation. It will be configured to stop as soon as an uncertainty trigger occurs.
    -   **Uncertainty Callback:** The key mechanism will be a way for the running simulation to report back to the Python code. This could be done by having LAMMPS write a special "UNCERTAIN" file, which the `_run_simulation` method monitors for. Upon detection, it stops the LAMMPS process and parses the structure.
    -   **`_extract_periodic_embedding(self, full_structure: ase.Atoms, center_atom_index: int) -> ase.Atoms`**: This helper method will implement the logic for carving out a small, periodic sub-cell from a larger simulation cell.

**`common/akmc.py`**:
A new file to house the optimised kMC logic.
-   **`run_akmc_loop_numba(...)`**: A Numba-jitted function that contains the performance-critical parts of the adaptive kMC algorithm, such as calculating the rates of all possible events and selecting the next event using the Gillespie algorithm.

## 4. Implementation Approach

1.  **Uncertainty Metric Selection:** Choose a robust uncertainty metric compatible with the ACE model framework. A common approach is to train a committee of models and use the variance of their predictions as the uncertainty score.
2.  **LAMMPS Integration:** Develop the interface for running LAMMPS. This involves writing functions to:
    *   Generate a LAMMPS input script (`in.lammps`) from the `FullConfig`.
    *   Convert the trained MLIP into a format LAMMPS can read (a potential file).
    *   Run LAMMPS as a subprocess.
    *   Monitor the simulation for the uncertainty signal.
    *   Parse the final LAMMPS trajectory to extract the triggering structure.
3.  **Implement the OTF Loop:** Implement the main `run_active_learning_loop` in the `SimulationEngine`. This is the high-level orchestration logic.
4.  **Implement Structure Extraction:** Write the `_extract_periodic_embedding` helper function. This will require careful geometric manipulations to handle periodic boundary conditions correctly.
5.  **kMC Engine Implementation:** Implement the adaptive kMC logic in `common/akmc.py`. Profile the pure Python version to identify bottlenecks and then apply Numba optimisations to the critical loops.
6.  **Update Workflow Orchestrator:** The main orchestrator needs a significant update. It must be able to manage the new, cyclical data flow. Instead of a linear sequence, it will now manage a loop: `... -> D -> E -> (if uncertain) -> C -> D -> E -> ...`

## 5. Test Strategy

Testing Cycle 4 is centered on verifying that the active learning loop behaves as expected.

**Unit Testing Approach (Min 300 words):**
We will write unit tests for the individual components of the `SimulationEngine`. The uncertainty calculation will be tested by providing a mock model committee and asserting that the calculated variance is correct. For the `_extract_periodic_embedding` function, we will create a large `ase.Atoms` object and test the extraction of a sub-cell, asserting that the resulting cell is periodic and has the correct dimensions and atomic content. The Numba-optimised kMC functions will also have unit tests. We will run a simple kMC simulation with known rates and assert that the event sequence and timings match the expected analytical results. We will also have a performance test that asserts the Numba-jitted kMC loop is significantly faster than a pure Python equivalent.

**Integration Testing Approach (Min 300 words):**
The key integration test for this cycle is to prove that the OTF active learning loop works end-to-end. To do this reliably and quickly, we will create a scenario designed to fail. The test will proceed as follows:
1.  **Train a "bad" MLIP:** We will intentionally train an MLIP on a very small and narrow dataset (e.g., only a single, static crystal structure of silicon).
2.  **Start the Active Learning Loop:** We will run the `run_active_learning_loop` with this bad MLIP and configure it to run a short MD simulation at a high temperature.
3.  **Trigger Uncertainty:** The high temperature will cause the atoms to move into configurations that are very different from the single training structure. The MLIP's uncertainty score is therefore guaranteed to spike almost immediately.
4.  **Assertions:** The test will assert the following:
    *   The `_run_simulation` method stops prematurely and reports that an uncertain structure was found.
    *   The `run_active_learning_loop` correctly identifies this and sends the new structure to a *mocked* `LabelingEngine`. We use a mock here to avoid running a real DFT calculation in a test.
    *   The test verifies that the system attempts to retrain the model with the new data.
    *   The test asserts that one full cycle of the loop (simulate -> detect -> extract -> retrain) has been completed.
This test provides strong evidence that the entire feedback mechanism is correctly implemented.
