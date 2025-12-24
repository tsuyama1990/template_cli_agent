# CYCLE04 Specification: Active Learning and Simulation Integration

## 1. Summary

CYCLE04 represents the culmination of the project's core vision: creating a fully autonomous, self-improving system. This cycle introduces the **Simulation Engine (Module E)** and the "On-the-fly" (OTF) active learning loop. While previous cycles built a powerful pipeline for generating an initial MLIP, this cycle gives that MLIP the ability to discover its own weaknesses and iteratively improve itself without any human intervention.

The primary objective is to close the feedback loop between the trained model and the data generation process. The `SimulationEngine` will use the newly trained MLIP to run realistic, long-time-scale simulations (e.g., MD or Kinetic Monte Carlo). Crucially, during this simulation, it will continuously monitor the model's prediction uncertainty. If the simulation encounters an atomic configuration that the model is not confident about (an "out-of-sample" or extrapolation region), the engine will automatically pause, extract that challenging structure, send it back to the DFT Labeling Engine (Module C) for accurate calculation, and then retrain an improved version of the MLIP using this new data point.

This active learning loop is the system's most powerful feature, transforming it from a static "train-once" tool into a dynamic, "always-learning" platform. This cycle will also involve the integration of a robust simulation code like **LAMMPS**, which is necessary for running the large-scale, long-duration simulations where rare events and model failures are most likely to be discovered. Furthermore, to enable the study of phenomena like diffusion and phase transitions, this cycle will include the implementation of an Adaptive Kinetic Monte Carlo (kMC) routine. The performance-critical components of this kMC loop will be accelerated with **`Numba`** to ensure the system can efficiently simulate events over microsecond timescales and beyond.

## 2. System Architecture

The architecture in CYCLE04 evolves from a linear pipeline into a true loop. The new `SimulationEngine` (Module E) acts as both the final consumer of the MLIP and the initiator of a new data generation cycle.

The complete, looped data flow is now:
1.  **Initial Model Generation (from Cycles 1-3)**: The workflow begins by generating the first version of the MLIP (`MLIP v1`) using the full pipeline developed previously (Structure Generation -> Surrogate Exploration -> DFT Labeling -> Training).
2.  **Simulation (Module E - New)**: The `SimulationEngine` takes `MLIP v1` and begins a large-scale simulation (e.g., an MD run at a high temperature). The simulation is executed by an external code like LAMMPS, controlled by the engine.
3.  **Uncertainty Monitoring (Module E - New)**: At every step of the simulation, the engine calculates an "uncertainty metric". This could be the variance between predictions from a committee of models, or a specific extrapolation grade provided by the model itself (e.g., in the ACE framework).
4.  **Trigger and Pause**: If the uncertainty metric exceeds a predefined threshold, it signals that the model is unreliable for the current atomic structure. The engine immediately pauses the LAMMPS simulation.
5.  **Structure Extraction (Module E - New)**: The engine extracts the problematic structure from the simulation. To avoid surface effects from simple cluster cutouts, it will use a "Periodic Embedding" technique, extracting a small, periodic simulation cell centered on the high-uncertainty region.
6.  **Re-Labeling (Module C)**: This extracted structure is sent back to the `LabelingEngine`. Module C calculates its DFT properties just like any other structure.
7.  **Re-Training (Module D)**: The new, high-value data point is added to the training database. The `TrainingEngine` is called again to produce an improved model, `MLIP v2`.
8.  **Resume Simulation**: The `SimulationEngine` replaces `MLIP v1` with the new, more robust `MLIP v2` and resumes the LAMMPS simulation from the point where it was paused.

This loop (Simulate -> Detect Uncertainty -> Pause -> Extract -> Re-Label -> Re-Train -> Resume) continues until the simulation can proceed for very long times without encountering any more high-uncertainty configurations, indicating that the potential has become robust for the conditions being simulated.

## 3. Design Architecture

This cycle's main development effort is in `e_simulation_engine.py` and its integration with LAMMPS. It also requires a robust mechanism for communication and control between the Python orchestrator and the external simulation binary.

**Key Files and Classes:**

*   `src/mlip_autoprope/modules/e_simulation_engine.py`:
    *   `SimulationEngine` class: The main orchestrator.
    *   `run(mlip_model, config)`: The main public method.
    *   `_prepare_lammps_input(mlip_path, sim_params)`: A method to generate the LAMMPS input script, including commands to load the MLIP, set up the simulation box, and define the MD or kMC run.
    *   `_run_lammps()`: The method that executes the LAMMPS binary as a subprocess. It will need to monitor the output of LAMMPS in real-time to detect signals for pausing.
    *   `_monitor_uncertainty(lammps_output)`: This function will parse the live output from a specially configured LAMMPS run to extract the uncertainty values at each step.
    *   `_extract_periodic_embedding(structure_file)`: A method to read a LAMMPS dump file and carve out the small, periodic cell needed for re-labeling.
*   `src/mlip_autoprope/modules/d_training_engine.py`: The `TrainingEngine` will be modified to support incremental training. When a new data point is added, it should be able to load the existing MLIP and fine-tune it with the new data, rather than re-training from scratch every time.
*   `src/mlip_autoprope/modules/common/uncertainty.py`: A new file to house the logic for calculating uncertainty. This might include wrappers for specific model features (like ACE's `get_extrapolation_grade`) or committee model logic.
*   `src/mlip_autoprope/modules/common/akmc.py`:
    *   `AdaptiveKMC` class: An implementation of the adaptive kMC algorithm.
    *   `run_step()`: A `Numba`-jitted function to perform a single kMC step (calculating rates, selecting an event, advancing time). This is a performance-critical loop.
*   `src/mlip_autoprope/config/core.py`: The `ExecConfig` will be updated with a `simulation` section to control the active learning loop. Parameters will include the uncertainty threshold, the type of simulation (MD, kMC), and settings for the simulation itself (temperature, pressure, duration).

**New Dependencies in `pyproject.toml`:**
*   A dependency on a specific LAMMPS installation will not be in `pyproject.toml`, but the system will be configured to expect a `lammps` executable in the PATH. The project build process will now include instructions or a script to download and compile a compatible version of LAMMPS.

## 4. Implementation Approach

1.  **LAMMPS Integration**: The first step is to establish basic control over LAMMPS. Implement `_prepare_lammps_input` and `_run_lammps` to the point where the Python code can launch a simple MD simulation of a known system (e.g., silicon) using the MLIP trained in previous cycles.
2.  **Uncertainty Communication**: Modify the LAMMPS input script and the MLIP pair style to print per-atom uncertainty metrics to the standard output at each step. Implement the `_monitor_uncertainty` function in Python to parse this live output from the running LAMMPS subprocess.
3.  **Pause/Resume Mechanism**: Implement the control loop. The Python orchestrator needs to be able to send a `SIGSTOP` or similar signal to the LAMMPS process when high uncertainty is detected, and `SIGCONT` to resume it later. The state (atom positions, velocities) must be saved by LAMMPS before pausing so it can be reloaded.
4.  **Structure Extraction**: Implement the `_extract_periodic_embedding` function. This involves reading a LAMMPS dump file, identifying the region of interest, and using tools (like those in ASE) to wrap the atoms into a new periodic cell.
5.  **Loop Integration**: Connect all the pieces. The `SimulationEngine.run()` method will orchestrate the full loop: call LAMMPS -> monitor output -> pause on trigger -> extract structure -> pass to Module C -> re-train with Module D -> resume LAMMPS with the new potential.
6.  **kMC Implementation**:
    *   Implement the core kMC logic in pure Python first to ensure correctness. This involves saddle point searches (e.g., using NEB or Dimer method, possibly via ASE) and rate calculations (Harmonic Transition State Theory).
    *   Once the logic is correct, apply `Numba`'s `@jit` decorator to the main event loop and rate calculation functions.
    *   Benchmark the `Numba`-optimised kMC implementation to verify a significant performance increase.
7.  **Incremental Training**: Refactor the `TrainingEngine` to support loading an existing model and fine-tuning it with a small number of new data points, which is much faster than retraining from scratch.

## 5. Test Strategy

Testing for CYCLE04 focuses on the robustness of the active learning loop and the correctness of the simulation and kMC components.

### Unit Testing Approach

*   **LAMMPS Input Generation**: Write a unit test that calls `_prepare_lammps_input` with a given configuration and asserts that the resulting script string matches a known-good LAMMPS input file.
*   **Uncertainty Parsing**: Create a text file that mimics the live output of a LAMMPS run, including lines with uncertainty data. Write a unit test that feeds this file to `_monitor_uncertainty` and asserts that it correctly parses and returns the numerical uncertainty values.
*   **Periodic Embedding**: Create a sample dump file of a large cell with atoms near the boundaries. Write a unit test that calls `_extract_periodic_embedding` and asserts that the resulting ASE `Atoms` object is periodic and contains the correct atoms, properly wrapped.
*   **kMC Logic**: Create a simple 2D potential energy surface with known minima and saddle points. Unit test the kMC algorithm's ability to correctly identify the escape paths and calculate the transition rates. Compare the pure Python and `Numba` versions for numerical consistency.

### Integration Testing Approach

*   **Pause/Resume Test**: This is a critical integration test. It will run the `SimulationEngine` on a simple system. The test will mock the uncertainty monitor to artificially trigger a pause signal after a few hundred steps. The test will verify that:
    1.  The LAMMPS process is successfully paused.
    2.  The engine correctly simulates the "re-training" step.
    3.  The LAMMPS process is successfully resumed and continues the simulation.
    4.  The final trajectory file is continuous and does not show a disjoint jump.
*   **End-to-End Active Learning Loop Test**: This is the main UAT for the cycle. It will use a deliberately "bad" initial MLIP (e.g., trained on only low-temperature data).
    1.  Start a high-temperature MD simulation using this bad MLIP.
    2.  Assert that the `SimulationEngine` quickly detects high uncertainty as the system explores new configurations.
    3.  Verify that the system automatically pauses, extracts a structure, runs a DFT calculation, and retrains the model.
    4.  Verify that the database contains a new entry corresponding to the structure extracted from the simulation.
    5.  Verify that a new, improved MLIP file (`model_v2.yace`) is created.
    6.  The test will assert that the simulation is able to run for a longer period with the new model before triggering another uncertainty event, proving that the system is learning and improving.
