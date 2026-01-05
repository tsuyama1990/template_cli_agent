# Cycle 2 Specification: Foundational Exploration Engine

## 1. Summary

Cycle 2 builds directly upon the foundational pipeline established in Cycle 1. The primary objective of this cycle is to introduce the core scientific computing capability of the project: a foundational Molecular Dynamics (MD) exploration engine. This transforms the tool from a simple structure generator into a powerful data generation pipeline that can explore the potential energy surface of a material.

The key deliverables for this cycle are:
- The implementation of an `ExplorationEngine` capable of managing multiple, concurrent MD simulations using Python's `ProcessPoolExecutor`. This is critical for achieving high throughput.
- A basic but functional MD runner (`md_mc_explorer.py`) that can take an `ase.Atoms` object and perform a standard NVT (constant volume and temperature) simulation using an ASE-compatible calculator.
- The integration of a real Machine Learning Interatomic Potential (MLIP), specifically MACE, as the calculator for computing forces and energies during the MD simulation. This is the first step in making the tool genuinely useful for its intended purpose.
- Enhancements to the `AseDBWrapper` to support the new data generated during this stage. This includes methods for retrieving structures that need to be explored and for storing metadata about the completed simulation runs, such as the file path to the output trajectory and the simulation parameters used.

By the end of Cycle 2, the MLIP-AutoPipe will have a significantly more powerful workflow. A user will be able to generate initial seed structures (as in Cycle 1) and then automatically subject each of these structures to an MD simulation in parallel. The results—full atomic trajectories showing how the structures evolve over time—will be saved to disk, and the database will be updated to link the initial structures to their corresponding simulation outputs. This cycle delivers the computational heart of the application.

## 2. System Architecture

This cycle introduces several new components and modifies existing ones to support the exploration stage. The files to be created or modified are marked in bold.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py
├── **config.py**               # Pydantic models updated for exploration
├── database/
│   ├── __init__.py
│   └── **ase_db_wrapper.py**   # Updated for exploration metadata
├── domain/
│   ├── __init__.py
│   ├── models.py
│   └── interfaces.py
├── engines/
│   ├── __init__.py
│   ├── generation_engine.py
│   ├── **exploration_engine.py** # **New implementation**
│   ├── sampling_engine.py    # (stubbed)
│   └── storage_engine.py
├── **explorers/**
│   ├── **__init__.py**
│   └── **md_mc_explorer.py**     # **New implementation** for MD logic
├── generators/
│   ├── __init__.py
│   ├── base_generator.py
│   └── alloy_generator.py
├── utils/
│   └── __init__.py
└── **workflow_orchestrator.py**  # Updated to call the ExplorationEngine
```

**File Blueprints:**

*   **`config.py`**:
    *   The `MainConfig` Pydantic model will be extended to include a new section for exploration parameters.
    *   A new `ExplorationConfig` model will be created to hold settings such as the simulation temperature, duration (number of steps), time step, and the name of the MLIP model to be used.

*   **`database/ase_db_wrapper.py`**:
    *   The `AseDBWrapper` class will be significantly enhanced.
    *   A new method, `get_unexplored_structures()`, will be added. This method will query the database for structures that have been generated but do not yet have associated exploration metadata, returning them as a list of `ase.Atoms` objects.
    *   A new method, `update_with_exploration_results(structure_id, trajectory_path)`, will be added. This method will allow the `ExplorationEngine` to update the database record for a specific structure, adding key-value pairs that store the path to the output XYZ trajectory file and other run metadata.

*   **`explorers/md_mc_explorer.py`**:
    *   This is a new file containing the core MD logic.
    *   It will feature a primary function, `run_md_simulation(atoms, config)`, which takes a starting `ase.Atoms` object and a configuration object.
    *   Inside this function:
        1.  An MLIP calculator (e.g., MACE) will be instantiated. **Crucially, the calculator instantiation will happen inside this worker function**, not in the main process. This avoids issues with pickling large, complex model objects when sending them to worker processes.
        2.  The calculator will be attached to the `atoms` object.
        3.  An ASE dynamics object, `ase.md.langevin.Langevin`, will be used to run the NVT simulation.
        4.  An `ase.io.trajectory.Trajectory` logger will be attached to the dynamics to save the state of the system at regular intervals to a file.
        5.  The simulation will be run for the configured number of steps.
        6.  The function will return the path to the saved trajectory file.

*   **`engines/exploration_engine.py`**:
    *   This file will contain the `ExplorationEngine` class, which orchestrates the parallel simulations.
    *   Its `run()` method will:
        1.  Use the `AseDBWrapper` to fetch the list of structures to be simulated.
        2.  Set up a `concurrent.futures.ProcessPoolExecutor`.
        3.  For each structure, it will submit a job to the pool executor, calling the `run_md_simulation` function from `md_mc_explorer.py`.
        4.  As each job completes, it will collect the resulting trajectory file path.
        5.  It will then call `db_wrapper.update_with_exploration_results()` to record the successful completion of the simulation in the database.

*   **`workflow_orchestrator.py`**:
    *   The `WorkflowOrchestrator`'s `run()` method will be updated to call the new `ExplorationEngine` in the correct sequence, after the `GenerationEngine`.

## 3. Design Architecture

The design of the exploration engine prioritizes performance, robustness, and avoidance of common pitfalls in parallel programming, particularly those related to large, non-serializable objects like ML models.

**Pydantic Models (`config.py`):**

The configuration schema will be expanded to provide user control over the MD simulation.

```python
from pydantic import BaseModel, Field
from typing import Optional

# ... existing AlloyGeneratorConfig ...

class MLIPConfig(BaseModel):
    """Configuration for the MLIP calculator."""
    model_name: str = Field("MACE_small", description="Name of the MLIP model to use.")
    # Add other model-specific args here in the future

class ExplorationConfig(BaseModel):
    """Configuration for the Exploration stage."""
    temperature_k: float = Field(300.0, gt=0, description="Simulation temperature in Kelvin.")
    num_steps: int = Field(1000, gt=0, description="Total number of MD steps.")
    time_step_fs: float = Field(1.0, gt=0, description="Time step for MD integration in femtoseconds.")
    mlip: MLIPConfig

class MainConfig(BaseModel):
    """Main configuration model updated for Cycle 2."""
    generator: AlloyGeneratorConfig
    exploration: Optional[ExplorationConfig] = None # Make exploration optional for now
    db_path: str = "structures.db"
```

**Design Principles:**

*   **Process Safety and Late Binding:** The most critical design choice is the "late binding" of the MLIP calculator. The main process, which manages the `ProcessPoolExecutor`, will *never* instantiate the MACE model. The model is large and complex and often cannot be "pickled" (serialized) to be sent to child processes. Instead, the `run_md_simulation` function, which executes *inside* the worker process, is responsible for creating its own instance of the calculator. This completely bypasses serialization issues and is a robust pattern for parallelizing ML-driven simulations.
*   **State Management:** The database is the single source of truth for the state of the workflow. The `ExplorationEngine` is designed to be restartable. If the process is killed midway, on the next run the `get_unexplored_structures()` method will simply return the structures that were not yet processed, and the engine will pick up where it left off. This is achieved by separating the act of fetching tasks (unexplored structures) from the act of recording their completion.
*   **Decoupling:** The `ExplorationEngine` is decoupled from the actual MD logic. Its job is process management and database communication. The scientific logic is encapsulated entirely within the `explorers` module. This separation of concerns makes the code easier to test and maintain. For example, we can test the parallelization logic of the engine by providing it with a simple mock function that just sleeps for a second and returns a fake file path, without needing to run a real, computationally expensive MD simulation.

## 4. Implementation Approach

The implementation will follow a logical progression, starting with the core MD logic and then building the parallel orchestration layer around it.

1.  **Update Configuration:**
    *   Modify `config.py` to include the new `ExplorationConfig` and `MLIPConfig` Pydantic models. Update `MainConfig` to include the exploration section.

2.  **Core MD Functionality (`md_mc_explorer.py`):**
    *   Create the `explorers` directory and the `md_mc_explorer.py` file.
    *   Implement the `run_md_simulation` function.
    *   Inside, add the logic to instantiate the MACE calculator using `ase.calculators.mace.MACE`. Note: this requires `mace-torch` to be added as a dependency.
    *   Set up the Langevin dynamics and attach the trajectory logger. The trajectory file should have a unique name, e.g., based on the structure's database ID (`md_run_{id}.traj`).
    *   Run the dynamics using `dyn.run(steps)`.

3.  **Enhance Database Wrapper (`ase_db_wrapper.py`):**
    *   Add the `get_unexplored_structures` method. This will involve using `db.select()` with a filter to find rows that do not have a `trajectory_path` key in their `key_value_pairs`.
    *   Add the `update_with_exploration_results` method. This will use the `db.update(id, key_value_pairs=...)` method to add the new metadata to a specific database row.

4.  **Implement the Engine (`exploration_engine.py`):**
    *   Create the `ExplorationEngine` class.
    *   The `run` method will first call the DB wrapper to get the work items (structures).
    *   It will then create a `ProcessPoolExecutor`.
    *   It will loop through the structures and use `executor.submit()` to schedule the `run_md_simulation` calls.
    *   It will iterate through the completed futures using `concurrent.futures.as_completed` to process results as they become available.
    *   For each result, it will call the DB wrapper's update method.

5.  **Integrate into Orchestrator:**
    *   Modify `WorkflowOrchestrator` to instantiate and run the `ExplorationEngine` after the `GenerationEngine`.

## 5. Test Strategy

Testing the exploration engine involves verifying both the scientific correctness of the MD simulation and the technical correctness of the parallel execution.

**Unit Testing Approach (Min 300 words):**

*   **`md_mc_explorer.py`:** Directly testing the `run_md_simulation` function is computationally expensive. We will instead test it with a very fast, simple, and deterministic calculator, like `ase.calculators.emt.EMT`, instead of MACE. The test will:
    1.  Create a simple `ase.Atoms` object.
    2.  Call `run_md_simulation` with a small number of steps (e.g., 10).
    3.  Assert that a trajectory file is created.
    4.  Read the trajectory file back using `ase.io.read`.
    5.  Assert that the number of frames in the trajectory is correct.
    6.  Check that the positions of the atoms in the final frame are different from the initial frame, confirming that the dynamics simulation actually ran.
*   **`ExplorationEngine`:** We will test the engine's parallelization and database logic without running real simulations. We will use `unittest.mock.patch` to replace the `run_md_simulation` function with a mock.
    1.  The test will set up a mock database containing a few "unexplored" structures.
    2.  The mock `run_md_simulation` will be configured to simply return a predictable fake file path for each structure ID it receives.
    3.  We will run the `ExplorationEngine`.
    4.  We will assert that the mock simulation function was called the correct number of times (once for each structure).
    5.  We will inspect the mock database object to assert that the `update_with_exploration_results` method was called with the correct structure IDs and the fake file paths returned by the mock simulation function. This isolates the test to the engine's orchestration logic.

**Integration Testing Approach (Min 300 words):**

The main integration test will verify the entire workflow from generation through exploration.
*   **Test Setup:**
    1.  The test will run in a temporary directory.
    2.  A Hydra YAML configuration will be created for a simple alloy.
    3.  The exploration section of the config will be enabled, but configured for a very short run: perhaps 2 structures, each for only 5-10 MD steps. This keeps the test runtime manageable. We will use a real but small MACE model to ensure the integration is working.
*   **Execution:**
    1.  The test will invoke the full pipeline via the CLI runner.
*   **Assertions:**
    1.  Assert the CLI command completes successfully.
    2.  Connect to the output database.
    3.  Verify that the database contains the 2 initial structures.
    4.  Check that the `key_value_pairs` for *each* of these structures now contains a `trajectory_path` key.
    5.  Verify that the file paths stored in the database actually exist on the filesystem.
    6.  Read one of the trajectory files and confirm it contains the correct number of frames. This end-to-end test provides high confidence that the configuration, database, parallel engine, and core MD logic are all correctly integrated.
