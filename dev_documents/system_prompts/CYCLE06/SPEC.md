# Specification: CYCLE06 - On-the-Fly Simulation & Active Learning

## 1. Summary

CYCLE06 introduces the dynamic, self-improving core of the MLIP-AutoPipe pipeline: `Module E: Simulation Engine` and the On-the-Fly (OTF) active learning loop. Until now, the workflow has been linear: generate data, then train a model. This cycle closes the loop. The primary goal is to use the MLIP trained in earlier stages to run a large-scale MD simulation and, crucially, to detect when the simulation enters an atomic configuration for which the model is uncertain. When uncertainty is detected, the system will automatically trap this "difficult" structure and send it back to the labelling queue, thus initiating a new cycle of learning and refinement.

The scope of this cycle is to implement the foundational OTF simulation capability. This involves integrating the trained MLIP with an MD engine (like LAMMPS, or initially, ASE's own dynamics), running a production-style simulation, and implementing an uncertainty quantification (UQ) metric. For ACE models, this UQ is often based on the variance of predictions from a committee of models or a dropout-based approach. The engine will monitor this uncertainty for every frame of the simulation. If the uncertainty exceeds a predefined threshold, the simulation is paused, the high-uncertainty structure is extracted, and it is added to the database for re-labelling by Quantum Espresso. The successful completion of this cycle transforms the pipeline from a static, one-shot generator into a dynamic, "living" system that actively seeks out its own weaknesses and systematically improves itself over time.

## 2. System Architecture

This cycle implements the `SimulationEngine`, a major new component that consumes the trained MLIP model and produces new candidate structures for labelling.

**File Structure (CYCLE06 Focus):**

Files to be created or modified are in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── **cli.py**              # Add a 'run-active-learning' command
│       ├── **workflow.py**         # Orchestrator manages the active learning loop
│       ├── **config.py**           # Add config for OTF simulation (threshold, etc.)
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── structure_generator.py
│           ├── explorer_sampler.py
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── **simulation_engine.py** # Module E - Main implementation
└── tests/
    ├── conftest.py
    ├── unit/
    │   └── **test_simulation_engine.py** # New test file
    └── integration/
        └── **test_active_learning_loop.py** # New integration test
```

**Component Breakdown:**

*   **`pyproject.toml`**: Dependencies for the production MD engine may be added, for example, if an interface to LAMMPS is used.
*   **`config.py`**: A new Pydantic model, `ActiveLearningConfig`, will be added.
    *   `ActiveLearningConfig`: Will contain fields like `md_engine: str` (default "ase"), `simulation_steps: int` (default 100,000), `uncertainty_threshold: float` (a fixed value for now), and `max_generations: int` (the number of times the active learning loop should run).
*   **`modules/simulation_engine.py`**: This file will house the `SimulationEngine` class.
    *   It will be initialized with the `ActiveLearningConfig`, the path to the trained MLIP model, and a reference to the `AseDBWrapper`.
    *   Its main public method, `run_otf_md()`, will be the entry point.
    *   Inside `run_otf_md`, it will load the trained MLIP and attach it as an ASE calculator. This calculator will be wrapped in a custom UQ calculator.
    *   The custom `UQCalculator` will override the default `get_potential_energy` or `get_forces` methods. Inside, it will call the underlying MLIP to get the forces and also compute an uncertainty score for the current atomic configuration.
    *   The MD simulation loop will be run step-by-step. In each step, after the forces are calculated, the uncertainty score is checked. If `uncertainty > threshold`, the loop is terminated, and the current `atoms` object is returned.
    *   If the loop completes without exceeding the threshold, it returns `None`.
*   **`workflow.py`**: The `WorkflowOrchestrator` will be significantly updated to manage the active learning loop.
    *   A new main method, `run_active_learning()`, will be created.
    *   This method will implement a `for` loop that runs for `max_generations`. Inside the loop, it will:
        1.  Train a new model based on all currently available labelled data (`TrainingEngine`).
        2.  Start a simulation with this new model (`SimulationEngine`).
        3.  If the simulation returns a high-uncertainty structure, add it to the database for labelling.
        4.  Label all new structures in the database (`LabelingEngine`).
        5.  Repeat.
*   **`cli.py`**: A new command, `cdd active-learn`, will be added to trigger the `run_active_learning` method in the orchestrator.

## 3. Design Architecture

The design of this cycle focuses on creating a robust and modular simulation loop with a clear mechanism for uncertainty feedback.

*   **`UQCalculator` Wrapper Design:**
    *   The core of the design is a custom ASE calculator that wraps the MLIP calculator. This is a powerful ASE pattern that allows for modular extension of functionality.
    *   `class UQCalculator(Calculator):`
        *   `__init__(self, mlip_calculator)`
        *   `calculate(self, atoms, properties, system_changes)`: This is the main method called by ASE.
            *   It calls `self.mlip_calculator.calculate(...)` to get the actual energy and forces.
            *   It then computes the uncertainty. For an ACE model, this might involve querying a committee of potentials and calculating the standard deviation of the predicted forces.
            *   It stores the computed uncertainty score, e.g., `self.results['uncertainty'] = score`.
    *   This design decouples the simulation logic from the UQ logic. The MD loop simply asks the calculator for its results, and the `UQCalculator` is responsible for providing both the physical properties and the uncertainty metric.

*   **Active Learning Loop Logic:**
    *   The orchestrator's `run_active_learning` method is the "brain" of the process. Its state management is critical.
    *   It will need to keep track of the current generation number.
    *   The loop's exit condition is either reaching `max_generations` or a state of convergence (e.g., a full simulation run completes without finding any uncertain structures).
    *   **Data Flow:** The loop creates a circular data flow:
        1.  The `TrainingEngine` reads from the DB.
        2.  The `SimulationEngine` uses the trained model and writes a new structure *back* to the DB.
        3.  The `LabelingEngine` reads this new structure from the DB and enriches it.
        4.  This enriched data is now available for the `TrainingEngine` in the next iteration.

## 4. Implementation Approach

The implementation will begin by creating the core `SimulationEngine` and its UQ logic, then wrapping it with the orchestrator's active learning loop.

1.  **Configuration:** Update `config.py` with the `ActiveLearningConfig` Pydantic model.
2.  **TDD for `UQCalculator`:** Create `test_simulation_engine.py`. The first test will focus on the `UQCalculator`.
    *   Create a mock MLIP calculator that returns a fixed energy and forces.
    *   The mock will also have a method like `get_force_variance()` that returns a predefined value.
    *   Wrap this mock in the `UQCalculator`. Run a single-point calculation (`atoms.get_forces()`).
    *   Assert that the returned forces match the mock's forces.
    *   Assert that you can access the uncertainty score from the calculator's results dictionary (`calc.results['uncertainty']`).
3.  **Implement `UQCalculator`:** Implement the `UQCalculator` to make the test pass.
4.  **TDD for `SimulationEngine`:** Write a test for the `run_otf_md` method.
    *   Use the tested `UQCalculator` (wrapping a mock MLIP).
    *   Configure the mock so that the uncertainty will be low for the first few steps, but exceed the threshold on the 5th step.
    *   Run `run_otf_md`.
    *   Assert that the method returns an `ase.Atoms` object.
    *   Check the positions of the returned atoms to ensure they correspond to the 5th step of the simulation, where the threshold was breached.
5.  **Implement `SimulationEngine`:** Implement the `run_otf_md` method, including the step-by-step MD loop and the uncertainty check, to make the test pass.
6.  **Implement Active Learning Loop:** In `workflow.py`, implement the main `run_active_learning` loop in the `WorkflowOrchestrator`. This will involve calling the already-existing `TrainingEngine` and `LabelingEngine` modules in a `for` loop, orchestrated around the new `SimulationEngine`.
7.  **New Integration Test:** Create a new test file, `test_active_learning_loop.py`. This will be a high-level test.
    *   It will mock all three engines: `TrainingEngine`, `LabelingEngine`, and `SimulationEngine`.
    *   The test will assert that the engines are called in the correct order within a loop. For example, in one loop iteration, it should assert: `train` -> `run_otf_md` -> `label`.
    *   The mock `SimulationEngine` will be configured to return a "high uncertainty" atom on the first call, and `None` on the second call, simulating convergence. The test will assert that the loop runs twice and then terminates.

## 5. Test Strategy

Testing this cycle is complex as it involves a stateful loop. The strategy is to unit test the components and then use a high-level integration test to verify the loop orchestration itself.

**Unit Testing Approach (Min 300 words):**
The unit tests will focus on the deterministic behavior of the `SimulationEngine` given a predictable (mocked) uncertainty signal.

*   **`test_simulation_engine.py`**:
    *   **`test_stops_on_uncertainty`**: As described above, this is the key test. It provides a mock MLIP where the uncertainty is programmed to spike at a specific timestep. The test asserts that the engine correctly stops at that exact frame and returns it. This validates the core trapping mechanism.
    *   **`test_runs_to_completion`**: This test provides a mock MLIP where the uncertainty *never* exceeds the threshold. The test asserts that `run_otf_md` runs for the full number of configured steps and returns `None`, indicating a successful and stable simulation run.
    *   **`test_calculator_setup`**: This test will verify that the `run_otf_md` method correctly loads the specified MLIP model file and wraps it in the `UQCalculator` before attaching it to the `ase.Atoms` object. This is tested by mocking the model loading function and the `UQCalculator`'s `__init__` method and asserting they were called with the correct arguments.

**Integration Testing Approach (Min 300 words):**
The integration test in `test_active_learning_loop.py` will verify the multi-generation orchestration without running any real science.

*   **`test_two_generation_loop`**: This test will simulate a two-generation active learning run.
    1.  **Setup**: An empty database and mocks for all three engines.
    2.  **Mock Behavior**:
        *   `TrainingEngine.train`: Mocked to do nothing.
        *   `LabelingEngine.label`: Mocked to do nothing.
        *   `SimulationEngine.run_otf_md`: This mock is crucial. We use `mocker.patch('...SimulationEngine.run_otf_md', side_effect=[high_uncertainty_atoms, None])`. This means the first time it's called, it will return a fake `ase.Atoms` object. The second time, it will return `None`.
    3.  **Execution**: Call `orchestrator.run_active_learning()` with `max_generations=5`.
    4.  **Assertions**:
        *   Check the call counts of the mocks. `TrainingEngine.train` should be called twice. `SimulationEngine.run_otf_md` should be called twice. `LabelingEngine.label` should be called once (only after the first simulation found something to label).
        *   We will also inspect the mock `AseDBWrapper`. We'll assert that its `add_atoms` method was called exactly once, with the `high_uncertainty_atoms` object that was returned by the first simulation run.
        *   The test will assert that the loop terminated after the second generation (when the simulation returned `None`), even though `max_generations` was 5. This proves the convergence logic is working correctly.
