# SPEC.md - CYCLE 02

## 1. Summary

This specification document outlines the technical requirements for Cycle 02 of the **MLIP-AutoPipe** project. Building upon the foundational framework established in Cycle 01, this cycle is dedicated to implementing the sophisticated core functionalities of the application: the **Exploration** and **Sampling** stages. The primary goal of this cycle is to transform the simple structure generator into a powerful tool capable of producing diverse, high-quality datasets for training Machine Learning Interatomic Potentials (MLIPs).

The key deliverable for this cycle is the implementation of a robust, parallelized exploration engine that uses Molecular Dynamics (MD) and Monte Carlo (MC) methods to discover a wide range of atomic configurations. This engine will be capable of using MLIP models as calculators via the Atomic Simulation Environment (ASE) and will feature advanced logic for automatically selecting the appropriate simulation ensemble. Following the exploration stage, we will implement an intelligent sampling module, featuring Farthest Point Sampling (FPS), to select a maximally diverse subset of structures from the vast number of configurations generated during exploration.

By the end of Cycle 02, the MLIP-AutoPipe tool will be feature-complete with respect to its core pipeline (Generate -> Explore -> Sample -> Store). A user will be able to define a complex workflow in a configuration file, including MD simulation parameters and sampling strategies, and execute it via the command-line interface to produce a refined, diverse, and ready-to-use database for MLIP training. This cycle represents the "brains" of the operation, turning a simple structure generator into an intelligent data-generation engine.

## 2. System Architecture

Cycle 02 involves creating new modules for the exploration and sampling stages and integrating them into the existing pipeline established in Cycle 01. The file structure will be expanded, and existing modules like `config.py`, `pipeline.py`, and `cli.py` will be modified to support the new functionality.

**File Structure (Cycle 02 Create/Modify):**

```
.
└── src
    └── mlip_autopipec
        ├── **cli.py**                # Modified to add new CLI options
        ├── **config.py**             # Modified to add exploration/sampling config
        ├── **exploration**
        │   ├── **__init__.py**
        │   └── **engine.py**         # **ExplorationEngine** implementation
        ├── **sampling**
        │   ├── **__init__.py**
        │   ├── **base.py**           # **BaseSampler** ABC
        │   ├── **random_sampler.py** # **RandomSampler** implementation
        │   └── **fps.py**            # **FarthestPointSampler** implementation
        ├── **pipeline.py**           # **Modified** to run the full 4-stage pipeline
        └── **services**
            └── **workflow.py**       # **Modified** to manage the full pipeline
```

**Component Blueprint:**

*   **`config.py`**: The Pydantic models will be significantly expanded.
    *   `ExplorationConfig`: A new model to hold all MD/MC parameters, including `temperature`, `pressure`, `timestep`, `n_steps`, the MLIP `model_path`, and settings for MC `swap_frequency`.
    *   `SamplingConfig`: A new model defining the sampling strategy. It will include a `method` field (e.g., a `Literal["random", "fps"]`), the number of structures to select, and parameters for FPS like the SOAP descriptor settings.
    *   `FullConfig`: The top-level model will be updated to include these new `ExplorationConfig` and `SamplingConfig` sections.

*   **`exploration/engine.py`**: This is the most critical and complex new module.
    *   **`ExplorationEngine` class**: The main class responsible for the exploration stage.
    *   `__init__`: Will accept the `ExplorationConfig` and the list of seed structures from the generation stage.
    *   `run()`: The main public method that orchestrates the parallel simulations. It will use `concurrent.futures.ProcessPoolExecutor` to distribute the simulation tasks across multiple CPU cores.
    *   `_run_single_simulation(atoms: ase.Atoms)`: A private method that constitutes a single simulation job. Its responsibilities are:
        1.  **Late Calculator Binding**: It will instantiate the ASE MLIP calculator (e.g., MACE) *inside* this worker process, not in the main process. This is a crucial optimization to avoid pickling large model objects.
        2.  **Ensemble Switching**: It will analyze the input `atoms` object to detect the presence of a vacuum and accordingly choose between an `NVT` or `NPT` ASE dynamics object.
        3.  **Trajectory Saving**: It will attach an `ase.io.Trajectory` logger to the dynamics object to save frames of the simulation to disk at a regular interval.
        4.  **MD/MC Loop**: It will run the main simulation loop for `n_steps`. If MC moves are configured, it will periodically interrupt the MD run to attempt an atom swap or other MC move.
    *   The engine will be designed to be robust, catching common simulation errors and saving partial results.

*   **`sampling/`**: This new package will house the logic for structure selection.
    *   **`base.py`**: Will define the `BaseSampler` abstract base class with a single method `sample(trajectories: list[Path]) -> list[ase.Atoms]`.
    *   **`random_sampler.py`**: A simple implementation that reads all trajectory files, concatenates all the frames, and then randomly selects the desired number of structures.
    *   **`fps.py`**: The `FarthestPointSampler` implementation. It will use a library like `dscribe` to compute SOAP descriptors for each structure. It will then implement the iterative FPS algorithm to select a subset of structures that are maximally distant from each other in the SOAP feature space.

*   **`pipeline.py` / `services/workflow.py`**: These modules will be updated to integrate the new stages.
    *   The `WorkflowOrchestrator`'s `execute` method will be expanded to a four-stage process:
        1.  Call the `Generator` to get seed structures.
        2.  Instantiate and run the `ExplorationEngine` with the seed structures. This will produce a set of trajectory files.
        3.  Instantiate the appropriate `Sampler` (based on the config) and call its `sample` method with the trajectory file paths. This returns the final, selected list of `ase.Atoms`.
        4.  Call the `AseDBWrapper` to write the final, sampled structures to the database.

## 3. Design Architecture

The design of Cycle 02 continues the schema-first, modular approach. The new components are designed as data-in, data-out processors, which makes them highly testable and reusable.

**Pydantic Schema Design (`config.py` additions):**

*   **`ExplorationConfig(BaseModel)`**:
    *   `model_path: str`: Path to the trained MLIP model file.
    *   `temperature: float = Field(gt=0)`: Simulation temperature in Kelvin.
    *   `pressure: Optional[float] = Field(default=None, gt=0)`: Simulation pressure in bars. If `None`, an NVT ensemble is implied.
    *   `n_steps: int = Field(gt=100)`: Total number of MD steps.
    *   `save_interval: int = Field(gt=10)`: Frequency at which to save frames to the trajectory.
    *   `enable_mc_swaps: bool = False`: Flag to enable/disable atom swapping.
    *   **Producers/Consumers**: Produced by the `FullConfig` parser, consumed by the `ExplorationEngine`.
    *   **Validation**: A `@model_validator` will ensure that `pressure` is only set for 3D (non-slab) systems, enforcing physical correctness.

*   **`SamplingConfig(BaseModel)`**:
    *   `method: Literal["random", "fps"]`: The sampling algorithm to use.
    *   `n_select: int = Field(gt=0)`: The final number of structures to select for the database.
    *   `fps_soap_config: Optional[dict] = None`: A dictionary of parameters for the SOAP descriptor calculation, only relevant if `method` is "fps".
    *   **Validation**: A `@model_validator` will check that if `method` is "fps", then `fps_soap_config` is not `None`.

**Decoupling and Dependency Inversion:**

The architecture will strictly enforce the separation of concerns.
*   **Engines are Pure**: The `ExplorationEngine` and `Sampler` classes are pure computational components. They do not interact with the database directly. Their responsibility is limited to processing input data (structures/trajectories) and producing output data (trajectories/structures).
*   **Orchestrator Handles State**: The `WorkflowOrchestrator` is the *only* component responsible for managing state and orchestrating I/O. It fetches data for the engines and persists their results to the database. This pattern makes the core logic stateless and easier to reason about and test. For example, `ExplorationEngine` doesn't know or care that its output trajectories will be used by a sampler; it just writes them to a specified directory.

**Late Binding of Heavy Objects:**

To ensure performance and scalability, the design explicitly incorporates a "late binding" pattern for the MLIP calculator.
*   The main process, managed by `WorkflowOrchestrator`, will **not** load the MLIP model (which can be hundreds of megabytes).
*   Instead, the `model_path` (a simple string) is passed to the `ExplorationEngine`.
*   The `ExplorationEngine` passes this path to the `ProcessPoolExecutor` when it submits a new simulation task.
*   The worker process (`_run_single_simulation`) receives the string path and only then loads the model and instantiates the ASE calculator. This avoids the massive overhead and potential errors of pickling and sending the large PyTorch model object between processes.

## 4. Implementation Approach

The implementation will focus on building the new components and then integrating them into the existing workflow.

1.  **Update Configuration (`config.py`)**:
    *   Implement the `ExplorationConfig` and `SamplingConfig` Pydantic models.
    *   Add them to the `FullConfig` model.
    *   Implement the new `@model_validator` logic for cross-field validation.

2.  **Implement Samplers (`sampling/`)**:
    *   Start with the simpler `RandomSampler`. Create the `BaseSampler` ABC first. Then, implement the `RandomSampler` which will involve reading ASE trajectory files, collecting all atoms objects, and using Python's `random.sample` to select the final subset.
    *   Next, implement the `FarthestPointSampler`. This will be more involved:
        *   Add `dscribe` as a project dependency.
        *   The `sample` method will first load all structures.
        *   It will then instantiate `dscribe.descriptors.SOAP` using the parameters from the config.
        *   It will call the descriptor's `create()` method to get the feature vectors for all structures.
        *   Finally, it will implement the FPS algorithm: start with a random point, then iteratively find the point farthest from the currently selected set and add it, until `n_select` points are chosen.

3.  **Implement Exploration Engine (`exploration/engine.py`)**:
    *   Create the `ExplorationEngine` class structure.
    *   Implement the `_run_single_simulation` worker function first, as it is the core of the logic.
        *   Use a placeholder calculator like `ase.calculators.emt.EMT` initially for rapid testing.
        *   Write the logic to detect a vacuum (e.g., by checking if the cell's z-dimension is much larger than the atomic coordinates' z-range).
        *   Based on the vacuum detection, choose between `ase.md.langevin.Langevin` (for NVT) or `ase.md.npt.NPT` (for NPT).
        *   Attach `ase.io.Trajectory` and `ase.md.MDLogger` to the dynamics object.
        *   Run the dynamics with `dyn.run(n_steps)`.
    *   Implement the main `run` method using `ProcessPoolExecutor`, mapping the `_run_single_simulation` function over the list of seed structures.
    *   Finally, replace the placeholder EMT calculator with the logic to load a general MLIP calculator (like MACE) from the `model_path` provided in the config.

4.  **Integrate into Workflow (`pipeline.py`, `services/workflow.py`)**:
    *   Modify `WorkflowOrchestrator.execute` to implement the full four-stage logic.
    *   Create a temporary directory for the trajectories.
    *   Pass the list of seed structures to the `ExplorationEngine`.
    *   Pass the path to the trajectory directory to the `Sampler`.
    *   Pass the final list of sampled atoms to the `AseDBWrapper`.

5.  **Update CLI (`cli.py`)**:
    *   Add new optional command-line arguments to override specific exploration or sampling parameters (e.g., `--temperature`, `--n-steps`). This allows for quick experiments without editing the YAML file. Use `typer.Option` and ensure they have sensible defaults (`None`), so they only override the config file if provided.

## 5. Test Strategy

Testing in Cycle 02 must be carefully designed to handle the complexity and long-running nature of simulations.

**Unit Testing Approach (Min 300 words):**

*   **`test_exploration.py`**: Unit testing the `ExplorationEngine` will rely heavily on mocking. We will not run actual simulations.
    *   We will test the **late binding** of the calculator by patching the ASE calculator loading function (e.g., `mace.calculators.mace_mp`) and asserting that it is called with the correct `model_path` from the config.
    *   We will test the **ensemble switching logic** in isolation. We will create mock `ase.Atoms` objects, one representing a bulk crystal and one representing a slab with a large vacuum layer. We will pass these to a helper function within the engine and assert that it correctly returns "NPT" for the bulk system and "NVT" for the slab.
    *   We will use `unittest.mock.patch` to mock the entire ASE `dynamics.run()` method. This allows us to test the setup logic of the `_run_single_simulation` function (e.g., that `Trajectory` and `MDLogger` are attached correctly) without actually executing a long simulation.

*   **`test_samplers.py`**:
    *   The `RandomSampler` will be tested by creating a dummy trajectory file with a known number of frames (e.g., 100) and asserting that the sampler returns a list of atoms with the correct length (`n_select`) and that the selected atoms are indeed a subset of the original ones.
    *   The `FarthestPointSampler` will be tested with a small, 2D dataset where the "farthest" points are visually obvious. We will mock the SOAP descriptor calculation to return these pre-determined 2D points and assert that the FPS implementation correctly picks the expected sequence of points. This validates the algorithm's logic without depending on a complex descriptor library.

**Integration Testing Approach (Min 300 words):**

The integration test will validate the full, four-stage pipeline, using a fast-running classical potential to make the test feasible in a CI environment.

*   **`test_full_pipeline.py`**:
    *   **Setup**: The test will create a `config.yaml` that specifies the full workflow. Crucially, for the `exploration` section, instead of pointing to a real MLIP model, it will specify the use of a built-in, fast ASE calculator like `EMT`. The simulation will be configured to run for a very small number of steps (e.g., 200). The sampling method will be set to `random` with `n_select: 5`.
    *   **Execution**: It will use `typer.testing.CliRunner` to execute the main CLI command with this configuration.
    *   **Validation**: The test will assert the following sequence of events:
        1.  The command exits with a success code.
        2.  A temporary directory for trajectories is created and contains trajectory files.
        3.  The output database file is created.
        4.  The final database contains exactly 5 structures (as specified by `n_select`).
        5.  The energy and forces of the structures in the database are reasonable (i.e., not `NaN` or infinity), confirming that the simulation ran successfully.
    This test provides strong evidence that all components—from configuration parsing to generation, exploration, sampling, and storage—are correctly integrated and that the data flows between them as expected.
