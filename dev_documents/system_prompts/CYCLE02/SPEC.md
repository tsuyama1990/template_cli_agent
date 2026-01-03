# CYCLE02/SPEC.md

## 1. Summary

This specification outlines the technical details for Cycle 2 of the MLIP-AutoPipe project. Building upon the foundational structure generation capabilities developed in Cycle 1, this cycle introduces the dynamic and intelligent core of the application. The primary goals are to implement the **Exploration** and **Sampling** stages of the pipeline, and to provide a user-friendly **Web UI** for interactive workflow management. This cycle will transform the tool from a static structure generator into a powerful, automated framework for discovering diverse and physically relevant atomic configurations. The Exploration stage will feature a sophisticated hybrid Molecular Dynamics (MD) and Monte Carlo (MC) engine, designed to efficiently sample the potential energy surface of a material. The Sampling stage will implement intelligent algorithms to select the most informative structures from the vast amount of data generated during exploration.

Key technical features to be developed include a robust MD engine capable of parallel execution, the integration of MLIP calculators with a classical ZBL potential for short-range stability, and an automatic ensemble switching mechanism (NVT/NPT) based on the physical geometry of the system. This ensures that simulations are both stable and physically appropriate. For the Sampling stage, we will implement both a simple random sampler and a more advanced Farthest Point Sampler (FPS) using SOAP (Smooth Overlap of Atomic Positions) descriptors to maximize the structural diversity of the final dataset. The `PipelineRunner` developed in Cycle 1 will be extended to incorporate these new stages, creating the full, four-stage data generation workflow. Finally, a Web UI will be developed using a framework like Gradio or Streamlit, providing a graphical interface for users to configure, run, and monitor the pipeline, as well as visualize the results. This will significantly lower the barrier to entry and make the tool accessible to a wider audience.

## 2. System Architecture

Cycle 2 integrates the Exploration and Sampling components into the existing architecture and adds a new entry point for the Web UI.

**File Structure (Cycle 2 Focus):**

Files to be created or modified in this cycle are marked in bold.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py
├── **web.py**                  # Web User Interface entry point using Gradio/Streamlit
├── config/
│   ├── __init__.py
│   └── **schemas.py**          # Add Pydantic models for Exploration & Sampling
├── pipeline/
│   ├── __init__.py
│   └── **runner.py**           # Extend PipelineRunner to a 4-stage workflow
├── generators/
│   └── ...                 # (No changes in this cycle)
├── explorers/
│   ├── __init__.py
│   └── **md_engine.py**        # Hybrid MD/MC exploration engine
├── samplers/
│   ├── __init__.py
│   ├── **base.py**             # Abstract BaseSampler class
│   ├── **factory.py**          # SamplerFactory implementation
│   ├── **random.py**           # RandomSampler implementation
│   └── **fps.py**              # FPSSampler implementation
├── storage/
│   └── ...                 # (No changes in this cycle)
└── utils/
    ├── __init__.py
    └── **physics.py**          # Add vacuum detection logic
```

**Code Blueprints:**

-   **`pipeline/runner.py` (Extended)**:
    ```python
    from mlip_autopipec.config.schemas import FullConfig
    from mlip_autopipec.generators.factory import GeneratorFactory
    from mlip_autopipec.explorers.md_engine import MDEngine
    from mlip_autopipec.samplers.factory import SamplerFactory
    from mlip_autopipec.storage.db_wrapper import AseDBWrapper

    class PipelineRunner:
        def __init__(self, config: FullConfig):
            self.config = config

        def run(self):
            # Stage 1: Generation (from Cycle 1)
            generator = GeneratorFactory.create_generator(self.config.system)
            seed_structures = generator.generate()

            # Stage 2: Exploration
            explorer = MDEngine(self.config.exploration)
            trajectories = explorer.run_simulations(seed_structures)

            # Stage 3: Sampling
            sampler = SamplerFactory.create_sampler(self.config.sampling)
            final_structures = sampler.sample(trajectories)

            # Stage 4: Storage (from Cycle 1)
            db_wrapper = AseDBWrapper(filepath=self.config.storage.db_path)
            db_wrapper.write_structures(final_structures)
    ```

-   **`explorers/md_engine.py`**:
    ```python
    from typing import List
    from ase import Atoms
    from ase.md.langevin import Langevin
    # ... other ASE imports
    from mlip_autopipec.config.schemas import ExplorationConfig
    from mlip_autopipec.utils.physics import detect_vacuum
    import concurrent.futures

    class MDEngine:
        def __init__(self, config: ExplorationConfig):
            self.config = config

        def _get_calculator(self):
            # Logic to instantiate MLIP + ZBL mixed potential
            # This is a "late binding" to avoid pickling issues
            pass

        def _run_single_md(self, atoms: Atoms) -> List[Atoms]:
            # 1. Attach the calculator
            # 2. Use detect_vacuum to choose NVT or NPT ensemble
            # 3. Set up the MD simulation (e.g., Langevin dynamics)
            # 4. Run the simulation for configured number of steps
            # 5. Periodically save frames to a trajectory list
            # 6. Return the trajectory
            trajectory = []
            # ... MD logic ...
            return trajectory

        def run_simulations(self, structures: List[Atoms]) -> List[List[Atoms]]:
            all_trajectories = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.config.max_workers) as executor:
                future_to_atoms = {executor.submit(self._run_single_md, s): s for s in structures}
                for future in concurrent.futures.as_completed(future_to_atoms):
                    all_trajectories.append(future.result())
            return all_trajectories
    ```

-   **`samplers/fps.py`**:
    ```python
    from typing import List
    from ase import Atoms
    from dscribe.descriptors import SOAP
    # FPS implementation algorithm from a library or custom

    from .base import BaseSampler

    class FPSSampler(BaseSampler):
        def _calculate_soap_descriptors(self, trajectory: List[Atoms]):
            # Use dscribe library to compute SOAP vectors for each frame
            pass

        def sample(self, trajectories: List[List[Atoms]]) -> List[Atoms]:
            # 1. Flatten all trajectories into a single list of frames
            # 2. Calculate SOAP descriptors for all frames
            # 3. Apply the FPS algorithm to select a diverse subset of frames
            # 4. Return the selected Atoms objects
            pass
    ```

## 3. Design Architecture

This cycle introduces new Pydantic schemas for the Exploration and Sampling stages, which will be added to the main `FullConfig` model.

**Configuration Schemas (`config/schemas.py` additions):**

-   **`MLIPConfig`**: Defines the machine learning potential.
    -   `model_path` (str): Path to the trained MLIP model file (e.g., a MACE .model file).
    -   `calculator_type` (Literal["mace", "sevennet"]): The type of MLIP to use.

-   **`MCConfig`**: Defines Monte Carlo move settings.
    -   `enable_swap_moves` (bool): Whether to perform atom swaps.
    -   `swap_frequency` (int): Perform a swap attempt every N steps. Constraint: `gt=0`.

-   **`ExplorationConfig(BaseModel)`**:
    -   `mlip: MLIPConfig`
    -   `mc: Optional[MCConfig]`
    -   `temperature_k` (float): Simulation temperature in Kelvin. Constraint: `gt=0`.
    -   `pressure_gpa` (Optional[float]): Simulation pressure in GPa. If None, NVT is used unless vacuum is detected.
    -   `num_steps` (int): Number of MD steps to run for each structure. Constraint: `gt=0`.
    -   `max_workers` (int): Number of parallel processes to use. Constraint: `ge=1`.

-   **`SamplingConfig(BaseModel)`**:
    -   `method` (Literal["random", "fps"]): The sampling algorithm to use.
    -   `num_samples` (int): The target number of structures to select for the final dataset.

-   **`FullConfig(BaseModel)` (Extended)**:
    -   `system: ...`
    -   `validation: ...`
    -   `exploration: ExplorationConfig`  # New
    -   `sampling: SamplingConfig`      # New
    -   `storage: ...`

**Data Flow and Interactions:**

The data flow is extended from Cycle 1. The list of seed `Atoms` objects from the `Generator` is now the input for the `MDEngine`. The `MDEngine` produces a list of trajectories (a list of lists of `Atoms` objects). This complex data structure is then consumed by the `Sampler`, which flattens it and selects a final, smaller list of `Atoms` objects. This final list is then passed to the `AseDBWrapper`, as before.

The **Web UI (`web.py`)** will be a new producer of configuration. It will use the `FullConfig` Pydantic model to dynamically generate a web form. When the user submits the form, the web application will construct the `FullConfig` object and pass it to the `PipelineRunner` to execute the workflow, likely in a background thread or process to avoid blocking the UI.

## 4. Implementation Approach

The implementation will be phased to tackle the most complex parts first.

1.  **Extend Configuration**: First, update `config/schemas.py` with the new Pydantic models: `MLIPConfig`, `MCConfig`, `ExplorationConfig`, and `SamplingConfig`. Integrate them into the main `FullConfig` model.
2.  **MD Engine Implementation**: This is the most complex part of the cycle.
    -   Implement the `MDEngine` class in `explorers/md_engine.py`.
    -   Begin with the parallel execution logic using `concurrent.futures.ProcessPoolExecutor`.
    -   Implement the `_get_calculator` method. This will involve using the ASE interface for MACE (or other MLIPs) and combining it with the `ase.calculators.mixing.MixedCalculator` to add the ZBL potential.
    -   Implement the `_run_single_md` method. Add the logic to call `detect_vacuum` (from `utils/physics.py`, to be implemented next) and choose between `ase.md.velocityscaling.VelocityVerlet` (for NVT) and `ase.md.npt.NPT` (for NPT) dynamics.
3.  **Vacuum Detection**: Implement the `detect_vacuum` function in `utils/physics.py`. This will involve creating a 3D grid over the simulation cell, marking voxels that are close to an atom, and then checking for contiguous regions of "empty" voxels, particularly along the cell axes, which would indicate a vacuum slab.
4.  **Sampler Implementation**:
    -   Create the `BaseSampler` abstract class in `samplers/base.py`.
    -   Implement the `SamplerFactory` in `samplers/factory.py`.
    -   Implement the `RandomSampler` as a simple baseline. It will flatten all trajectories and randomly select `num_samples` frames.
    -   Implement the `FPSSampler`. This will require adding `dscribe` as a dependency. The implementation will first calculate SOAP descriptors for all frames in the trajectories and then apply an FPS algorithm to select the most diverse subset.
5.  **Pipeline Integration**: Modify `pipeline/runner.py`. Extend the `run` method to instantiate and call the `MDEngine` and the chosen `Sampler` between the Generation and Storage stages.
6.  **Web UI Development**: Once the backend pipeline is complete and tested, develop the Web UI in `web.py`.
    -   Choose a framework (e.g., Gradio).
    -   Use the `FullConfig` Pydantic schema to generate the input form fields.
    -   Write a function that takes the user inputs from the form, constructs the `FullConfig` object, and calls `PipelineRunner.run()`.
    -   Provide a text area for logging output and some form of progress indication. A simple visualization of the final structures using `matplotlib` or a similar library can be added as a final touch.

## 5. Test Strategy

Testing in Cycle 2 must account for the complexity and stochasticity of MD simulations.

**Unit Testing Approach (Min 300 words):**

Unit tests will heavily rely on mocking the expensive computational parts.
-   **MD Engine**: Testing the `MDEngine` will not involve running actual MD. Instead, we will mock the ASE dynamics objects. For the automatic ensemble switching, we will create two hand-crafted `Atoms` objects: one representing a bulk crystal and one with a large vacuum gap. We will mock `utils.physics.detect_vacuum` to return `False` for the bulk and `True` for the slab. We will then assert that `_run_single_md` attempts to instantiate the correct ASE ensemble class (`NPT` for bulk, `NVT` for slab). We will also test the late binding of the calculator by asserting that `_get_calculator` is called within the subprocess, not in the main thread.
-   **Samplers**: The `RandomSampler` will be tested by providing a known list of 100 `Atoms` objects and requesting 10 samples. We will assert that it returns exactly 10 objects and that all returned objects are members of the original list. For the `FPSSampler`, we will create a small, deterministic set of `Atoms` objects with pre-calculated SOAP descriptors. We will run the FPS algorithm and assert that it selects the expected subset of structures based on their known diversity. This will validate the correctness of the sampling logic without relying on expensive descriptor calculations in the test suite.

**Integration Testing Approach (Min 300 words):**

Integration tests will cover the full four-stage pipeline, using a fast classical potential to make the tests feasible in a CI environment.
-   **Full Pipeline Run**: We will create a test configuration that uses a fast calculator like `ase.calculators.emt.EMT`. The test will configure a very short MD run (e.g., 2 structures, 10 steps each). The test will use `CliRunner` to execute the full pipeline. After completion, it will inspect the final database. The primary assertion will be that the database contains the number of structures requested in the `sampling` configuration (e.g., `num_samples: 5`). We will also read the structures and check their energy, asserting that it is a plausible value as calculated by the EMT potential. This end-to-end test validates that all components are correctly wired together.
-   **Web UI Smoke Test**: We will write a simple "smoke test" for the Web UI. This test will not perform complex UI interactions but will ensure the application can be launched without crashing. Using a library like `requests` or `playwright`, the test will start the web server in a separate process, make an HTTP request to the main page, and assert that it receives a 200 OK status code. This provides a basic check that the UI and its dependencies are correctly installed and configured.
