# Specification: CYCLE 02 - Advanced Exploration, Sampling, and UI

## 1. Summary

This document provides the detailed technical specifications for the second and final development cycle of the MLIP-AutoPipe project. This cycle builds directly upon the robust and stable foundation established in Cycle 1. Its primary objective is to implement the scientifically advanced and computationally intensive core of the application, transforming the foundational pipeline into a powerful, intelligent, and user-friendly tool capable of exploring a vast chemical phase space to find diverse and valuable training data for Machine Learning Interatomic Potentials (MLIPs). This phase represents the realization of the project's core scientific value proposition.

The key deliverables for this cycle are threefold and represent the most complex components of the system: the `ExplorationEngine`, the `SamplingEngine`, and a user-friendly Web UI. The `ExplorationEngine` will be the computational heart of the system, responsible for running massively parallelized molecular dynamics (MD) and hybrid Monte Carlo (MC) simulations. This engine will not just be a simple wrapper around a simulation library; it will incorporate sophisticated, domain-specific logic, such as the automatic switching of thermodynamic ensembles (NVT/NPT) based on the physical system being simulated and the integration of a ZBL potential to handle high-energy atomic interactions gracefully and prevent simulation crashes.

Following the exploration stage, the `SamplingEngine` will be responsible for intelligently processing the vast amount of trajectory data generated. A key feature of this engine will be a Farthest Point Sampling (FPS) implementation, which will select a subset of the most structurally diverse configurations from the tens of thousands of generated frames. This is a critical step to ensure the final dataset is rich in information and not redundant, maximizing the efficiency of the MLIP training process. Finally, to make this powerful backend accessible to a wider audience, including those who are not command-line experts, a Web UI will be developed. This graphical interface will provide an intuitive way for users to configure, launch, and monitor the data generation pipeline, as well as to visualize the resulting atomic structures. By the end of this cycle, MLIP-AutoPipe will be a feature-complete, scientifically robust, and highly usable tool, ready to be deployed to the materials science community.

## 2. System Architecture

The architecture in Cycle 2 undergoes a significant expansion. It introduces the new, computationally intensive modules for exploration and sampling and integrates them seamlessly into the main pipeline. A new user-facing entry point, the Web UI, is also added, which will interact with the existing backend pipeline.

**File Structure (Cycle 2 Focus):**
The files to be created or modified in this cycle are marked in bold. This represents the primary development effort for this cycle.

```
.
├── src
│   └── mlip_autopipec
│       ├── **main_gui.py**           # **[CREATE]** Web UI entry point using Streamlit.
│       │
│       ├── pipeline
│       │   └── **runner.py**         # **[MODIFY]** Integrate new Exploration and Sampling stages.
│       │
│       ├── exploration
│       │   ├── **__init__.py**
│       │   └── **engine.py**         # **[CREATE]** ExplorationEngine, with parallel MD/MC logic.
│       │
│       ├── sampling
│       │   ├── **__init__.py**
│       │   ├── **base.py**           # **[CREATE]** Abstract BaseSampler class.
│       │   ├── **random.py**         # **[CREATE]** Simple RandomSampler implementation.
│       │   ├── **fps.py**            # **[CREATE]** Advanced FPSSampler implementation.
│       │   └── **factory.py**        # **[CREATE]** Factory to select the appropriate sampler.
│       │
│       └── ... (existing modules from C1)
│
└── ...
```

**Component Interaction Blueprint:**

1.  **`pipeline/runner.py` (Modified `PipelineRunner`):**
    *   The `run()` method, which was a simple two-stage process in Cycle 1, will be expanded to orchestrate the full four-stage process: **Generation -> Exploration -> Sampling -> Storage**.
    *   After the initial structures are generated, it will now instantiate and invoke the `ExplorationEngine`.
    *   The raw trajectories (represented by a list of file paths) produced by the `ExplorationEngine` will be passed as input to the `SamplingEngine`.
    *   The final, curated list of sampled structures will then be passed to the `StorageEngine` for archival in the database.
    *   **Code Blueprint (Conceptual):**
        ```python
        # In PipelineRunner.run()
        # Stage 1: Generation (from C1)
        initial_structures = self.generator.generate()
        self.storage_engine.save(initial_structures, "initial")

        # Stage 2: Exploration
        exploration_engine = ExplorationEngine(self.config.exploration)
        trajectory_paths = exploration_engine.run(initial_structures)

        # Stage 3: Sampling
        sampler = SamplerFactory.get_sampler(self.config.sampling)
        sampled_structures = sampler.sample(trajectory_paths)

        # Stage 4: Storage
        self.storage_engine.save(sampled_structures, "final_dataset")
        ```

2.  **`exploration/engine.py` (`ExplorationEngine`):**
    *   This component is the computational workhorse of the application. It will be designed to manage the parallel execution of many independent MD/MC simulations.
    *   It will leverage Python's `concurrent.futures.ProcessPoolExecutor` to run these simulations concurrently, allowing it to take full advantage of multi-core processors.
    *   It will implement the critical "late-binding" calculator pattern. This means that the large MLIP calculator object (which can be hundreds of megabytes) will be instantiated *inside* each worker process, rather than being created in the main process and pickled to the workers. This is a crucial optimization that avoids significant inter-process communication overhead and potential serialization errors.
    *   The core logic for the hybrid MD/MC simulations will reside here, including the logic for periodically attempting Monte Carlo atom swap moves during a molecular dynamics run.
    *   **Code Blueprint (Conceptual):**
        ```python
        from concurrent.futures import ProcessPoolExecutor
        from ase.atoms import Atoms
        from pathlib import Path
        # ... other imports

        def run_single_md_process(structure: Atoms, config: ExplorationConfig, output_dir: Path) -> Path:
            """Worker function to run one MD simulation in a separate process."""
            # 1. Instantiate the MLIP calculator HERE (late-binding pattern)
            #    calculator = MACE(model_path=config.mlip_model_path)
            # 2. Setup the ASE dynamics object (e.g., Langevin dynamics) with the calculator.
            # 3. Implement the main simulation loop.
            # 4. Inside the loop, periodically attempt an MC swap move.
            # 5. Attach a trajectory logger to save frames to a unique file in output_dir.
            # 6. Return the path to the completed trajectory file.
            pass

        class ExplorationEngine:
            def __init__(self, config: ExplorationConfig):
                self.config = config

            def run(self, structures: list[Atoms]) -> list[Path]:
                with ProcessPoolExecutor() as executor:
                    # Create a temporary directory for trajectories
                    with tempfile.TemporaryDirectory() as tmpdir:
                        output_dir = Path(tmpdir)
                        # Farm out the work to the process pool
                        futures = [executor.submit(run_single_md_process, s, self.config, output_dir) for s in structures]
                        # Gather the paths of the completed trajectory files
                        results = [f.result() for f in futures]
                return results
        ```

## 3. Design Architecture

The design in Cycle 2 is focused on managing the significant increase in complexity associated with the scientific simulation parameters and the corresponding workflows. The Pydantic schemas defined in `config.py` will be extended to provide a robust and validated configuration interface for the new `exploration` and `sampling` stages. This ensures that the powerful new features are exposed to the user in a clear, predictable, and safe manner.

**`config.py` - Pydantic Schema Extensions:**
*   **`ExplorationConfig`:** This new Pydantic model will define the complete schema for the `exploration` section of the user's YAML file. It will include fields for all parameters that control the MD/MC simulations.
    *   `md_type`: A string field constrained to an `Enum` containing valid ASE thermodynamic ensembles, such as `"nvt"` (Langevin dynamics) and `"npt"` (NPT dynamics).
    *   `temperature_k`: This field will be a `float` or a `list[float]`, allowing the user to specify either a constant temperature or a temperature ramp for simulated annealing protocols. It will be validated to be greater than zero.
    *   `pressure_gpa`: An optional `float`, which is only applicable when `md_type` is `"npt"`.
    *   `mc_swap_probability`: A `float` constrained to be between 0.0 and 1.0 using `Field(ge=0, le=1)`. This controls the frequency of Monte Carlo swap moves.
    *   `mlip_model_path`: A `Path` that points to the trained MLIP model file. Pydantic will validate that this path exists.
*   **`SamplingConfig`:** This model will define the schema for the `sampling` section.
    *   `sampler_type`: An `Enum` to allow the user to select a valid sampling algorithm, e.g., `"random"` or `"fps"`.
    *   `num_samples`: A `PositiveInt` that specifies the desired number of final structures to be selected from the trajectories.
    *   `fps_soap_params`: A nested Pydantic model that defines all the parameters for the SOAP descriptor calculation required by the FPS algorithm. This includes fields like `n_max`, `l_max`, and `r_cut`. This nested structure helps to organize the configuration and make it more readable.

**Producers and Consumers of Data:**
*   The `ExplorationEngine` is a major **producer** of data, consuming a small set of initial structures and producing a large volume of raw trajectory data, which is stored on the filesystem.
*   The `SamplingEngine` is both a **consumer** and a **producer**. It consumes the raw trajectory file paths and produces the final, curated, and much smaller list of `ase.Atoms` objects.
*   The `StorageEngine` is the final **consumer** in the pipeline, taking the curated list of structures and archiving it.

**Extensibility and Scientific Robustness:**
The design of the `ExplorationEngine` must be scientifically robust. The logic for automatic ensemble switching (i.e., detecting a vacuum slab in the initial structure and forcing the use of the NVT ensemble to prevent simulation artifacts) is a critical piece of embedded domain knowledge that will be implemented as a preparatory step before the main simulation run. Similarly, the `BaseSampler` abstract class in the sampling module is a key design pattern for extensibility. It will define a simple `sample(trajectories)` interface, which will ensure that adding new and innovative sampling algorithms in the future is a straightforward process that does not require any changes to the core pipeline runner.

## 4. Implementation Approach

The implementation will be staged to tackle the complex backend logic first, ensuring it is well-tested and robust before building the user-facing Web UI on top of it.

1.  **Extend Configuration (`config.py`):** The first step is to implement the new `ExplorationConfig` and `SamplingConfig` Pydantic models. This schema-first approach ensures that all the required parameters for the new stages are well-defined, validated, and documented before the implementation of the logic that uses them begins.
2.  **Implement Samplers (`sampling/`):** Implement the `BaseSampler` abstract class, which will define the `sample()` method. Then, implement the `RandomSampler` as a simple baseline. Following that, implement the more complex `FPSSampler`. This will involve integrating a third-party library, such as `dscribe`, to perform the heavy lifting of computing the SOAP descriptors that are required as input for the FPS algorithm. Finally, create the `SamplerFactory` to select and instantiate the appropriate sampler based on the configuration.
3.  **Implement Exploration Engine (`exploration/engine.py`):** This is the most complex and critical implementation task for this cycle.
    *   First, develop the `run_single_md_process` worker function. This function will be the core of the simulation logic. It will encapsulate everything needed for a single simulation run: setting up the ASE `Atoms` object with the correct calculator, defining the ASE dynamics object (e.g., `Langevin`), attaching a trajectory logger, and running the simulation loop.
    *   Next, implement the `ExplorationEngine` class that uses `ProcessPoolExecutor` to call this worker function in parallel for each of the initial structures.
    *   Careful management of file I/O is critical. The engine must ensure that each worker process writes to a unique temporary trajectory file within a managed temporary directory to prevent race conditions.
4.  **Integrate into Pipeline (`pipeline/runner.py`):** Modify the `PipelineRunner` to incorporate the new `ExplorationEngine` and `SamplingEngine`. The `run()` method will be updated to execute the full four-stage data flow. The runner will need to correctly manage the intermediate data (i.e., the list of paths to the trajectory files) that is passed between the new stages.
5.  **Implement Web UI (`main_gui.py`):**
    *   The Web UI will be implemented using Streamlit, chosen for its simplicity and rapid development capabilities, which are ideal for data-centric applications.
    *   The UI will be organized into logical sections using tabs or expanders for configuring the `system`, `exploration`, and `sampling` parameters. The interactive widgets (e.g., sliders for temperature, dropdowns for enums) will dynamically update a configuration object in the session state.
    *   A prominent "Run Pipeline" button will trigger the execution of the `PipelineRunner`. To prevent the UI from freezing during the long-running pipeline, the runner will be executed in a background thread or a separate process.
    *   A status area (e.g., a `st.text_area` or `st.progress`) will be used to display real-time log updates and progress from the pipeline, providing essential feedback to the user.
    *   Upon completion, a results section will be populated. It will allow the user to visualize some of the final generated structures directly in the browser using a library like `py3Dmol`, providing immediate and gratifying visual confirmation of the results.
6.  **Testing:** As each new component is developed, a corresponding suite of unit and integration tests will be written.

## 5. Test Strategy

Testing in Cycle 2 is crucial and more complex than in Cycle 1 due to the computationally intensive nature and non-determinism of the simulations. The strategy will rely heavily on mocking external processes and focusing on the correctness of the data flow and the internal logic of the components.

**Unit Testing Approach:**
*   **`test_exploration.py`:** Testing the `ExplorationEngine` requires a specialized approach that avoids running actual, time-consuming MD simulations.
    *   We will heavily mock the `run_single_md_process` function. The mock will be configured to bypass the real simulation and instead simply create a dummy trajectory file with a few known structures and return its path.
    *   The test for `ExplorationEngine.run()` will then call the method with a list of `Atoms` objects and assert that the mocked worker function was called the correct number of times (once per input structure). This effectively verifies the parallelization logic without the computational cost.
    *   The worker function itself (`run_single_md_process`) will be tested in isolation. In this test, the actual call to the ASE dynamics object (`dyn.run()`) will be mocked. The test will focus on asserting that the ASE dynamics object was initialized with the correct parameters (temperature, timestep, friction, etc.) that were derived from the Pydantic configuration object. This verifies that the configuration is being correctly translated into the simulation setup.
*   **`test_sampling.py`:**
    *   The `RandomSampler` is simple to test: we will give it a trajectory file and assert that it returns the correct number of randomly selected structures.
    *   Testing the `FPSSampler` is more involved but critical. We will create a test fixture with a small, known set of `Atoms` objects where some structures are very similar to each other, and one is very different. We will mock the expensive SOAP descriptor calculation to return pre-computed, deterministic feature vectors for these structures. The test will then call `FPSSampler.sample()` and assert that the algorithm correctly selects the most diverse structures (i.e., the ones whose feature vectors are farthest apart in the metric space), based on the known ground truth.

**Integration Testing Approach:**
*   **`test_full_pipeline.py`:** A new integration test will be created to cover the entire, end-to-end, four-stage pipeline.
    *   This test will be carefully configured to be as fast and lightweight as possible. It will run on a very small and simple system (e.g., a 2-atom H2 molecule) and for a very short simulation time (e.g., only 10 MD steps). Crucially, it will use the built-in, fast EMT potential from ASE instead of a real, slow MLIP.
    *   The test will use the `Typer.testing.CliRunner` to invoke the pipeline with a configuration that enables all four stages, including exploration and FPS sampling.
    *   The final assertions will check the integrity of the database. It will verify that a `final_dataset` group has been created and that it contains the expected number of structures after the full workflow has been completed. This test, while slower than a unit test, is absolutely essential for confirming that the data contracts between the different stages are correctly implemented and that the entire system works in concert. It will be marked as a slow test in `pytest` to be run selectively.
