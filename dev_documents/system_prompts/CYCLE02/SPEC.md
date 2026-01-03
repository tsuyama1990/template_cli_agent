# CYCLE 02: SPECIFICATION

## 1. Summary

This document provides the detailed technical specifications for the second and final implementation cycle of the MLIP-AutoPipe project. This cycle is designed to build upon the robust architectural foundation laid in Cycle 01, introducing the advanced, intelligent, and performance-critical features that will elevate the system from a simple automation tool to a fully autonomous, self-improving learning platform. The central and overriding objective of this cycle is to implement and validate the active learning loop. This transforms the linear, "one-shot" pipeline of Cycle 01 into a dynamic, cyclical process where the system can intelligently identify and correct its own deficiencies. This is the core value proposition of the entire project.

The implementation in this cycle will focus on bringing the remaining, previously stubbed-out modules to life: the `Explorer & Sampler` (Module B) and the `SimulationEngine` (Module E). The integration of a powerful, pre-trained universal surrogate model like MACE into the `Explorer & Sampler` is a key feature, as it will enable the system to perform vast explorations of a material's potential energy surface at a trivial computational cost compared to DFT. This will be paired with the DIRECT sampling algorithm to ensure that the initial datasets are of the highest possible quality and diversity. The centerpiece of the cycle, however, is the "on-the-fly" (OTF) `SimulationEngine`. This module will deploy the MLIPs trained by the system into live molecular dynamics simulations, constantly monitor their predictive uncertainty, and trigger the crucial feedback loop by sending novel, high-uncertainty atomic configurations back for DFT labeling and subsequent model retraining.

Performance is a non-negotiable requirement for this cycle. The new modules involve computationally intensive tasks, particularly the calculation of atomic environment descriptors for the DIRECT sampler, which must be performed for millions of simulation frames. To address this, these performance-critical sections of the code will be aggressively optimized using Numba, a high-performance JIT compiler, to achieve speeds comparable to compiled languages like C or Fortran. The `WorkflowOrchestrator` will also undergo a significant refactoring to manage the new, stateful, and cyclical data flow of the active learning process. By the end of Cycle 02, the MLIP-AutoPipe system will be feature-complete, fulfilling the initial vision of a platform that can autonomously generate a high-fidelity MLIP from a minimal user input by intelligently and efficiently navigating the vast space of atomic configurations.

## 2. System Architecture

The system architecture in Cycle 02 evolves from the linear pipeline of Cycle 01 into a sophisticated, closed-loop system. This evolution is primarily reflected in the promotion of the `Explorer & Sampler` and `SimulationEngine` from stubs to fully implemented components, and in the significant enhancement of the `WorkflowOrchestrator`'s logic to manage the new, cyclical workflow. The file structure remains consistent with the blueprint laid out in the previous cycle, but the newly implemented files are now marked in **bold**.

**File Structure Update:**

```
mlip_autopipec/
├── pyproject.toml
├── uv.lock
└── src/
    └── mlip_autopipec/
        ├── __init__.py
        ├── cli.py
        ├── config/
        │   ├── __init__.py
        │   ├── models.py      # **Updated with ActiveLearningConfig sub-model**
        │   └── loader.py
        ├── data/
        │   ├── __init__.py
        │   └── models.py
        ├── database/
        │   ├── __init__.py
        │   └── ase_db_wrapper.py
        ├── workflows/
        │   ├── __init__.py
        │   └── orchestrator.py  # **Refactored to manage the active learning loop**
        ├── interfaces/
        │   ├── __init__.py
        │   └── engines.py     # **Updated with IExplorerSampler and ISimulationEngine ABCs**
        └── engines/
            ├── __init__.py
            ├── structure_generator.py
            ├── **explorer_sampler.py**    # **Module B: Full implementation with MACE and DIRECT**
            ├── labeling_engine.py
            ├── training_engine.py
            └── **simulation_engine.py**   # **Module E: Full implementation of OTF simulation**
```

**Architectural Flow for Cycle 02:**

The workflow is now a stateful, multi-generation loop, representing the core intelligence of the system.

1.  **Initial Seeding and "Generation 0" Training**: The process can now start in one of two ways, determined by the `FullConfig`. It can use the basic `StructureGenerator` from Cycle 01, or, more powerfully, it can use the newly implemented `ExplorerSampler` (Module B). When using the latter, the system will first run a large-scale, high-temperature molecular dynamics simulation using the fast MACE surrogate model. The resulting multi-million-frame trajectory is then passed to the DIRECT sampler, which uses Numba-accelerated descriptor calculations and clustering to select a maximally diverse and informative set of initial structures. This initial set of structures is then passed to the `LabelingEngine` (Module C) for DFT calculation, and the results are used by the `TrainingEngine` (Module D) to create the initial, "Generation 0" MLIP.

2.  **The Active Learning Loop**: This is the main, iterative phase of the workflow, managed by the `WorkflowOrchestrator`.
    a. **Deployment**: The orchestrator takes the latest trained MLIP (e.g., the "Generation 0" model) and passes its file path to the `SimulationEngine` (Module E).
    b. **On-the-Fly Simulation**: The `SimulationEngine` loads the MLIP and begins a new molecular dynamics simulation. At regular intervals (e.g., every 10 steps), it calculates an uncertainty metric for the current atomic configuration.
    c. **Uncertainty Trigger**: The calculated uncertainty is compared against a dynamic threshold. This threshold is not a fixed value but is intelligently calculated based on the distribution of uncertainties on the model's *current* training set (e.g., the 95th percentile).
    d. **Yielding**: If the uncertainty exceeds this threshold, the `SimulationEngine` pauses the simulation and `yields` the high-uncertainty `AtomsModel` back to the orchestrator. This generator-based approach allows the orchestrator to process the new structure without terminating the simulation process.
    e. **Labeling**: The orchestrator receives the yielded structure and immediately sends it to the `LabelingEngine` to get its ground-truth DFT energy and forces.
    f. **Augmentation**: The orchestrator adds this new, labeled data point to the master dataset within the database.
    g. **Retraining**: Once the simulation for a given generation concludes (e.g., after a fixed number of MD steps), the orchestrator gathers the *entire* augmented dataset from the database and calls the `TrainingEngine` again. This produces a new, improved "Generation 1" MLIP.
    h. **Iteration**: The process repeats from step 2a, this time deploying the improved "Generation 1" model. This iterative refinement loop continues until a user-defined stopping criterion is met, such as a maximum number of generations or the rate of finding new, uncertain structures dropping below a certain frequency.

## 3. Design Architecture

The design for Cycle 02 extends the rigorous, schema-first foundation of Cycle 01 to the new, more dynamic components. The Pydantic models are expanded to provide a validated and type-safe configuration for the active learning process, and the abstract interfaces are updated to define the contracts for the new engines. This ensures that the new, complex features are integrated in a way that maintains the overall robustness and modularity of the system.

**`config.models.py` - Expanded Configuration Schema:**
The `FullConfig` model will be updated to include a new, nested Pydantic model that specifically controls the active learning workflow. This encapsulates all the relevant parameters in a single, validated namespace.
-   `ActiveLearningConfig`: This new `BaseModel` will be added as a field within the `FullConfig`. It will contain all the necessary knobs to tune the active learning process, ensuring they are type-safe and validated. Its fields will include:
    -   `strategy: Literal["direct_then_active", "active_only"] = "direct_then_active"`: A validated string literal to control the overall workflow.
    -   `surrogate_model: str = "mace_mp"`: The identifier for the surrogate model to be used by the `ExplorerSampler`.
    -   `max_generations: int = Field(10, gt=0)`: The maximum number of retraining cycles, validated to be a positive integer.
    -   `uncertainty_threshold_strategy: Literal["dynamic_95percentile", "fixed"] = "dynamic_95percentile"`: Defines the method for setting the uncertainty trigger.
    -   `fixed_uncertainty_threshold: Optional[float] = None`: An optional fixed threshold value.
    -   `md_steps_per_generation: int = Field(10000, gt=100)`: The number of MD steps to run in the `SimulationEngine` for each generation, validated to be a reasonable number.

**`interfaces.engines.py` - Expanded Engine Contracts:**
The `interfaces` module will be updated with new Abstract Base Classes to define the roles of the new engines, ensuring they can be seamlessly integrated into the `WorkflowOrchestrator`.
-   `IExplorerSampler(ABC)`: This will define the contract for the initial, surrogate-based sampling. It will have a single abstract method: `@abstractmethod def sample_structures(self, config: FullConfig) -> list[AtomsModel]: pass`.
-   `ISimulationEngine(ABC)`: This will define the contract for the on-the-fly simulation engine. Its design will be asynchronous and generator-based to handle the yielding of structures. The core abstract method will be: `@abstractmethod def run_simulation(self, model_path: Path, config: FullConfig) -> AsyncGenerator[AtomsModel, None]: pass`. This signature clearly defines that the engine takes a trained model and the configuration, and it asynchronously yields `AtomsModel` objects back to the caller.

**Data Flow and Component Responsibilities in Cycle 02:**
-   **ExplorerSampler (Module B)**: This engine is a consumer of the `FullConfig` and a high-volume producer of initial `AtomsModel` objects. Its responsibility is to encapsulate all the logic for interacting with the MACE surrogate model API, running the exploratory MD simulation, and then executing the performance-critical DIRECT sampling algorithm. The descriptor calculation and clustering logic is entirely contained within this module.
-   **SimulationEngine (Module E)**: This engine is a consumer of a trained MLIP model file and the `FullConfig`. It is a producer of a stream of high-uncertainty `AtomsModel` objects. Its responsibilities are extensive: it must correctly interface with the chosen MD backend (e.g., ASE's built-in integrators), load the custom-trained MLIP as the force calculator, implement the chosen uncertainty quantification algorithm (e.g., committee disagreement for an ensemble model), continuously apply the dynamic thresholding logic, and manage the generator-based communication with the orchestrator.
-   **WorkflowOrchestrator**: The orchestrator's role becomes significantly more complex and stateful. It is now responsible for managing the entire multi-generation active learning loop. It must maintain the state of the current generation number, orchestrate the initial "Generation 0" training, and then enter the main loop. Inside the loop, it is the component that receives structures yielded by the `SimulationEngine` and immediately routes them to the `LabelingEngine`, waits for the result, and persists it to the database. It is also responsible for deciding when to end a generation's simulation and trigger a full model retraining with the `TrainingEngine`. Finally, it enforces the `max_generations` stopping condition.

## 4. Implementation Approach

The implementation of Cycle 02 is a multi-stage process that focuses on building out the new intelligent engines and then weaving them into the existing workflow by refactoring the orchestrator to handle a cyclical, stateful process.

1.  **Configuration and Interface Expansion**: The first step is to update the foundational contracts. The `ActiveLearningConfig` Pydantic model will be added to `config/models.py`. The `ConfigExpander` will be updated to populate this new sub-model with sensible, physically-based defaults. Simultaneously, the new abstract base classes, `IExplorerSampler` and `ISimulationEngine`, will be added to `interfaces/engines.py`.
2.  **Implementation of the Explorer & Sampler Engine (Module B)**:
    -   **MACE Integration**: A new class, `ExplorerSampler`, will be created. It will use the `mace-torch` Python library to load the specified pre-trained universal model.
    -   **Exploratory MD**: The engine will use ASE's MD libraries (e.g., `VelocityVerlet`) to run a high-temperature simulation using the MACE model as the calculator. This will generate a long trajectory of atomic configurations.
    -   **Numba-Optimized Descriptors**: A standalone Python function will be written to calculate an atomic environment descriptor (e.g., SOAP). This function will be heavily decorated with `@numba.jit(nopython=True, parallel=True)` and will operate directly on NumPy arrays for maximum performance. This is the most performance-critical part of the module.
    -   **DIRECT Sampling**: The engine will call the Numba-optimized function to get descriptors for the entire trajectory. It will then use a standard clustering algorithm from a library like `scikit-learn` (e.g., a mini-batch k-means) to select a diverse subset of structures from the high-dimensional descriptor space.
3.  **Implementation of the Simulation Engine (Module E)**:
    -   **MD Backend**: A new class, `SimulationEngine`, will be created. It will use ASE's `VelocityVerlet` as the primary MD integrator for simplicity and direct Python integration.
    -   **Custom Calculator**: The engine will have logic to load the custom-trained ACE model file produced by the `TrainingEngine` and wrap it in an ASE-compatible `Calculator` object.
    -   **Uncertainty Quantification (UQ)**: A method will be implemented to calculate the model's uncertainty for a given `Atoms` object. If the ACE model is an ensemble, this method will calculate the standard deviation of the force predictions among the committee members. This UQ function will be called at a regular interval during the MD simulation.
    -   **Dynamic Thresholding**: A helper method will be implemented that takes the current training set, calculates the uncertainty for all structures in it, and returns a high percentile (e.g., the 95th) of that distribution.
    -   **Asynchronous Generator**: The main `run_simulation` method will be structured as an `async def` function that `yields` results. It will contain the main MD loop. Inside the loop, it will call the UQ function, compare the result to the dynamic threshold, and if it's exceeded, it will `yield` the current `AtomsModel`.
4.  **Refactoring the WorkflowOrchestrator**:
    -   A new public method, `run_active_learning_pipeline()`, will be created.
    -   This method will first call the necessary components (either `StructureGenerator` or `ExplorerSampler`, then `LabelingEngine`, then `TrainingEngine`) to produce the initial "Generation 0" model.
    -   The main active learning `for` loop will be implemented, iterating from 0 to `max_generations - 1`.
    -   Inside the loop, the orchestrator will `await` the results from the `simulation_engine.run_simulation` generator, iterating through the yielded structures.
    -   For each structure that is yielded, the orchestrator will call the `labeling_engine` and update the database. This process is sequential and blocking for each new structure.
    -   After the inner loop over the generator is complete, the orchestrator will call the `training_engine` with the full, augmented dataset to produce the next-generation model.
5.  **Updating the CLI**: The `cli.py` file will be modified to expose the new functionality. A new option or flag, such as `--active-learning`, will be added to the `run` command. The CLI will inspect this flag to decide whether to call the old `run_linear_pipeline` method or the new `run_active_learning_pipeline` method on the `WorkflowOrchestrator`.

## 5. Test Strategy

The testing strategy for Cycle 02 must be significantly more sophisticated than for Cycle 01, as it needs to cover not only the functionality of the new components but also the complex, stateful, and asynchronous interactions of the active learning feedback loop.

**Unit Testing Approach (Min 300 words):**
-   **Explorer & Sampler**: The expensive MACE model evaluation will be mocked to ensure tests are fast. The primary testing focus will be on the DIRECT sampler's logic. We will create a fixture that provides a large but deterministic, pre-computed trajectory (as a simple NumPy array). The test will then call the sampler with this trajectory and assert that the selection logic works as expected. For example, we will assert that it returns the correct number of structures and that the selected structures are geometrically distinct and far apart in the trajectory. A separate test, marked as a benchmark, will be created for the Numba-optimised descriptor function. This benchmark will run the function on a large input and assert that its execution time is below a certain threshold, ensuring our performance optimizations are effective.
-   **Simulation Engine**: This is the most challenging component to unit test effectively. A real MD simulation will not be run. Instead, the MD integrator will be mocked to yield a pre-defined, finite sequence of `AtomsModel` objects. The MLIP model's uncertainty calculation will also be mocked to return a corresponding sequence of pre-determined float values. For example, we can program the mocks to return four structures with corresponding uncertainties of `[0.1, 0.2, 5.0, 0.3]`. The test will then initialize the `SimulationEngine` with a fixed uncertainty threshold of `1.0`. We will then execute the `run_simulation` generator and assert that it *only* yields the third structure (the one with the `5.0` uncertainty). This isolates and verifies the correctness of the core uncertainty-triggering logic. A separate unit test will focus exclusively on the dynamic thresholding function, providing it with a list of uncertainty values and asserting that it correctly calculates the 95th percentile.
-   **Orchestrator Logic**: The new, cyclical logic in the `WorkflowOrchestrator` will be tested with all engine dependencies fully mocked. This allows us to test the high-level state management logic. We will configure the run for `max_generations: 3`. We will then run the pipeline and assert that the `training_engine.train` method is called exactly 3 times. We will also use spies to inspect the arguments passed between the mocked engines, for example, asserting that the model path returned by the `TrainingEngine` from the first generation is the same one that is passed to the `SimulationEngine` for the second generation.

**Integration Testing Approach (Min 300 words):**
-   **Full Active Learning Loop Integrity**: The most critical integration test for this cycle will verify the data flow through the entire feedback mechanism. To make this test fast and deterministic, we will use a real `Orchestrator` but with mocked versions of the most expensive and non-deterministic engines. The test will be configured as follows:
    1.  Use a mock `SimulationEngine` that is programmed to yield one specific, pre-defined `AtomsModel` and then immediately terminate its simulation.
    2.  Use a mock `LabelingEngine` that is programmed to return a specific, pre-defined `DFTResult` object only when it is passed the exact structure from the mock `SimulationEngine`.
    3.  Use a mock `TrainingEngine`.
    The test will then execute the active learning pipeline for a single generation. The core of the test will be a series of assertions that verify the sequence of events: we will assert that `simulation_engine.run_simulation` is called first, then we assert that `labeling_engine.run_calculation` is called with the *exact* structure object that the simulation engine yielded, and finally, we assert that `training_engine.train` is called with a dataset that now contains the *exact* result object from the labeling engine. This provides strong evidence that the data is flowing correctly and without corruption through the entire feedback loop.
-   **Configuration Propagation to Loop Behavior**: Another key integration test will ensure that the user's configuration choices in the new `ActiveLearningConfig` section are respected by the pipeline. This test will use a real `ConfigExpander` and a real `Orchestrator`, but with mocked engines. We will create an `input.yaml` that sets a non-default value, for example, `max_generations: 2`. The test will then run the full pipeline. The assertion will be simple: we will verify that the mock `TrainingEngine`'s `train` method was called exactly twice. This confirms that the configuration is not just being parsed correctly but is actively controlling the behavior of the main orchestrator loop, which is critical for user control over the system.
