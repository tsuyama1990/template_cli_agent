# SPEC.md: Cycle 2 - Advanced Simulation, Sampling, and UI

## 1. Summary

This document outlines the technical specifications for Cycle 2 of the MLIP-AutoPipe project. Building upon the stable core framework established in Cycle 1, this cycle is dedicated to implementing the advanced, high-value features that distinguish MLIP-AutoPipe as an intelligent and powerful tool for materials science research. The primary objective is to transform the basic pipeline into a sophisticated workflow capable of generating highly diverse and physically realistic training data for state-of-the-art Machine Learning Interatomic Potentials (MLIPs) like MACE.

The major technical advancements in this cycle are centered around the **Exploration Engine**. The simple EMT-based MD simulation from Cycle 1 will be replaced with a far more powerful and physically nuanced simulation engine. This includes integrating real MLIP models (e.g., MACE) as the calculator, which requires implementing a "late-binding" pattern to manage the significant memory footprint of these models in a parallel processing environment. Furthermore, the engine will be enhanced to perform **hybrid Molecular Dynamics / Monte Carlo (MD/MC)** simulations, introducing Monte Carlo moves like atom swaps to more effectively explore the compositional and configurational space of materials. A critical feature to be added is the **automatic thermodynamic ensemble switching**, where the engine intelligently detects the presence of a vacuum slab and chooses between NVT and NPT ensembles accordingly, a crucial step for realistic simulations of surfaces and interfaces.

In addition to the exploration capabilities, this cycle will introduce a more intelligent **Sampling Module**. The basic random sampler will be augmented with an implementation of **Farthest Point Sampling (FPS)**. This technique uses SOAP (Smooth Overlap of Atomic Positions) descriptors to select a subset of structures that are maximally diverse, which is known to produce more efficient and robust training datasets. The **Generator** module will also be expanded to include an `IonicGenerator` for creating charge-neutral ionic crystal structures.

Finally, to improve accessibility and user experience, Cycle 2 includes the development of a proof-of-concept **Web User Interface (Web UI)**. This will provide a graphical front-end for configuring and launching pipeline runs, making the tool accessible to a broader range of users who may not be comfortable with a command-line-only interface. The development will involve creating a FastAPI backend to orchestrate the pipeline and a simple, intuitive front-end. Comprehensive testing will accompany each new feature to ensure its correctness and integration into the existing framework.

## 2. System Architecture

The architecture in Cycle 2 expands upon the foundation from Cycle 1, introducing new components and enhancing existing ones. The modular design allows these new features to be integrated cleanly.

**File Structure (Cycle 2 Focus):**

Files to be created or modified are in **bold**. The focus is on enhancing the core library and adding a new `web/` application directory.

```
.
└── src/
    └── mlip_autopipec/
        ├── cli/
        │   └── main.py              # (Minor modifications for new configs)
        ├── core/
        │   ├── orchestrator.py      # (Modified to support new components)
        │   └── **models.py**            # (Expanded with new config options)
        ├── generators/
        │   ├── __init__.py
        │   ├── base.py
        │   ├── alloy.py
        │   └── **ionic.py**             # New IonicGenerator
        ├── exploration/
        │   └── **engine.py**            # Heavily modified for advanced features
        ├── sampling/
        │   ├── __init__.py
        │   ├── base.py
        │   ├── random_sampler.py
        │   └── **fps.py**               # New Farthest Point Sampler
        ├── storage/
        │   └── database.py
        └── **web/**
            ├── **__init__.py**
            ├── **main.py**              # FastAPI backend application
            └── **static/**              # (Directory for HTML/JS/CSS files)
                └── **index.html**
```

**Code Blueprints:**

*   **`core/models.py`**:
    *   The Pydantic models will be expanded.
    *   `ExplorationConfig` will be updated:
        *   `calculator: Literal["emt", "mace"]`: Add support for MACE.
        *   `mc_moves: bool = False`: Add a flag to enable/disable hybrid MD/MC.
        *   `mc_swap_probability: float = Field(0.1, ge=0, le=1)`: Probability for atom swap moves.
    *   `SamplingConfig` will be updated:
        *   `method: Literal["random", "fps"]`: Add support for FPS.
        *   `fps_n_features: int = 100`: Parameters specific to the FPS algorithm.
    *   A new `WebServerConfig` model might be added to configure the web server (e.g., host, port).

*   **`exploration/engine.py`**:
    *   The `ExplorationEngine` will be significantly refactored.
    *   **MLIP Integration:** It will include a `_get_calculator` method that dynamically imports and instantiates the MLIP model (e.g., `mace.calculators.mace_mp`). This method will be called *within* the worker process to avoid pickling large model objects. This is the "late-binding" pattern.
    *   **Hybrid MD/MC:** The main simulation loop within the worker process will be modified. On certain steps (controlled by a probability), instead of a standard MD step, it will perform a Monte Carlo move. This will involve writing a function to perform random atom swaps on the `ase.Atoms` object.
    *   **Auto Ensemble Switching:** A helper function `_detect_vacuum(atoms: Atoms) -> bool` will be implemented. This function will analyze the atomic density along each axis to determine if a vacuum slab exists. The main simulation setup logic will then use this function's output to choose between `ase.md.NVTBerendsen` (for NVT, if vacuum is detected) and `ase.md.NPTBerendsen` (for NPT, for bulk systems).
    *   **ZBL Mixing (Optional but Recommended):** The `_get_calculator` method will also be responsible for wrapping the MLIP calculator with ASE's `ComboCalculator` to mix in a ZBL potential for handling short-range repulsive interactions gracefully.

*   **`generators/ionic.py`**:
    *   A new `IonicGenerator` class inheriting from `BaseStructureGenerator` will be implemented.
    *   Its `generate()` method will take into account not just composition but also oxidation states of the elements to produce charge-neutral structures, which is critical for realistic ionic material simulations.

*   **`sampling/fps.py`**:
    *   A new `FPSSampler` class inheriting from `BaseSampler` will be created.
    *   Its `sample()` method will be more complex than the random sampler.
    *   It will first read all frames from the trajectories.
    *   For each frame, it will compute a feature vector, specifically the SOAP descriptor, using a library like `dscribe`.
    *   It will then implement the Farthest Point Sampling algorithm: iteratively select the data point that is farthest from the set of already-selected points, until the desired number of samples is reached.

*   **`web/main.py`**:
    *   This will be a new file containing a `FastAPI` application.
    *   It will have a main endpoint `/` that serves a static HTML file (`web/static/index.html`).
    *   It will have a POST endpoint, `/run_pipeline`, that accepts a JSON payload representing the pipeline configuration.
    *   This endpoint will parse the JSON into the `FullConfig` Pydantic model.
    *   It will then instantiate the `PipelineOrchestrator` and call its `run()` method as a background task using `fastapi.BackgroundTasks`. This is crucial to avoid blocking the HTTP request while the long-running simulation is in progress.
    *   It will immediately return a JSON response like `{"status": "started", "job_id": "..."}`.
    *   (Optional) A GET endpoint `/status/{job_id}` could be implemented to poll for the status of the run.

## 3. Design Architecture

The design in Cycle 2 focuses on enhancing modularity and introducing advanced scientific computing concepts. The schema-first approach with Pydantic continues to be central, providing the validation and structure needed to manage the increased complexity of the configuration.

**Pydantic Schema Evolution:**

The evolution of the Pydantic models is a key aspect of the design. The changes are designed to be backward-compatible where possible.

*   **`ExplorationConfig` Evolution:**
    *   The `calculator` field is changed from `Literal["emt"]` to `Literal["emt", "mace"]`. This is a non-breaking change for old configs that still use "emt".
    *   New fields like `mc_moves` are added with default values (`False`). This ensures that existing configuration files that do not contain this key will still parse correctly and default to the old behavior (no MC moves).
    *   **Consumer:** The `ExplorationEngine` is the primary consumer of this evolved model. It will now contain conditional logic based on these new fields (e.g., `if config.mc_moves: ...`).

*   **`SamplingConfig` Evolution:**
    *   Similar to the above, the `method` field is expanded to `Literal["random", "fps"]`.
    *   **Consumer:** The `PipelineOrchestrator` will consume this. It will contain a factory pattern: based on the value of `config.sampling.method`, it will decide whether to instantiate `RandomSampler` or the new `FPSSampler`. This demonstrates the extensibility of the architecture.

*   **`IonicGeneratorConfig` (Potential New Model):** A new model could be introduced to hold parameters specific to the `IonicGenerator`, such as a dictionary of oxidation states. This would then be nested within the `GenerationConfig`.

**Architectural Patterns:**

*   **Factory Pattern:** The `PipelineOrchestrator` will now act as a factory for the sampler component. This decouples the orchestrator from the concrete implementations of the samplers, making it easy to add new sampling methods in the future. `if method == 'random': return RandomSampler() elif method == 'fps': return FPSSampler()`.
*   **Strategy Pattern:** The exploration engine's ability to switch between MD and hybrid MD/MC, and between NVT and NPT ensembles, can be viewed as an implementation of the Strategy pattern. The specific "strategy" for the simulation is chosen at runtime based on the configuration and an analysis of the system's properties.
*   **Background Task Processing:** The Web UI's backend (`web/main.py`) will use this pattern to handle long-running simulation jobs. By offloading the `orchestrator.run()` call to a background task, the web server remains responsive and can immediately acknowledge the user's request.

## 4. Implementation Approach

The implementation of Cycle 2 will be tackled by feature, ensuring that each new major capability is built, tested, and integrated before moving to the next.

1.  **Evolve Pydantic Models:**
    *   First, update the models in `core/models.py` with the new fields and `Literal` options.
    *   Update the unit tests for the models to cover the new validation rules and default values.

2.  **Advanced Exploration Engine:** This is the most complex part of the cycle.
    *   **Refactor `ExplorationEngine`:** Modify the engine to support the "late-binding" calculator pattern. This involves moving the calculator instantiation logic into a separate function that will be called by the worker processes.
    *   **Integrate MACE:** Add the logic to the new calculator function to instantiate a MACE calculator if `config.exploration.calculator == "mace"`. This will require adding `mace-models` as a project dependency.
    *   **Implement Auto-Ensemble Logic:** Write and unit-test the `_detect_vacuum` function. The main engine logic will then be updated to call this function and select the appropriate ASE dynamics class.
    *   **Implement Hybrid MD/MC:** Add the conditional logic for performing Monte Carlo atom swaps within the simulation loop. Write a dedicated unit test for the swap logic itself.

3.  **Farthest Point Sampler:**
    *   Add `dscribe` as a dependency.
    *   Implement the `FPSSampler` class in `sampling/fps.py`.
    *   The `sample()` method will first have a loop to compute SOAP descriptors for all frames in the input trajectories.
    *   A second loop will implement the iterative FPS logic.
    *   Write a unit test with a small, predictable set of feature vectors to verify that the FPS algorithm selects the correct points.

4.  **Ionic Generator:**
    *   Implement the `IonicGenerator` in `generators/ionic.py`.
    *   Write unit tests to assert that the generated structures are charge-neutral based on a predefined set of oxidation states.

5.  **Update Orchestrator:**
    *   Modify the `PipelineOrchestrator` to implement the factory pattern for samplers. It should read the `config.sampling.method` and instantiate the correct sampler class.
    *   Update the orchestrator's unit tests to verify this new factory logic.

6.  **Web UI Development:**
    *   Create the `web/` directory.
    *   Implement the FastAPI backend in `web/main.py`. Test the API endpoints using an HTTP client library like `httpx` in the test suite. Ensure the `/run_pipeline` endpoint correctly starts a background task.
    *   Create a simple `index.html` file with a form for submitting the key configuration parameters. Add basic JavaScript to send the form data as a JSON payload to the `/run_pipeline` endpoint.

## 5. Test Strategy

The test strategy for Cycle 2 expands significantly to cover the new complex functionalities, with a continued emphasis on isolating logic in unit tests and verifying component interactions in integration tests.

**Unit Testing Approach (Min 300 words):**

The unit tests in Cycle 2 will target the new, advanced logic. For the **`ExplorationEngine`**, we will write specific tests for the new decision-making capabilities. A key test will focus on the automatic ensemble switching. We will create two `ase.Atoms` objects in the test: one representing a bulk material and another with a large vacuum gap. We will then pass these to a method in the engine and assert that it correctly returns "NPT" for the bulk system and "NVT" for the slab system. The hybrid MD/MC logic will be tested by creating a small atoms object, running a single step of the simulation with the MC probability set to 1.0, and asserting that the atom positions or species have changed in a way consistent with a swap move, rather than a standard dynamics step. These tests verify the *control flow* and *logic* without running expensive simulations.

For the **`FPSSampler`**, the unit test will be crucial for validating the algorithm's correctness. We will not test the SOAP descriptor generation itself (as that's `dscribe`'s job). Instead, we will create a mock trajectory with a predefined set of simple feature vectors (e.g., 2D vectors). For instance, we can create two clusters of points. The test will then run the FPS algorithm and assert that the selected indices correspond to one point from each cluster, demonstrating that the sampler prioritizes diversity. This validates the core FPS implementation.

The new **`IonicGenerator`** will be tested by providing it with a target composition (e.g., MgO) and oxidation states (`{'Mg': 2, 'O': -2}`). The test will then call `generate()` and assert on the returned `ase.Atoms` object that the total charge is zero. This confirms that the generator is correctly enforcing the charge neutrality constraint, which is its primary responsibility.

**Integration Testing Approach (Min 300 words):**

The integration tests for Cycle 2 will be more comprehensive and realistic. The main end-to-end CLI test from Cycle 1 will be extended to validate the new features. We will create new configuration files (`config_mace.yml`, `config_fps.yml`) specifically for these tests.

A new integration test will be created to validate the **MACE calculator integration**. This test will run the full pipeline using `config_mace.yml`. While it will still use a small system and short simulation time to be practical for CI, its primary purpose is to confirm that the "late-binding" calculator pattern works correctly in the `ProcessPoolExecutor` environment. A successful run without pickling errors or other multiprocessing-related crashes is the main success criterion. After the run, we will inspect the output database to ensure that energies and forces have been computed, proving the MLIP was called successfully.

Another integration test will focus on the **FPS sampler**. It will run the pipeline with `config_fps.yml`, which specifies `method: "fps"`. The test will be constructed with an initial structure that is known to produce a trajectory with low diversity. The goal is to verify that the `FPSSampler` is correctly integrated into the `PipelineOrchestrator`. After the run, while it's hard to assert the "diversity" of the output, we can at least verify that the pipeline ran to completion without error and produced a database with the correct number of samples. This confirms the orchestrator's factory logic and the sampler's integration.

Finally, the **Web UI backend** will have its own set of integration tests. These tests will use a library like `httpx` to make live API calls to the FastAPI test server. A test will POST a valid configuration to the `/run_pipeline` endpoint, assert that it receives a `200 OK` response with a job ID, and verify that a (mocked) `PipelineOrchestrator` was instantiated and its `run` method was called in the background. This confirms that the web interface can correctly trigger the core application logic.
