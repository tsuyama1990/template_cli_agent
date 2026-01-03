# CYCLE 02: SPECIFICATION - Advanced Exploration, Sampling, and Web UI

## 1. Summary

This document provides the detailed technical specification for Cycle 2 of the MLIP-AutoPipe project. Building upon the foundational command-line pipeline established in Cycle 1, this cycle introduces the advanced scientific features and the user-friendly graphical interface that distinguish MLIP-AutoPipe as a powerful, expert-level tool. The primary focus is on transforming the simplified exploration engine into a sophisticated, hybrid Molecular Dynamics / Monte Carlo (MD/MC) engine capable of more effective and diverse sampling of the potential energy surface. This includes the implementation of critical physics-based logic such as automatic thermodynamic ensemble switching and the mixing of MLIPs with classical potentials to prevent simulation artifacts.

Furthermore, this cycle will enhance the data selection process by implementing a Farthest Point Sampling (FPS) strategy, moving beyond the basic random sampling of Cycle 1 to ensure that the final dataset is maximally diverse and information-rich. Support for a wider range of material systems, including ionic and covalent crystals, will be added by creating new generator modules. The culminating feature of this cycle will be the development of a Web-based User Interface (Web UI). This interface will provide an intuitive, interactive way for users to configure and launch pipeline runs, monitor their progress, and visualize the generated atomic structures directly in their browser. This will significantly lower the barrier to entry for the software and broaden its user base.

## 2. System Architecture

Cycle 2 involves enhancing existing modules with more complex logic and adding a new top-level component for the Web UI. The core architecture remains the same, demonstrating its extensibility.

**File Structure (Cycle 2 Focus):**
The files and directories to be created or modified in this cycle are marked in **bold**. Existing files from Cycle 1 that will be modified are also marked.

```
src/mlip_autopipec/
├── cli/
│   └── main.py              # No major changes, will benefit from new features
├── web/
│   └── **app.py**               # Web UI entry point (FastAPI)
│   └── **templates/**           # HTML templates for the UI
│       └── **index.html**
│   └── **static/**              # CSS/JS files
├── pipeline/
│   └── orchestrator.py      # Minor changes to support new generator/sampler types
│   └── interfaces.py        # Add IStructureSampler interface
│   └── factories.py         # Update to include new components
├── generators/
│   ├── factory.py           # **Modified** to select new generator types
│   └── **ionic.py**             # IonicGenerator implementation
│   └── **covalent.py**          # CovalentGenerator implementation
├── exploration/
│   └── **engine.py**            # **Heavily Modified** to implement advanced MD/MC logic
│   └── **physics.py**           # **New file** for vacuum detection and other physics logic
├── sampling/
│   ├── **__init__.py**
│   ├── **base.py**              # **New file**: Abstract BaseSampler class
│   ├── **factory.py**           # **New file**: SamplerFactory
│   ├── **fps.py**               # **New file**: Farthest Point Sampler implementation
│   └── **random.py**            # **Refactored** from orchestrator into its own file
├── storage/
│   └── database.py          # No major changes
├── config/
│   └── **schema.py**            # **Modified** to add config options for new features
└── common/
    └── atoms_utils.py       # May be expanded with new utilities
tests/
├── unit/
│   ├── **test_exploration_engine.py** # New tests for advanced logic
│   ├── **test_physics.py**
│   ├── **test_fps_sampler.py**
│   └── **test_web_app.py**
└── integration/
    └── **test_web_pipeline.py**
```

**Architectural Blueprint:**

1.  **Advanced Exploration Engine (`exploration/engine.py` & `physics.py`):**
    *   The `LabelingEngine` will be significantly upgraded. It will now manage a more complex simulation state.
    *   **Hybrid MD/MC Logic:** The `run` method will be modified to include a loop that alternates between MD steps and MC move attempts (e.g., one MC attempt every 10 MD steps).
    *   **MC Moves:** Internal methods like `_attempt_atom_swap()` and `_attempt_vacancy_hop()` will be implemented. These methods will respect physical constraints, such as a `ChargeSafety` check to prevent swapping of ions with different charges in ionic crystals.
    *   **Auto Ensemble Switching:** A new helper module, `physics.py`, will be created. It will contain a `detect_vacuum` function that uses a grid-based algorithm to determine if an `Atoms` object represents a bulk material or a slab/surface. Before starting a simulation, the `LabelingEngine` will call this function. If a vacuum is detected, it will use an NVT ensemble; otherwise, it will use NPT. This prevents artificial compression/expansion of the simulation box's vacuum layer.
    *   **ZBL Potential Mixing:** The `_get_calculator` internal method will be enhanced. It will load the primary MLIP calculator (e.g., MACE) and, if configured, will wrap it in `ase.calculators.mixing.MixedCalculator`. The second calculator in the mix will be a purely repulsive potential like ZBL, which will handle the physics of atoms getting extremely close, preventing the MLIP from having to learn these rare, high-energy events and making the simulation more stable.

2.  **Farthest Point Sampler (`sampling/`):**
    *   A new `sampling` module will be created.
    *   `base.py` will define the `IStructureSampler` interface with a `sample(atoms_list: list[Atoms]) -> list[Atoms]` method.
    *   `random.py` will encapsulate the simple random sampling logic.
    *   `fps.py` will contain the `FarthestPointSampler`. Its `sample` method will first compute a SOAP descriptor for every `Atoms` object in the input list (using a library like `dscribe`). It will then implement the iterative FPS algorithm to select a subset of structures that are maximally distant from each other in this descriptor space.
    *   The `WorkflowOrchestrator` will be modified. The sampling logic will be removed from its main loop and delegated to a sampler instance created by a new `SamplerFactory`.

3.  **Expanded Generators (`generators/`):**
    *   New classes `IonicGenerator` and `CovalentGenerator` will be created, inheriting from `BaseStructureGenerator`. They will contain the specific logic for creating their respective crystal structures, including handling charge balance for ionic systems.
    *   The `GeneratorFactory` will be updated to recognize new `system_type` values in the configuration and instantiate the correct generator class.

4.  **Web User Interface (`web/app.py`):**
    *   A new `web` module will be created.
    *   `app.py` will be a FastAPI application. It will serve a main HTML page and provide a RESTful API.
    *   **API Endpoints:**
        *   `GET /`: Serves the main `index.html`.
        *   `POST /run_pipeline`: Accepts a JSON payload containing the pipeline configuration. This endpoint will validate the JSON using the `FullConfig` Pydantic model. It will then start the `WorkflowOrchestrator` in a background thread or process and immediately return a job ID to the user.
        *   `GET /status/{job_id}`: A polling endpoint that the frontend can use to get the current status of a running job (e.g., "Generating", "Exploring: 25% complete", "Finished").
        *   `GET /results/{job_id}`: An endpoint that returns the final results from the database, perhaps as a list of JSON-serialized `Atoms` objects.
    *   **Frontend (`templates/index.html`):** A single-page application using modern HTML, CSS, and vanilla JavaScript. It will feature:
        *   A form with input fields that correspond to the `FullConfig` Pydantic model. JavaScript will provide client-side validation for a better user experience.
        *   An "Run Pipeline" button that sends the form data to the `POST /run_pipeline` endpoint.
        *   A status display area that polls the `/status/{job_id}` endpoint to show the live progress of the run.
        *   A 3D molecular viewer (e.g., `nglview` embedded or a similar JS library like `3Dmol.js`) in a results area. When the job is finished, the frontend will fetch data from `/results/{job_id}` and display the generated structures in the viewer.

## 3. Design Architecture

The design for Cycle 2 continues the schema-first, loosely-coupled philosophy.

**Pydantic-based Schema Design:**

The `config.schema.FullConfig` model will be extended to support the new features.

```python
# In src/mlip_autopipec/config/schema.py

class SystemType(str, Enum):
    ALLOY = "alloy"
    IONIC = "ionic"
    COVALENT = "covalent" # New

class MCConfig(BaseModel):
    swap_probability: float = Field(0.1, ge=0, le=1)
    # ... other MC parameters

class ExplorationConfig(BaseModel):
    # ... existing fields
    use_mc: bool = False
    mc_config: MCConfig | None = None
    use_zbl_mixing: bool = True

class SamplingMethod(str, Enum):
    RANDOM = "random"
    FPS = "fps" # New

class FPSSamplerConfig(BaseModel):
    soap_descriptor_config: dict # Configuration for the SOAP descriptor

class SamplingConfig(BaseModel):
    method: SamplingMethod
    n_samples: int = Field(..., gt=0)
    fps_config: FPSSamplerConfig | None = None

class FullConfig(BaseModel):
    # ... existing fields
    # The nested models are updated
```

*   **Producers/Consumers:** The Web UI's frontend JavaScript will now also be a producer of this configuration data (in JSON format). The FastAPI backend (`web/app.py`) will be the consumer, validating the incoming JSON.
*   **Constraints:** The use of `Union` types and optional models (e.g., `mc_config: MCConfig | None = None`) allows for a flexible configuration where settings for a specific feature are only required if that feature is enabled (`use_mc: True`). Pydantic's validators will be used to enforce these co-dependencies (e.g., if `method` is `FPS`, then `fps_config` must be provided).

**Decoupling the Web UI from the Core Logic:**

The Web UI (`web/app.py`) will **not** contain any business logic. It is purely a presentation and API layer. It will interact with the core pipeline in the same way the CLI does: by importing and using the `WorkflowOrchestrator`. To avoid blocking the web server, the orchestrator will be run in a separate thread or process managed by a task queue (like Celery, or more simply, FastAPI's `BackgroundTasks`).

*   `web/app.py` -> (uses) `pipeline/orchestrator.py`
*   `cli/main.py` -> (uses) `pipeline/orchestrator.py`

This ensures that the core logic is developed and tested once and used by multiple interfaces, adhering to the DRY (Don't Repeat Yourself) principle.

## 4. Implementation Approach

1.  **Update Configuration (`config/schema.py`):** First, extend the Pydantic models to include all the new configuration options for the advanced MD/MC, FPS sampling, and new generator types.
2.  **Implement Physics Helpers (`exploration/physics.py`):** Create the new `physics.py` module. Implement and unit-test the `detect_vacuum` function in isolation.
3.  **Upgrade Exploration Engine (`exploration/engine.py`):** This is the most complex step.
    a. Refactor the `run` method to incorporate the main MD/MC loop.
    b. Add the auto-ensemble switching logic, calling `detect_vacuum`.
    c. Add the ZBL potential mixing logic.
    d. Write extensive unit tests for these new, complex behaviors. Use mock objects for calculators and fixed random seeds to test the MC logic deterministically.
4.  **Implement Sampler Module (`sampling/`):**
    a. Create the `IStructureSampler` interface and the `sampling` directory.
    b. Refactor the existing random sampling logic into `random.py`.
    c. Implement the `FarthestPointSampler` in `fps.py`. This may require adding a new dependency (like `dscribe`) to `pyproject.toml`. Unit test it with a small, fixed dataset.
    d. Create the `SamplerFactory`.
5.  **Update Orchestrator (`pipeline/orchestrator.py`):**
    a. Modify the `WorkflowOrchestrator` to remove the old sampling logic.
    b. In its place, it will now accept a sampler from the factory in its constructor and call it at the appropriate point in the workflow.
6.  **Implement New Generators (`generators/`):** Create `IonicGenerator` and `CovalentGenerator`. Update the `GeneratorFactory`.
7.  **Develop Web UI Backend (`web/app.py`):**
    a. Set up a basic FastAPI application.
    b. Create the Pydantic models for API requests and responses.
    c. Implement the `/run_pipeline`, `/status/{job_id}`, and `/results/{job_id}` endpoints. Use FastAPI's `BackgroundTasks` to run the orchestrator asynchronously.
8.  **Develop Web UI Frontend (`web/templates/index.html`):**
    a. Create the HTML structure for the configuration form.
    b. Write JavaScript to handle form submission (using `fetch` to call the backend API), status polling (`setInterval`), and results visualization (integrating with a 3D viewer library).
9.  **Integration Testing:**
    a. Create a new integration test file (`test_web_pipeline.py`) that uses FastAPI's `TestClient`.
    b. The test will simulate a frontend by making HTTP requests to the backend endpoints and asserting that the responses are correct and that the pipeline runs successfully in the background.

## 5. Test Strategy

**Unit Testing Approach (Min 300 words):**
Unit tests for Cycle 2 will focus on the new, complex logic modules.
*   **Physics Module (`test_physics.py`):** The `detect_vacuum` function will be tested with carefully constructed `ase.Atoms` objects. One object will represent a bulk crystal (e.g., face-centered cubic), and the test will assert `detect_vacuum` returns `False`. Another object will have a large empty space added in the z-direction, and the test will assert the function returns `True`.
*   **Exploration Engine (`test_exploration_engine.py`):** Testing the engine's new logic is critical. We will mock the ASE calculators and MD runners. To test auto-ensemble switching, we will patch `exploration.engine.detect_vacuum` to return `True` or `False` and assert that the correct ASE dynamics object (NVT or NPT) is instantiated and used. For MC moves, we will test the acceptance/rejection criteria. For example, we can provide an `Atoms` object, mock the random number generator to produce a value that guarantees acceptance, and assert that the returned structure has two atoms swapped.
*   **FPS Sampler (`test_fps_sampler.py`):** We will create a list of 3 mock `Atoms` objects. We will then patch the SOAP descriptor calculation function to return a fixed set of vectors, for example, `[[0, 0], [1, 0], [10, 0]]`. If we ask the `FarthestPointSampler` to select 2 samples, the test must assert that it selects the first and third objects, as they are farthest apart.
*   **Web Backend (`test_web_app.py`):** Using FastAPI's `TestClient`, we will make mock HTTP requests. For the `/run_pipeline` endpoint, we will patch the `WorkflowOrchestrator` so it doesn't actually run. The test will send a valid JSON config and assert a `200 OK` response with a job ID. It will also send invalid JSON and assert a `422 Unprocessable Entity` response with a detailed error message from the Pydantic validation layer.

**Integration Testing Approach (Min 300 words):**
Integration tests will verify that the new and modified components work together and with the existing pipeline.
*   **Advanced CLI Pipeline (`test_cli_pipeline.py`):** The existing CLI integration test from Cycle 1 will be expanded. We will create a new configuration file that enables the advanced features, for example, `use_mc: true` and `sampling: {method: fps}`. The test will run the full pipeline via the CLI and perform assertions on the final database. For example, it will check that the number of atoms of each species is conserved (verifying the swap logic) and that the final number of structures in the database matches the `n_samples` requested by the FPS sampler.
*   **Web UI End-to-End Test (`test_web_pipeline.py`):** This will be a "headless" integration test of the web interface. It will not use a real browser but will interact with the application at the HTTP level using the `TestClient`. The test will:
    1.  POST a valid configuration to `/run_pipeline`.
    2.  Assert that it gets a job ID back.
    3.  Enter a loop that calls the `/status/{job_id}` endpoint, simulating the frontend polling. The test will assert that the status changes from "Running" to "Finished".
    4.  Once finished, it will call the `/results/{job_id}` endpoint and assert that it receives a valid JSON response containing the structural data.
    5.  Finally, it will inspect the filesystem to ensure the database file was actually created and populated correctly. This provides end-to-end confidence, from the user's HTTP request through the background processing to the final data persistence.
