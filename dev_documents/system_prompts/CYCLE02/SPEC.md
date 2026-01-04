# SPECIFICATION: Cycle 2 - Advanced Features and Web UI

## 1. Summary

This document provides the detailed technical specification for the second development cycle of the MLIP-AutoPipe project. Cycle 2 represents a significant leap forward in both scientific sophistication and user accessibility, building upon the solid foundation of the core CLI pipeline established in Cycle 1. The primary objectives of this cycle are twofold. First, it aims to substantially upgrade the exploration engine with more advanced, physically-motivated simulation techniques, transforming it from a simple MD runner into an intelligent phase space explorer. Second, it will deliver a rich, interactive web-based user interface (Web UI), making the full power of the pipeline accessible to a broader audience without requiring command-line expertise. This cycle is about moving beyond a functional MVP to create a truly powerful, intelligent, and user-friendly scientific tool.

The key technical enhancements are centered on the `MDMCExplorer` service. A headline feature is the implementation of a **hybrid MD/MC (Molecular Dynamics / Monte Carlo)** algorithm. While MD is excellent at exploring the local potential energy surface, it can struggle to overcome high energy barriers. The introduction of MC "moves," such as atomic swaps and vacancy hops, allows the simulation to make non-physical but chemically plausible jumps, which is a highly effective technique for efficiently exploring the vast configurational and compositional space of complex materials like multi-component alloys. A second major innovation is the development of an **automatic thermodynamic ensemble switching** mechanism. This intelligent feature will analyze the geometry of the simulation cell to detect the presence of a vacuum layer (a characteristic of surface or slab simulations). Based on this detection, it will automatically switch the simulation from the NPT (constant pressure) to the NVT (constant volume) ensemble, a critical step to prevent the unphysical collapse of the vacuum layer, thereby ensuring the physical realism of surface simulations. Furthermore, a **mixed potential model combining the MLIP with the classical ZBL potential** will be integrated. This addresses a common failure mode in high-temperature simulations where atoms can get too close, leading to extremely high repulsive forces that MLIPs often struggle to model. The ZBL potential provides a robust, physically-correct repulsive wall at short distances, dramatically improving the stability of the simulations. Finally, this cycle re-introduces the **Farthest Point Sampling (FPS)** algorithm, a more intelligent method for curating the final dataset to maximize structural diversity.

On the user-facing side, Cycle 2 will deliver a complete, browser-based graphical user interface (Web UI). This UI will serve as an alternative, interactive front-end to the entire pipeline. It will guide users through the process of building and validating simulation configurations via an intuitive web form, eliminating the need to manually edit YAML files. Users will be able to launch runs, monitor their progress in real-time through a live log viewer, and, once a run is complete, browse and visualize the results. The results page will feature an integrated 3D structure viewer, allowing users to inspect the generated atomic configurations interactively. This will dramatically lower the barrier to entry for new or less technical users and provide a powerful tool for visual analysis and interactive exploration for experts.

## 2. System Architecture

The architecture in Cycle 2 expands significantly to accommodate the new Web UI and the enhanced backend services. A new presentation layer entry point, `web_ui.py`, is introduced, which will run a web server and interact with the existing core pipeline. The `MDMCExplorer` service and the `domain/models.py` file will undergo substantial modifications to incorporate the new, advanced features. However, a key architectural principle is to reuse the existing `PipelineOrchestrator` from Cycle 1 without modification. The Web UI will simply be another client that prepares a `FullConfig` object and passes it to the orchestrator, demonstrating the power of the modular design. This ensures that the core logic remains consistent and that all the robustness and testing from Cycle 1 are carried forward.

**File Structure (Cycle 2 Focus):**

Bold entries indicate new or significantly modified files for this cycle. The structure now includes top-level directories for static assets and HTML templates, standard for a web application.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py                        # (Existing)
├── **web_ui.py**                 # **NEW**: Main Web UI application, likely using FastAPI for its performance and Pydantic integration.
│
├── core/                         # (Existing)
│   ├── __init__.py
│   ├── pipeline_orchestrator.py
│   └── factories.py
│
├── domain/                       # (Modified)
│   ├── __init__.py
│   ├── **models.py**             # **MODIFIED**: Pydantic models will be extended to include configuration options for all new features.
│   └── interfaces.py
│
├── infrastructure/               # (Existing)
│   ├── __init__.py
│   ├── ase_db_wrapper.py
│   └── process_runner.py
│
└── services/                     # (Modified)
    ├── __init__.py
    ├── generation/
    │   ├── ... (existing generators)
    │   └── **specialized_generators.py** # **NEW**: Home for new generators like InterfaceGenerator and AdsorptionGenerator.
    ├── exploration/
    │   ├── __init__.py
    │   └── **md_mc_explorer.py**     # **MODIFIED**: Will receive major enhancements for hybrid MD/MC, auto-ensemble switching, and ZBL mixing.
    └── sampling/
        ├── ... (existing samplers)
        └── **fps_sampler.py**          # **NEW**: The intelligent Farthest Point Sampler will be implemented here.

**static/**                  # **NEW**: Directory to serve static files like CSS and JavaScript for the Web UI.
**templates/**               # **NEW**: Directory to store HTML templates (e.g., using Jinja2) for the Web UI.
```

**Component Blueprint:**

This blueprint details the responsibilities of the new and modified components for Cycle 2.

*   **`web_ui.py`**: This new file will contain a complete web application built using a modern Python framework like FastAPI, chosen for its excellent performance and native integration with Pydantic, which will simplify API validation. It will define several API endpoints:
    *   A main `/` endpoint that serves the HTML front-end, a single-page application that allows users to construct a `FullConfig` object through an interactive web form.
    *   An API endpoint `/api/run` that accepts the configuration as a JSON object. This endpoint will validate the incoming data against the Pydantic models and then launch the `PipelineOrchestrator` in a background process (e.g., using `multiprocessing` or a task queue) to avoid blocking the web server. It will return a unique `run_id` to the client.
    *   An API endpoint `/api/status/<run_id>` that the frontend can poll to get real-time progress updates on a specific pipeline run.
    *   API endpoints for the results page, allowing it to fetch a list of structures from a completed run's database and retrieve the specific coordinates for visualization.

*   **`services/exploration/md_mc_explorer.py`**: This existing file will be heavily modified to incorporate the advanced scientific logic.
    *   The core simulation method will be refactored to support a hybrid MD/MC loop. It will be configurable to perform MC moves (like atom swaps) at a specified frequency, interspersed with standard MD integration steps.
    *   A new private method, `_detect_vacuum()`, will be implemented. This function will be a key piece of new IP. It will analyze the geometry of an `ase.Atoms` object by creating a 3D grid over the cell and checking for large, contiguous regions with no atoms, which robustly identifies a vacuum slab.
    *   The main simulation logic will be updated to call `_detect_vacuum()` at the beginning of a run. If a vacuum is detected, it will force the simulation to use the NVT ensemble, overriding any user setting to prevent simulation artifacts.
    *   The calculator setup logic will be enhanced to create a mixed potential using `ase.calculators.mixing.MixedCalculator`. This will combine the forces from the primary MLIP with the forces from an `ase.calculators.zbl.ZBL` potential for short-range repulsion.

*   **`domain/models.py`**: The Pydantic models will be updated to expose the new features to the user.
    *   The `ExplorationConfig` model will be extended with new, well-documented fields: `enable_hybrid_mc: bool`, `mc_swap_frequency: int`, `auto_ensemble_switching: bool` (defaulting to `True`), and `use_zbl_repulsion: bool` (defaulting to `True`).
    *   The `SamplingConfig` model will be updated to allow the `method` field to be `Literal['random', 'fps']` and will include a new optional field `fps_soap_params: dict | None` for configuring the SOAP descriptors. A validator will be added to ensure that if `method` is 'fps', then `fps_soap_params` is not `None`.

*   **`services/sampling/fps_sampler.py`**: This new file will contain the implementation of the `FPSSampler`. It will depend on an external library (like `dscribe`) to compute the SOAP (Smooth Overlap of Atomic Positions) descriptors for each structure in the input trajectories. It will then implement the iterative Farthest Point Sampling algorithm to select a subset of structures that are maximally diverse in this high-dimensional descriptor space, providing a much more efficient sampling of the configurational landscape than the random approach.

## 3. Design Architecture

The design for Cycle 2 is an extension of the clean, modular architecture established in Cycle 1. The introduction of the Web UI is handled by adding a new component to the presentation layer, which communicates with the existing backend through a well-defined API. This maintains the strict separation of concerns and demonstrates the extensibility of the initial design. The core pipeline logic within the `PipelineOrchestrator` remains unchanged, acting as a stable foundation upon which these new features are built.

*   **Pydantic Model Enhancements and Backward Compatibility**: The new configuration fields added to `ExplorationConfig` and `SamplingConfig` in `domain/models.py` will be designed with default values (e.g., `enable_hybrid_mc: bool = False`). This is a crucial design choice that ensures **backward compatibility**. Configuration files written for Cycle 1, which do not contain these new keys, will still be valid and will run with the default (Cycle 1) behavior. Custom validators will be added for the new fields to ensure logical consistency. For instance, a validator will check that if `enable_hybrid_mc` is `True`, then the `mc_swap_frequency` must be a positive integer greater than zero.

*   **Web UI and Backend Interaction (Stateless API Design)**: The interaction between the frontend and backend will be designed around a stateless RESTful API. The Web UI will be a modern, single-page application (SPA) built with HTML, CSS, and JavaScript. All dynamic behavior will be handled by the frontend, which will make asynchronous API calls to the Python backend. The backend (`web_ui.py`) will be designed to be stateless. When it receives a request on the `/api/run` endpoint, it will validate the input, start the long-running `PipelineOrchestrator` in a separate background process, and immediately return a `run_id`. The state of the run (e.g., "running", "completed", "failed") will be tracked via simple status files on the disk, not in the memory of the web server. The `/api/status/<run_id>` endpoint will simply read these files to provide updates. This decoupled, stateless design is highly robust and scalable; it ensures that the web server remains responsive and is not tied up by the long-running scientific computations.

*   **Robust Vacuum Detection Algorithm**: The `_detect_vacuum` function will be a critical piece of new scientific logic. Its design will be more sophisticated than a simple check of the cell's aspect ratio. The proposed algorithm is as follows:
    1.  Input: An `ase.Atoms` object.
    2.  Grid Creation: Define a 3D grid of points (e.g., with a spacing of 1 Å) that spans the entire simulation cell.
    3.  Distance Calculation: For each point on the grid, calculate the distance to the nearest atom in the `Atoms` object.
    4.  Vacuum Identification: A grid point is considered to be "in a vacuum" if its distance to the nearest atom is greater than a specified threshold (e.g., a value slightly larger than a typical covalent bond length, ~3-4 Å).
    5.  Slab Detection: The algorithm will then check if there is a contiguous *plane* of these vacuum points along any of the three Cartesian axes (x, y, or z). The presence of such a plane is a robust indicator of a vacuum slab.
    6.  Output: The function will return `True` if a vacuum slab is detected, and `False` otherwise.
    This grid-based approach is highly robust and can correctly identify vacuum layers in complex geometries, such as those with tilted slabs or internal pores, where simpler methods would fail.

## 4. Implementation Approach

The implementation of Cycle 2 will be strategically divided into two parallel development tracks to allow for concurrent progress on the backend scientific features and the frontend user interface. This requires clear API definitions from the start to ensure smooth integration.

1.  **Track 1: Backend Scientific Enhancements**: This track focuses on upgrading the `MDMCExplorer` and adding the `FPSSampler`.
    *   **Step 1.1 (Models)**: The first step is to modify the Pydantic models in `domain/models.py`. The new configuration parameters for hybrid MC, ZBL repulsion, auto-ensemble switching, and FPS will be added, complete with default values and validators.
    *   **Step 1.2 (Services)**: Implement the new scientific logic within the services. This includes the `_detect_vacuum` function, the logic for creating a mixed ZBL/MLIP potential in the `MDMCExplorer`, and the implementation of the hybrid MD/MC loop. Simultaneously, the `FPSSampler` will be developed, including the integration of a third-party library for SOAP descriptor calculations.
    *   **Step 1.3 (Unit Testing)**: Each new piece of functionality will be accompanied by rigorous unit tests. For instance, the `_detect_vacuum` function will be tested with a variety of pre-defined bulk and slab structures to ensure its accuracy. The potential mixing logic will be tested to verify that the resulting forces are correct.

2.  **Track 2: Web UI Development**: This track focuses on creating the new user-facing application.
    *   **Step 2.1 (Web Server Setup)**: Set up the basic web server application in `web_ui.py` using FastAPI. Define the necessary HTML templates and create the basic API endpoints (`/api/run`, `/api/status`, `/api/results`). Initially, these endpoints can return mock data.
    *   **Step 2.2 (Frontend Development)**: Develop the frontend as a single-page application. This involves writing the HTML for the main configuration form, the run status page, and the results browser. JavaScript code will be written to handle user interactions, perform client-side validation, make AJAX calls to the backend API, and dynamically update the UI based on the responses.
    *   **Step 2.3 (3D Viewer Integration)**: Integrate a lightweight, browser-based 3D structure viewer (such as `nglview` or a similar JavaScript library) into the results page. The JavaScript will be responsible for fetching structure data from a backend API endpoint and rendering it in the viewer.

3.  **Track 3: Integration and End-to-End Testing**: Once significant progress has been made on both tracks, the final phase is to integrate them and perform comprehensive end-to-end testing.
    *   The mock API endpoints in the web server will be replaced with calls to the actual `PipelineOrchestrator`, connecting the UI to the backend engine.
    *   A suite of end-to-end tests will be developed using a browser automation tool like Playwright. These automated tests will script the entire user workflow: opening the web page, filling out the form to enable the new advanced features, launching a run, polling for its completion, and verifying that the results are displayed correctly. These tests are the ultimate guarantee that the full stack is working seamlessly.

## 5. Test Strategy

The testing strategy for Cycle 2 must expand significantly to cover the increased complexity of the backend and the entirely new Web UI stack. It will build upon the foundation of tests from Cycle 1.

**Unit Testing Approach (Min 300 words):**

The unit tests for Cycle 2 will be highly targeted, focusing on the new, complex pieces of logic to ensure they are correct in isolation.

*   **`MDMCExplorer` Advanced Features**: This service will receive the most extensive new set of unit tests.
    *   **Vacuum Detection**: A dedicated test suite will be created for the `_detect_vacuum` function. We will construct several static `ase.Atoms` objects in the test code that represent clear-cut physical cases: a perfect bulk FCC crystal, a surface slab with a large vacuum gap along the z-axis, a molecule placed in the center of a large, empty box, and perhaps a porous material. We will then write a separate test for each case, asserting that the function returns the expected boolean value (`False` for bulk, `True` for the slab and molecule).
    *   **Ensemble Switching Logic**: We will test the logic that *uses* the vacuum detection. This test will use `pytest-mock` to mock the `_detect_vacuum` method, forcing it to return `True` or `False`. We will then inspect the `ase.md` dynamics object that gets created and assert that it is of the correct type (`NVTBerendsen` when the mock returns `True`, and `NPTBerendsen` when it returns `False`).
    *   **Hybrid MC Logic**: The MC move proposal functions (e.g., `_propose_swap`) will be tested in isolation. A test will create a simple `Atoms` object with two different element types, call the swap function, and then assert that the atomic numbers of two specific atoms have been correctly interchanged.
    *   **ZBL Mixing**: The potential mixing logic will be tested by creating a simple two-atom system. We will calculate the forces using only the MLIP and only the ZBL potential. We will then calculate the forces using the mixed potential and assert that the result is the correct vector sum of the two individual force calculations.

*   **`FPSSampler`**: The `FPSSampler` will be tested with a carefully constructed dummy trajectory that contains clusters of nearly identical structures and a few unique, distinct structures. The test will mock the SOAP descriptor calculation to return deterministic, pre-computed fingerprints. The core assertion will be that the sampler correctly identifies and selects the set of unique structures, demonstrating that the diversity maximization logic is working.

*   **Web UI Backend (`web_ui.py`)**: We will use FastAPI's `TestClient` to perform in-memory unit testing of our API endpoints without needing to run a live web server. We will send POST requests with valid and invalid JSON payloads to the `/api/run` endpoint and assert that the server responds with the correct HTTP status codes (200 for valid, 422 for invalid). We will mock the `PipelineOrchestrator` to verify that it is called with the correct configuration when the API receives a valid request.

**End-to-End (E2E) Testing Approach (Min 300 words):**

With the introduction of a browser-based UI, E2E testing becomes absolutely critical. It is the only way to verify that the entire application stack, from the user's clicks in the browser down to the scientific calculations in the backend, is working together seamlessly. We will use the `pytest-playwright` framework for this purpose, which allows us to write our browser automation tests in Python.

*   **The "Golden Path" E2E Test**: The most important E2E test will be a "golden path" scenario that simulates a complete, successful user workflow through the Web UI. The automated test script will perform the following actions:
    1.  **Setup**: The test will start the `web_ui.py` application as a background process.
    2.  **Navigation**: Playwright will launch a real browser (e.g., Chromium), navigate to the application's local URL (e.g., `http://127.0.0.1:8000`).
    3.  **Form Interaction**: The script will programmatically interact with the web form's HTML elements. It will fill in the input fields for a simple Si-Ge alloy system and will be sure to click the checkboxes to enable the new advanced features like "Hybrid MC" and "ZBL Repulsion," and select "FPS" as the sampling method.
    4.  **Submission**: It will simulate a user click on the "Launch Run" button.
    5.  **Status Verification**: The script will then assert that the UI correctly transitions to the "Run in Progress" page and that a `run_id` is displayed. To avoid brittle tests that rely on UI timings, the script will then switch to polling the `/api/status/<run_id>` API endpoint directly in a loop, waiting for the JSON response to indicate that the status has changed from "Running" to "Completed".
    6.  **Results Verification**: Once the run is complete, the script will navigate the browser to the results page for that `run_id`. It will then assert that the page correctly displays the results, for example, by checking that the results table is present and contains the expected number of structures.
    7.  **Visualization Check**: As a final check, the script will simulate a click on one of the structures in the results table and assert that the 3D viewer element is successfully loaded into the page.

This single, comprehensive E2E test provides a very high degree of confidence that the integration between the frontend JavaScript, the web backend API, the pipeline orchestration logic, and the new advanced simulation features is correct and robust.
