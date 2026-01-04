# Specification: MLIP-AutoPipe Cycle 2

## 1. Summary

This document provides the detailed technical specification for the second development cycle of the MLIP-AutoPipe project. Building upon the solid foundation of the core command-line pipeline established in Cycle 1, Cycle 2 focuses on two primary objectives: significantly enhancing the **algorithmic sophistication** of the pipeline and dramatically improving **user accessibility** through the introduction of a web-based Graphical User Interface (Web UI). The goal of this cycle is to elevate the tool from a functional utility for experts to a powerful and user-friendly framework that is accessible to a much broader audience of scientists and researchers.

On the algorithmic front, this cycle will introduce more advanced and intelligent methods for both the **Exploration** and **Sampling** stages. The Exploration Engine will be upgraded to support a hybrid Molecular Dynamics/Monte Carlo (MD/MC) simulation scheme. This is a critical enhancement for efficiently exploring the complex potential energy surfaces of materials like multi-component alloys, allowing the system to overcome high energy barriers through MC "swap moves." Furthermore, the engine will be made "smarter" by incorporating automatic ensemble switching, enabling it to correctly handle simulations of surfaces and slabs without user intervention. The Sampling module will be similarly upgraded from simple random selection to the more sophisticated Farthest Point Sampling (FPS) algorithm. FPS uses structural descriptors to select a maximally diverse subset of configurations, leading to more efficient and robust MLIP training.

The second major pillar of Cycle 2 is the development of a Web UI. While the CLI is powerful, it can be intimidating for new users and cumbersome for visualising results. The Web UI will provide an intuitive, interactive, and graphical way to use the MLIP-AutoPipe. Users will be able to configure their entire workflow through a web browser, launch calculations, monitor their progress in real-time, and, most importantly, visualise the atomic structures that are being generated. This feature will transform the user experience, making the tool not just a data generation pipeline but also an interactive platform for computational materials exploration. By the end of this cycle, MLIP-AutoPipe will be a significantly more powerful, intelligent, and accessible tool, fulfilling its core mission of accelerating materials discovery by automating and simplifying the creation of high-quality MLIP datasets.

## 2. System Architecture

The architecture for Cycle 2 expands upon the existing structure from Cycle 1, introducing new components for the advanced algorithms and an entirely new entry point for the Web UI. The core modular design is maintained, ensuring that the new features are well-integrated and the existing functionality remains robust.

The file structure below highlights the key files to be created or modified in this cycle. Files marked in **bold** are the primary targets for creation and implementation.

```text
src/mlip_autopipec/
├── __init__.py
├── cli.py
├── **web_ui.py**               # **New file:** Web UI application using FastAPI
│
├── core/
│   ├── __init__.py
│   ├── models.py             # **Modified:** Add new Pydantic models for UI state and FPS config
│   ├── pipeline.py
│   └── interfaces.py
│
├── components/
│   ├── __init__.py
│   ├── generators/
│   │   ├── ... (existing files)
│   ├── exploration/
│   │   ├── __init__.py
│   │   └── md_engine.py        # **Modified:** Add hybrid MD/MC and auto-ensemble logic
│   ├── sampling/
│   │   ├── __init__.py
│   │   ├── random.py
│   │   └── **fps.py**            # **New file:** Farthest Point Sampling implementation
│   └── storage/
│       ├── __init__.py
│       └── database.py
│
└── utils/
    ├── __init__.py
    ├── config_loader.py
    └── **physics.py**          # **New file:** For vacuum detection logic
```

**Architectural Blueprint:**

1.  **`web_ui.py`**: This is the most significant new file in Cycle 2. It will contain a web application built with a modern async framework like FastAPI. This choice allows for the creation of a responsive, non-blocking backend that can handle long-running simulation tasks. The application will serve an HTML frontend (likely a simple, single-page application) and provide a RESTful API for the frontend to interact with the pipeline. It will have API endpoints to:
    *   Start a new pipeline run with a configuration received as a JSON payload.
    *   Query the status of a running pipeline.
    *   Retrieve and display results and visualisations.

2.  **Pydantic Models (`core/models.py`)**: The existing configuration models will be extended.
    *   `ExplorationConfig`: A new field, `mc_moves: bool`, will be added.
    *   `SamplingConfig`: The `method` field will be updated to `method: Literal["random", "fps"]` and a new nested model, `FPSConfig`, will be added to hold parameters specific to Farthest Point Sampling (e.g., descriptor type).
    *   We will also add new models for managing the state of the Web UI, such as `PipelineStatus`, which can be easily serialised to JSON to communicate with the frontend.

3.  **Hybrid MD/MC Engine (`components/exploration/md_engine.py`)**: The existing `MDEngine` will be significantly refactored.
    *   The core simulation loop will be modified to conditionally insert Monte Carlo steps in between the MD steps.
    *   New methods will be added to perform the MC "swap move," which involves randomly selecting two atoms of different species and attempting to swap their positions.
    *   The logic for automatic ensemble switching will be implemented. Before starting a simulation, the engine will call a new utility function in `utils/physics.py` to detect the presence of a vacuum slab. Based on the result, it will programmatically select either the NVT or NPT ASE dynamics object, overriding the user's configuration if necessary to prevent simulation artefacts.

4.  **Vacuum Detection (`utils/physics.py`)**: A new utility file will be created to house the physics-based logic for vacuum detection. The `detect_vacuum` function will take an `ase.Atoms` object and implement a grid-based algorithm to determine if there is a significant region of empty space, indicative of a surface or slab model.

5.  **Farthest Point Sampling (`components/sampling/fps.py`)**: A new file will be added to implement the `FPSSampler`. This class will adhere to the `BaseSampler` interface. It will use a library like `dscribe` to compute SOAP (Smooth Overlap of Atomic Positions) descriptors for every atom in every structure from the exploration trajectory. It will then implement the iterative FPS algorithm to select a subset of structures that are maximally distant from each other in this high-dimensional descriptor space. This ensures a structurally diverse dataset.

## 3. Design Architecture

The design for Cycle 2 continues to embrace the schema-first philosophy with Pydantic, while also carefully considering the architectural challenges of introducing an interactive Web UI alongside the existing batch-processing CLI.

**Pydantic-Based Schema Design:**

The new features will be integrated into our existing Pydantic configuration schema, ensuring type safety and validation.

1.  **`ExplorationConfig` Enhancement**:
    *   A new nested model `MCConfig` will be added, containing fields like `swap_frequency: int` and `use_charge_safety: bool`. This allows for fine-grained control over the Monte Carlo steps.
    *   The main `ExplorationConfig` will have an optional `mc_config: Optional[MCConfig] = None` field. If this is present in the user's YAML file, the hybrid MD/MC mode is activated.

2.  **`SamplingConfig` Enhancement**:
    *   The `method` field will now be a `Literal["random", "fps"]`.
    *   A `fps_config: Optional[FPSConfig] = None` field will be added. The `FPSConfig` model will define parameters like `descriptor: Literal["soap"]`, `n_species: int`, `r_cut: float`, etc., required for computing the SOAP descriptors.
    *   A root validator will be added to `SamplingConfig` to enforce the constraint that if `method` is "fps", then `fps_config` must not be `None`. This prevents invalid configuration states.

3.  **Web UI State Models**: To manage the communication between the backend and frontend of the Web UI, we will define specific Pydantic models:
    *   `PipelineJob`: Represents a single pipeline run initiated by the UI. It will have fields like `job_id: UUID`, `status: Literal["running", "completed", "failed"]`, `progress: float`, and `config: FullConfig`.
    *   `StructureVisualization`: A model to send atomic structure data to the frontend for rendering. It will contain `elements: list[str]` and `positions: list[list[float]]` in a format that is easy for a JavaScript 3D viewer (like `nglview`) to consume.

**Key Invariants and Constraints:**

*   **Asynchronous Execution**: The Web UI must remain responsive while a long-running simulation is in progress. The FastAPI backend will be designed to run the MLIP-AutoPipe pipeline in a separate background process or thread pool (using `asyncio.to_thread` or similar), preventing the main web server from blocking.
*   **State Management**: The Web UI needs to manage the state of multiple concurrent jobs. A simple in-memory dictionary or a more robust solution like a Redis cache will be used to store the `PipelineJob` objects, keyed by their `job_id`.
*   **Component Purity**: The new algorithmic components (`FPSSampler`, the MC logic) will be designed as "pure" functions or classes as much as possible. They will take data as input and return new data as output, without causing side effects like writing to files directly. This makes them highly testable. The `WorkflowOrchestrator` remains the sole manager of I/O.

**Consumers and Producers:**

*   The **Web UI Frontend (JavaScript)** is a new consumer of the FastAPI backend. It consumes JSON data representing job statuses and structure visualisations. It is a producer of new job requests, sending JSON payloads containing the user's desired configuration.
*   The **`web_ui.py` (FastAPI backend)** consumes these JSON requests, converts them into our Pydantic `FullConfig` model, and produces background jobs. It consumes the status updates from the running jobs and produces JSON responses for the frontend.
*   The **`FPSSampler`** consumes the large trajectory files from the exploration stage and produces a small, curated list of `ase.Atoms` objects.
*   The **MC logic within `MDEngine`** consumes `ase.Atoms` objects and produces modified `ase.Atoms` objects with swapped atomic species.

## 4. Implementation Approach

The implementation of Cycle 2 is split into two parallel streams: the backend algorithmic enhancements and the frontend Web UI development.

**Stream 1: Algorithmic Enhancements**

**Step 1: Implementing Vacuum Detection (`utils/physics.py`)**
We will start with the simplest new, isolated component. We will create the `detect_vacuum` function. The implementation will involve creating a 3D grid over the simulation cell, and for each grid point, calculating the distance to the nearest atom. If there is a contiguous region of grid points that are all further than a certain cutoff from any atom, this indicates a vacuum slab. This function will be thoroughly unit tested with known bulk and slab structures.

**Step 2: Upgrading the Exploration Engine (`components/exploration/md_engine.py`)**
This is a significant refactoring task. We will modify the main simulation loop in `MDEngine`.
*   First, we will add the call to our new `detect_vacuum` function at the beginning of a run. We will implement the logic to programmatically select the NVT or NPT ensemble based on the result.
*   Next, we will add a conditional block inside the main loop that checks if `mc_config` is enabled. If it is, after a specified number of MD steps, the MC logic will be triggered.
*   We will implement a `_perform_swap_move` method. This method will randomly select two atoms of different species, check if the swap is allowed (e.g., charge safety), and if so, swap their atomic numbers in the `ase.Atoms` object.

**Step 3: Implementing Farthest Point Sampling (`components/sampling/fps.py`)**
This is a new, self-contained component.
*   We will create the `FPSSampler` class, inheriting from the `BaseSampler` interface from Cycle 1.
*   The main `sample()` method will first load all structures from the exploration trajectory.
*   It will then use a library like `dscribe` to set up and compute the SOAP descriptors for all structures. This will result in a large matrix where each row represents a structure and the columns are the descriptor values.
*   Finally, it will implement the FPS algorithm. This is an iterative process:
    1.  Select a random structure to start.
    2.  In a loop, find the structure that is "farthest away" (in Euclidean distance in the descriptor space) from the set of already selected structures.
    3.  Add this structure to the selected set.
    4.  Repeat until the desired number of samples is reached.

**Stream 2: Web UI Development**

**Step 4: Setting up the FastAPI Backend (`web_ui.py`)**
We will create the `web_ui.py` file and set up a basic FastAPI application. We will define the Pydantic models for API communication (`PipelineJob`, etc.). We will create the main API endpoints:
*   `POST /run`: This endpoint will accept a JSON payload representing the `FullConfig`, validate it, create a `PipelineJob` object with a unique ID, store it in our job cache, and start the pipeline in a background process. It will immediately return the `job_id`.
*   `GET /status/{job_id}`: This endpoint will look up the job in the cache and return its current status (`running`, `completed`, etc.) and progress percentage.

**Step 5: Frontend Development (HTML/JavaScript)**
We will create a simple `index.html` file to be served by FastAPI. This file will contain:
*   A web form with input fields for all the key configuration parameters (elements, temperature, etc.).
*   A "Run" button that, when clicked, uses JavaScript's `fetch` API to send the form data to the `POST /run` endpoint.
*   A status display area that, after a run is started, periodically polls the `GET /status/{job_id}` endpoint to update the user on the progress.
*   A placeholder for a 3D structure viewer. We will integrate a library like `nglview.js` or `3Dmol.js` to render the atomic structures returned from the backend.

## 5. Test Strategy

The test strategy for Cycle 2 must cover the new complex algorithms and the full stack of the new Web UI. We will expand our existing `pytest` suite and introduce end-to-end browser testing.

**Unit Testing Approach (Min 300 words):**
We will continue our practice of writing isolated unit tests for all new business logic.

1.  **Vacuum Detection (`utils/physics.py`):** The `detect_vacuum` function is highly suitable for unit testing. We will create several `ase.Atoms` objects in our test suite: a standard bulk crystal (e.g., face-centered cubic), a slab with a vacuum layer on top, and a molecule in a large box. We will then call `detect_vacuum` on each and assert that it returns `False` for the bulk crystal and `True` for the other two. This ensures our geometric analysis is correct.

2.  **Hybrid MD/MC Logic (`components/exploration/md_engine.py`):** We will write specific unit tests for the MC moves. We will create a simple `Atoms` object with two different elements. We will call the internal `_perform_swap_move` method and then inspect the object to assert that the atomic numbers of two atoms have been swapped, while everything else (positions, cell) remains unchanged. We will also test the charge-safety feature by creating a mock ionic system and asserting that an invalid swap is correctly rejected.

3.  **FPS Sampler (`components/sampling/fps.py`):** This is a critical component to test. We cannot rely on the exact output of a real SOAP descriptor library in a unit test. Instead, we will mock the descriptor calculation. The test will provide the `FPSSampler` with a pre-defined NumPy array of mock "descriptors" where the distances are known. For example, we can create three points forming an equilateral triangle and a fourth point far away. We will then run the `sample(n=2)` method and assert that the algorithm correctly selects the two points that are farthest apart. This validates the correctness of our FPS implementation, independent of the descriptor generation.

**Integration & End-to-End Testing Approach (Min 300 words):**
For Cycle 2, we need to expand from backend-only integration tests to full-stack end-to-end (E2E) tests for the Web UI.

1.  **Backend Integration Tests:** We will add new integration tests for the CLI that enable the new features. For example, we will have a test that runs the full pipeline with `method: fps` in the configuration. This test will run a small, complete workflow and check that the final database is created. The main purpose is to ensure the `FPSSampler` is correctly integrated with the `PipelineRunner` and that the required dependencies (like `dscribe`) are handled correctly.

2.  **Web UI E2E Tests:** This is the most important new type of testing in Cycle 2. We will use a browser automation framework like **Playwright**.
    *   **Setup:** The test suite will use `pytest` fixtures to start the FastAPI web server as a background process at the beginning of the test session.
    *   **Test Execution:** The test function will use the Playwright API to automate a real web browser (e.g., Chrome).
    *   **Scenario:** A typical E2E test will perform the following steps:
        a. `page.goto("http://127.0.0.1:8000")`: Navigate to the web application.
        b. `page.fill("#elements", "Si,Ge")`: Use CSS selectors to find the input fields and fill out the form for a new calculation.
        c. `page.click("#run-button")`: Simulate a user clicking the "Run" button.
        d. The test will then assert that the UI updates to show a "Running" status. It will use `page.wait_for_selector` to wait for the status element to appear.
        e. The test will poll the status endpoint or check the UI periodically until the status changes to "Completed".
        f. Finally, it will check that a result is displayed, for example, by asserting that an `nglview` canvas element is now visible on the page.
    *   This E2E test provides the highest possible level of confidence, as it verifies the entire system, from the user's click in the browser, through the FastAPI backend, down to the execution of the pipeline, and back to the UI.
