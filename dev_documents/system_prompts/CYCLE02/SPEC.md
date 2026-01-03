# Specification: Cycle 2 - Advanced Features & Web UI

## 1. Summary

This document provides the detailed technical specification for Cycle 2 of the MLIP-AutoPipe project. Building upon the foundational CLI tool developed in Cycle 1, this cycle will focus on implementing the advanced scientific features and the user-friendly web interface that are core to the project's value proposition. The primary objective is to significantly enhance the quality and diversity of the generated datasets and to make the tool accessible to a broader audience.

Key technical deliverables for this cycle include a sophisticated exploration engine with a hybrid Molecular Dynamics/Monte Carlo (MD/MC) sampling capability, an advanced Farthest Point Sampling (FPS) method for intelligent data selection, and robust parallel processing to accelerate the computationally intensive exploration stage. Furthermore, this cycle will introduce a web-based User Interface (UI) that will allow users to configure and launch pipeline runs interactively. The web UI will provide a more intuitive user experience compared to the command line, lowering the barrier to entry. By the end of this cycle, MLIP-AutoPipe will be a feature-complete and powerful tool, capable of generating highly diverse and physically realistic datasets through both a powerful CLI and an intuitive web interface.

## 2. System Architecture

Cycle 2 of the MLIP-AutoPipe project builds upon the solid foundation established in Cycle 1, introducing a suite of advanced features that significantly enhance its scientific capabilities and user accessibility. The core architectural principles of modularity and separation of concerns are maintained, while existing components are extended and new ones are introduced. This cycle focuses on three primary enhancements: the implementation of a sophisticated hybrid MD/MC exploration engine, the introduction of an intelligent Farthest Point Sampling (FPS) algorithm, and the development of a user-friendly web interface. Furthermore, a critical performance optimization will be introduced by enabling parallel execution of the exploration stage. The architectural changes are carefully designed to be additive, ensuring that the existing functionality remains robust and that the new features are integrated seamlessly.

**File Structure (Changes in Cycle 2):**

The following file structure highlights the key additions and modifications for this cycle. The changes are concentrated in evolving the existing pipeline components and adding a new top-level package for the web application.

```
src/mlip_autopipec/
├── **web/**                  # New package for the web UI
│   ├── **__init__.py**       # To be created
│   └── **app.py**            # Streamlit application entry point
├── config/
│   └── **models.py**         # Modified to include new configuration options
├── explorers/
│   └── **md_engine.py**      # Heavily modified for MD/MC and parallel support
├── samplers/
│   ├── **fps.py**            # New module for the FPS algorithm
│   └── **random_sampler.py** # Unchanged from Cycle 1
├── pipeline/
│   └── **orchestrator.py**   # Modified to support parallelism and sampler selection
└── ... (other files remain structurally the same)
```

**Component Blueprints and Modifications:**

This section details the specific changes and additions to the system's components, providing a blueprint for the implementation.

**`config/models.py` (Modification)**: The Pydantic configuration schema will be extended to provide users with fine-grained control over the new features. The `ExplorationConfig` model will be updated to include an optional `MonteCarloConfig` sub-model, which will contain parameters such as `enabled` and `swap_frequency`. A `num_parallel` field will also be added to allow users to specify the number of parallel simulations. The `SamplingConfig` model will see its `method` field changed from a simple string to a `typing.Literal['random', 'fps']`, which provides static type checking and auto-completion benefits, ensuring that only valid sampler names can be provided. These changes will be implemented in a backward-compatible manner, with default values ensuring that configurations from Cycle 1 still work as expected.

**`explorers/md_engine.py` (Modification)**: This module will undergo the most significant refactoring. The monolithic `run_md` function will be evolved into a more sophisticated `ExplorationEngine` class. The core simulation loop will be modified to support the hybrid MD/MC scheme. This involves interleaving standard MD steps (using the `Langevin` dynamics) with Monte Carlo steps. After a set number of MD steps (determined by `swap_frequency`), the simulation will be paused, and a Monte Carlo swap move will be attempted using `ase.mc.AtomSwap`. This move proposes swapping the positions of two atoms of different species, and the move is accepted or rejected based on the standard Metropolis criterion, which evaluates the change in the system's potential energy. This introduces a powerful mechanism for exploring chemical configuration space, which is often inaccessible to pure MD.

**`samplers/fps.py` (New File)**: A new `FPSSampler` class will be implemented in this file, adhering to the `BaseSampler` interface. This class will provide a more intelligent method for selecting structures than the random baseline. Its `sample` method will be a multi-step process. First, it will require a method to "fingerprint" each atomic structure in the trajectory. For this, it will use the SOAP (Smooth Overlap of Atomic Positions) descriptor, a powerful tool for representing local atomic environments. A dependency on the `dscribe` library will be added to compute these SOAP vectors. Second, with a fingerprint vector for each structure, the core FPS algorithm will be implemented. This is an iterative process: it begins by selecting a random structure for the final dataset. Then, in each subsequent step, it finds the structure in the remaining pool that is "farthest" from all the structures already selected (i.e., its fingerprint has the largest minimum distance to the fingerprints of the selected set). This process continues until the desired number of samples is reached, guaranteeing a structurally diverse dataset.

**`pipeline/orchestrator.py` (Modification)**: The `WorkflowOrchestrator` will be enhanced to manage the new complexity. Firstly, it will be updated to use a factory pattern for the samplers. Instead of being hardcoded to use `RandomSampler`, it will inspect the `sampling.method` field from the configuration and instantiate the appropriate class (`RandomSampler` or `FPSSampler`). Secondly, and more importantly, it will be refactored to handle the parallel execution of the exploration stage. The simple loop over seed structures will be replaced by a `concurrent.futures.ProcessPoolExecutor`. The orchestrator will submit an `ExplorationEngine.run` task to the process pool for each seed structure. This will allow multiple MD/MC simulations to run concurrently, dramatically reducing the wall-clock time required for this computationally intensive stage. The orchestrator will be responsible for managing the pool, submitting the tasks, and collecting the trajectory results from each process as they complete, before aggregating them for the sampling stage. This change requires careful handling of data serialization between processes, reinforcing the "late binding" pattern for the ASE calculator to avoid pickling large objects.

**`web/app.py` (New File)**: This new module will introduce a completely new way for users to interact with the application. It will contain a web application built using the Streamlit framework. This choice allows for rapid, Python-only development of an interactive UI. The application will present the user with a clean, form-based interface with widgets like sliders, text inputs, and dropdown menus that correspond directly to the parameters in the `FullConfig` Pydantic model. When the user clicks a "Submit" button, the application will gather the state of these widgets, construct a dictionary, and use it to instantiate the `FullConfig` model, thereby reusing the same robust validation logic as the CLI. This validated configuration will be written to a temporary YAML file. The application will then use Python's `subprocess` module to invoke the existing CLI entry point (`mlip-autopipec run-pipeline`), passing the path to the temporary config file. A key feature will be the capturing of the subprocess's standard output and error streams in real-time and displaying them in a text area on the web page, providing the user with immediate feedback on the pipeline's progress. Upon successful completion, a download button will be generated, allowing the user to retrieve the final database file directly from their browser. This provides a gentle learning curve for new users and a convenient alternative for experienced users for routine tasks.

## 3. Design Architecture

The schema-first design continues in Cycle 2, with Pydantic models being updated to support the new functionality.

**Pydantic Schema Design (`config/models.py` Modifications):**

The existing models will be updated. The `ExplorationConfig` and `SamplingConfig` will see the most significant changes.

```python
# config/models.py (updated sections)

from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Literal

# ... SystemConfig remains the same ...

class MonteCarloConfig(BaseModel):
    """Configuration for Monte Carlo moves."""
    enabled: bool = Field(False, description="Enable hybrid MD/MC.")
    swap_frequency: int = Field(50, gt=0, description="Perform atom swap move every N steps.")
    # Add other MC move configs here if needed

class ExplorationConfig(BaseModel):
    """Configuration for the exploration (MD) stage."""
    temperature_k: float = Field(300.0, gt=0, description="Simulation temperature in Kelvin.")
    num_steps: int = Field(1000, gt=0, description="Number of MD steps to perform.")
    potential: str = Field("EMT", description="ASE-compatible potential to use (e.g., EMT).")
    mc_config: Optional[MonteCarloConfig] = Field(None, description="Monte Carlo settings.")
    num_parallel: int = Field(1, ge=1, description="Number of simulations to run in parallel.")

class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""
    method: Literal['random', 'fps'] = Field("random", description="Sampling method ('random' or 'fps').")
    num_samples: int = Field(100, gt=0, description="Number of structures to sample from the trajectory.")
    # FPS-specific settings could be added here in a sub-model if needed

class FullConfig(BaseModel):
    """The root configuration model."""
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    db_path: str = Field("mlip_data.db", description="Path to the output ASE database file.")
```

**Web UI Design (`web/app.py`):**
-   The web application will be built using Streamlit to allow for rapid development.
-   The UI will present a form with input fields that directly map to the `FullConfig` Pydantic model. Widgets like sliders for temperature, number inputs for counts, and dropdowns for methods will be used.
-   Upon submission, the UI will construct a dictionary from the form data, instantiate the `FullConfig` model to validate it, and then serialize it to a temporary YAML file.
-   It will then use a `subprocess` call to execute the `mlip-autopipec run-pipeline` command, passing the path to the temporary YAML file.
-   The UI will display the real-time stdout/stderr from the subprocess, providing progress feedback to the user.
-   After the run completes, it will offer a link to download the resulting database file.

## 4. Implementation Approach

The implementation will be staged to tackle the backend scientific enhancements first, followed by the user-facing web UI.

1.  **Update Configuration:** Modify the Pydantic models in `config/models.py` as detailed above. Ensure backward compatibility by making new settings optional where possible (e.g., `mc_config`).
2.  **Implement FPS Sampler:**
    -   In `src/mlip_autopipec/samplers/fps.py`, create the `FPSSampler` class, inheriting from `BaseSampler`.
    -   The implementation will require a dependency on a library that can compute SOAP descriptors (e.g., `dscribe`). This new dependency must be added to `pyproject.toml`.
    -   The `sample` method will first compute SOAP vectors for all atoms in all frames of the trajectory. It will then implement the iterative FPS algorithm to select the `num_samples` structures whose average SOAP vectors are furthest from each other.
3.  **Enhance Exploration Engine:**
    -   Refactor `explorers/md_engine.py`. The `run_md` function will be modified to become `run_exploration`.
    -   It will check if `config.exploration.mc_config.enabled` is true.
    -   If so, the main simulation loop will be modified to periodically pause the MD run and attempt a Monte Carlo atom swap move using `ase.mc.AtomSwap`. The move is accepted or rejected based on the Metropolis criterion.
4.  **Implement Parallel Processing:**
    -   In `pipeline/orchestrator.py`, the `WorkflowOrchestrator.run` method will be modified.
    -   It will use Python's `concurrent.futures.ProcessPoolExecutor` to parallelize the loop that runs the exploration on the seed structures. The number of workers for the pool will be set by `config.exploration.num_parallel`.
    -   This requires careful implementation to ensure that data (especially the large trajectory lists) is passed back from the worker processes to the main process efficiently. The "late binding" of the calculator inside the worker process (already a good practice) becomes critical here.
5.  **Update Orchestrator Logic:**
    -   The `WorkflowOrchestrator` will be updated to read the `config.sampling.method`.
    -   It will use a factory pattern to instantiate the correct sampler class (`RandomSampler` or `FPSSampler`) based on this value.
6.  **Develop Web UI:**
    -   Create the `web/app.py` file using Streamlit.
    -   Build the input form as described in the Design Architecture section.
    -   Implement the logic to generate the config file, run the CLI as a subprocess, and display the output.

## 5. Test Strategy

Testing in Cycle 2 will cover the new functionality and ensure that the existing functionality from Cycle 1 has not regressed.

**Unit Testing Approach:**
-   **`config/models.py`:** Add new tests to validate the `mc_config` and the new `method` literal in `SamplingConfig`. Test that invalid values (e.g., `method='invalid_method'`) raise a `ValidationError`.
-   **`samplers/fps.py`:** Unit testing the `FPSSampler` is critical. Create a small, deterministic "trajectory" of a few `ase.Atoms` objects. Pre-compute the expected SOAP vectors (or use a mock library). Run the FPS algorithm on this small dataset and assert that it selects the exact subset of structures that are known to be the most diverse. This will likely involve testing the core FPS logic on simple NumPy arrays first.
-   **`explorers/md_engine.py`:** Testing the hybrid MD/MC logic is challenging. The test will focus on the control flow. The ASE `AtomSwap` and `Langevin` dynamics objects will be mocked. The test will run a small number of steps and assert that the `swap.run()` method was called the correct number of times, corresponding to the `swap_frequency`.
-   **`web/app.py`:** The backend logic of the Streamlit app can be unit tested. Create a function that takes a dictionary of form inputs and generates the configuration. Test this function to ensure it produces the correct YAML output. Mock the `subprocess.run` call and assert that it is called with the correct command-line arguments.

**Integration Testing Approach:**
The integration test from Cycle 1 (`tests/integration/test_cli.py`) will be extended with new test cases.
-   **FPS Pipeline Test:** A new test will be created that uses a configuration file with `sampling.method: 'fps'`. The `explorers.md_engine.run_exploration` will be mocked as before. The `samplers.fps.FPSSampler.sample` method will *also* be mocked to avoid the expensive descriptor calculation. The test will assert that the `FPSSampler` was chosen by the orchestrator and that its `sample` method was called. This confirms the new logic path in the orchestrator is working.
-   **Parallel Execution Test:** A test will be created that sets `exploration.num_parallel: 2`. The test will mock `concurrent.futures.ProcessPoolExecutor` and assert that it was initialized with `max_workers=2`. It will also assert that the `executor.submit` method was called the correct number of times (equal to `num_structures`).

**End-to-End (E2E) Testing:**
-   A new E2E test will be created for the web UI using Playwright.
-   The test script will:
    1.  Start the Streamlit application as a separate process.
    2.  Use Playwright to navigate to the local web server.
    3.  Interact with the UI widgets to fill out the form (e.g., `page.fill`, `page.click`).
    4.  Click the "Run Pipeline" button.
    5.  The test will mock the `subprocess.run` call at the system level so that the full pipeline doesn't actually run. The mock will simulate a successful run.
    6.  The test will then assert that the UI displays the expected "Success" message.
This test validates that the user interface is correctly wired to the backend CLI and can trigger a pipeline run.
