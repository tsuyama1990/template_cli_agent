# CYCLE02/SPEC.md

## 1. Summary

Cycle 2 represents the most significant and scientifically innovative phase of the MLIP-AutoPipe project, building directly upon the robust foundational architecture established in Cycle 1. The central goal of this cycle is to implement the sophisticated, high-value features that truly differentiate this tool: the hybrid Molecular Dynamics/Monte Carlo (MD/MC) exploration engine, the intelligent Farthest Point Sampling (FPS) module, and an intuitive, browser-based User Interface (UI). The exploration engine lies at the very heart of the project's value proposition. It is designed to intelligently and efficiently traverse the potential energy surface of a given material, moving beyond simple, static structures to find a rich diversity of atomic configurations that are representative of the material's behavior under real-world conditions. The FPS module complements this by providing an intelligent, data-driven method for selecting the most valuable and non-redundant data points from the vast output of the exploration engine, ensuring the final dataset is both compact and information-rich. Finally, the web UI will democratize access to this powerful computational tool, allowing users who are not command-line experts to configure, run, and visualize the entire workflow through a simple and intuitive graphical interface.

By the end of this cycle, the MLIP-AutoPipe project will be transformed from a functional but basic command-line tool into a comprehensive, powerful, and user-friendly framework for advanced MLIP dataset generation. A user will be able to not only generate valid initial structures but also subject them to simulated thermodynamic conditions, exploring a wide and physically relevant range of configurations. They will then be able to automatically curate a high-quality, diverse dataset from these complex simulations with the click of a button. The addition of the web UI is a crucial component of this vision, as it significantly lowers the barrier to entry. This makes the advanced simulation techniques, which were previously the domain of computational experts, accessible to a much broader audience of materials scientists, chemists, and experimental researchers. The successful completion of this cycle will therefore fulfill the core vision of the project: to create a fully automated, physically realistic, and powerful tool that accelerates research and development in the materials science community. It will provide a platform for discovery, enabling the creation of the next generation of accurate and reliable Machine Learning Interatomic Potentials.

## 2. System Architecture

In Cycle 2, we will build out the remaining, more complex modules of the application. The architectural focus will be on implementing the `explorers`, `samplers`, and `web` components, and integrating them seamlessly into the existing pipeline established in Cycle 1. The file structure is designed to accommodate this new functionality in a clean and modular fashion.

**File Structure for Cycle 2:**

The ASCII tree below details the new files to be created and the existing files to be modified. The files marked in **bold** are the new implementation targets for this cycle. The `pipeline.py`, `config.py`, and `cli.py` files, which were created in Cycle 1, will be **modified** to incorporate the new functionality. This clear structure ensures that the new, complex logic is well-organized and decoupled from the core framework.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py                     # **MODIFIED**: Add `run-webui` command.
├── config.py                  # **MODIFIED**: Add `ExplorationConfig` and `SamplingConfig` models.
├── database.py
├── pipeline.py                # **MODIFIED**: `PipelineRunner` to include exploration and sampling stages.
├── generators/
│   ├── __init__.py
│   ├── base.py
│   ├── alloy.py
│   └── ionic.py
├── **explorers/**                 # **New Package**: Contains all simulation and exploration logic.
│   ├── **__init__.py**
│   └── **md_engine.py**       # **Implementation Target**: The `HybridMDEngine` class.
├── **samplers/**                  # **New Package**: Contains all sampling algorithms.
│   ├── **__init__.py**
│   ├── **base.py**            # **Implementation Target**: Defines the `BaseSampler` Abstract Base Class.
│   └── **fps.py**             # **Implementation Target**: The `FPSSampler` class.
└── **web/**                       # **New Package**: Contains the web user interface.
    └── **app.py**             # **Implementation Target**: FastAPI application for the web UI.
tests/
├── conftest.py
├── test_cli.py
├── test_config.py
├── test_generators.py
├── test_database.py
├── **test_explorers.py**          # **Implementation Target**: Unit tests for the `HybridMDEngine`.
└── **test_samplers.py**           # **Implementation Target**: Unit tests for the `FPSSampler`.
```
This architecture ensures a continued separation of concerns. The `HybridMDEngine` in `explorers/md_engine.py` will contain all the complex physics simulation code, depending on the ASE library but remaining independent of the data storage or pipeline orchestration logic. The `FPSSampler` will similarly encapsulate the sampling algorithm, taking a list of `ase.Atoms` objects and returning a subset, without needing to know where the atoms came from or where they are going. The modifications to `pipeline.py` are purely orchestrational—it will be updated to conditionally call these new components based on the configuration. This "pluggable" architecture, where new stages can be added to the pipeline without fundamentally changing the existing stages, is the key to the system's long-term maintainability. The `web/app.py` module will be a self-contained web server, communicating with the core application by calling the `PipelineRunner` in a separate process, ensuring that the web interface is a distinct layer that does not become entangled with the core scientific logic.

## 3. Design Architecture

The design in Cycle 2 continues to follow the rigorous schema-first, modular, and contract-driven approach established in Cycle 1. All new functionality will be built upon these principles.

**Pydantic-Based Schema Design (`config.py`):**

The `FullConfig` model will be extended to include detailed and validated configuration models for the new exploration and sampling components. This ensures that the powerful new features are exposed to the user in a clear, predictable, and safe manner.

*   `ExplorationConfig`: This model will contain all settings for the MD/MC engine.
    *   `run_md`: A `bool` to enable or disable this stage entirely.
    *   `temperature`: A `float` with `Field(gt=0)` validation to specify the simulation temperature in Kelvin.
    *   `pressure`: An optional `float` with `Field(ge=0)` for the NPT ensemble pressure in GPa. If `None`, an NVT ensemble will be used.
    *   `n_steps`: A `PositiveInt` for the total number of MD steps.
    *   `time_step`: A `PositiveFloat` for the MD time step in femtoseconds.
    *   `mc_swap_frequency`: An optional `PositiveInt` to control how often Monte Carlo atom swap moves are attempted (e.g., every 100 MD steps).

*   `SamplingConfig`: This model will define the sampling strategy.
    *   `method`: An `Enum` of available methods, e.g., `'random'`, `'fps'`.
    *   `n_samples`: A `PositiveInt` specifying the desired size of the final dataset.
    *   A custom `@model_validator` will be added to ensure that `n_samples` is not greater than the total number of structures that will be generated by the exploration phase.

*   The `FullConfig` model will be updated to include `exploration: ExplorationConfig` and `sampling: SamplingConfig`. Default values will be provided for all new fields to ensure backward compatibility with configurations from Cycle 1.

**Component Interfaces (Contracts):**

*   **`explorers.md_engine.HybridMDEngine`:**
    This class will be the public interface to the complex simulation logic.
    ```python
    from ase import Atoms
    from ..config import ExplorationConfig, SystemConfig

    class HybridMDEngine:
        def __init__(self, system_config: SystemConfig, exploration_config: ExplorationConfig): ...

        def run(self, initial_structures: list[Atoms]) -> list[Atoms]:
            """
            Runs the MD or hybrid MD/MC simulation for each initial structure.
            Returns a single, large list of Atoms objects representing the combined trajectory.
            """
            pass
    ```
    The engine will be initialized with the relevant configuration sub-models. Its `run` method accepts the list of seed structures from the Generation stage and returns the raw trajectory.

*   **`samplers.base.BaseSampler` (ABC):**
    This new Abstract Base Class will define the contract for all sampling algorithms, ensuring they are interchangeable.
    ```python
    from abc import ABC, abstractmethod
    from ase import Atoms

    class BaseSampler(ABC):
        @abstractmethod
        def sample(self, trajectory: list[Atoms], n_samples: int) -> list[Atoms]:
            """Selects a subset of n_samples Atoms objects from the trajectory."""
            pass
    ```
    This simple yet powerful abstraction allows the `PipelineRunner` to be completely agnostic to the details of the sampling method. It simply knows it has an object with a `sample` method that will return the desired number of structures. This makes adding new sampling algorithms in the future a trivial task from an architectural perspective.

## 4. Implementation Approach

The implementation of Cycle 2 will be carefully staged to manage its complexity. The backend simulation and sampling logic will be developed and tested first, followed by its integration into the pipeline, and finally, the development of the web UI as a separate layer on top.

1.  **Extend Configuration (`config.py`):** The first step is to update the schema. The `ExplorationConfig` and `SamplingConfig` Pydantic models will be implemented as designed. The main `FullConfig` model will be updated to include them. The corresponding unit tests in `tests/test_config.py` will be expanded to cover these new models, including tests for their validation logic.

2.  **Implement MD/MC Engine (`explorers/md_engine.py`):** This is the most substantial and complex task of the cycle.
    a. The core MD loop will be implemented using the powerful integrators available in the `ase.md` module (e.g., `Langevin` for NVT, `NPT` for constant pressure).
    b. The logic for attaching an MLIP calculator (e.g., MACE) to the `ase.Atoms` objects will be implemented.
    c. A key feature, the automatic ensemble switching, will be implemented. This will involve writing a helper function that can detect the presence of a significant vacuum layer in an `ase.Atoms` object, and using this to decide whether to use an NVT or NPT integrator for that specific simulation.
    d. The hybrid MC functionality will be implemented. An `ase.md.MDLogger`-like observer will be created that, every `mc_swap_frequency` steps, pauses the MD integration, attempts a random swap of two atoms of different species, and accepts or rejects the move based on the Metropolis criterion.
    e. Unit tests for the engine will be written in `tests/test_explorers.py`, using mock calculators to ensure the tests are fast and deterministic.

3.  **Implement Sampler Interface and FPS (`samplers/`):**
    a. The `BaseSampler` ABC will be created in `base.py`.
    b. The `FPSSampler` class will be implemented in `fps.py`. This will require adding a dependency on a library capable of calculating SOAP descriptors, such as `dscribe`. The `sample` method will first iterate through the entire input trajectory to compute and store the SOAP vector for each structure. It will then implement the iterative FPS algorithm to select the `n_samples` structures that are farthest from each other in this high-dimensional SOAP space.
    c. Unit tests will be added to `tests/test_samplers.py` using a small, pre-computed set of structures and descriptors to verify the correctness of the FPS implementation.

4.  **Update Pipeline (`pipeline.py`):** With the new backend components complete, the `PipelineRunner` will be updated to orchestrate the full workflow.
    a. Its `run` method will be modified to check the configuration and, if `exploration.run_md` is true, it will instantiate and run the `HybridMDEngine` after the generation stage.
    b. The output of the exploration stage (the raw trajectory) will then be passed to the sampler stage. The runner will use a factory to instantiate the correct sampler based on `sampling.method` and call its `sample` method.
    c. The final curated list of atoms from the sampler will be passed to the `AseDBWrapper` for storage. The integration tests for the CLI will be updated to test this full, four-stage pipeline.

5.  **Implement Web UI (`web/app.py` and `cli.py`):**
    a. A new command, `run-webui`, will be added to `cli.py`. This command will use `uvicorn` to programmatically start the FastAPI web server defined in `web/app.py`.
    b. The FastAPI application in `web/app.py` will be created. It will have at least two API endpoints:
        i. A `GET "/"` endpoint that serves a single static `index.html` file.
        ii. A `POST "/run"` endpoint that accepts a JSON payload conforming to the `FullConfig` Pydantic model. This endpoint will launch the `PipelineRunner.run` method in a background thread or separate process to avoid blocking the server. It will return a job ID to the client.
    c. A simple frontend (HTML, CSS, and vanilla JavaScript) will be created. It will provide a web form that allows the user to specify all the configuration parameters. A "Start Run" button will serialize the form data to JSON and `POST` it to the `/run` endpoint.

## 5. Test Strategy

The testing strategy for Cycle 2 must expand to cover the increased complexity, particularly the stochastic nature of MD simulations and the asynchronous, interactive nature of a web UI.

**Unit Testing Approach (Min 300 words):**
*   **MD/MC Engine (`test_explorers.py`):** Unit testing the `HybridMDEngine` requires a heavy reliance on mocking to ensure deterministic and fast tests. The MLIP calculator is the primary dependency to be mocked. We will use `pytest-mock`'s `mocker` fixture to replace the calculator object with a `MagicMock`. This mock can be configured to return predictable, deterministic values for energy and forces for any given `ase.Atoms` object. This allows us to test the engine's logic without performing any real, time-consuming calculations. We can write a test for the automatic ensemble switching by creating two `ase.Atoms` objects—one representing a bulk crystal and one with a large vacuum layer. We will then pass each to the engine and assert that the correct ASE dynamics object (`NPT` for bulk, `NVT` for the slab) was instantiated internally. For the MC swap logic, we can run a mocked 1-step simulation and assert that the positions of two different atom types have been swapped correctly.
*   **Samplers (`test_samplers.py`):** The `FPSSampler` must be tested for correctness. A robust unit test involves creating a fixture that provides a small, fixed list of `ase.Atoms` objects (e.g., 5-10 structures) where the structural diversity is known by construction. For example, we could have a set of structures with varying lattice constants. The test will also need a mocked version of the SOAP descriptor calculator that returns pre-computed, deterministic vectors for these specific structures. The test will then instantiate the `FPSSampler`, call `sample`, and assert that it returns the exact subset of structures that are known to be the most diverse. This validates the core logic of the FPS algorithm implementation.

**Integration Testing Approach (Min 300 words):**
*   **Full Pipeline via CLI:** A new integration test will be added to `test_cli.py` to validate the full, four-stage pipeline. To make this test runnable in a CI environment, it cannot use a real, slow MLIP. Instead, the test will be configured to use a very fast and simple classical potential provided by ASE, such as the Lennard-Jones potential for Argon. The test will configure a very small system (e.g., a 10-atom Argon cluster) and run a very short MD simulation (e.g., only 20 steps). After the `CliRunner` invokes the pipeline, the test will perform a series of assertions on the final database. It will assert that the database contains the correct *sampled* number of structures. Crucially, it will also assert that the average potential energy of the final structures is different from the energy of the initial structure, providing quantitative proof that the exploration stage actually ran and modified the system.
*   **Web UI (End-to-End):** Testing the web UI requires a dedicated end-to-end (E2E) testing strategy. We will use a framework like Playwright, which allows us to write Python scripts that launch and control a real web browser. A new test suite will be created in `tests/e2e/`. A typical test function would:
    1.  Use Python's `subprocess` module to start the FastAPI web server in the background.
    2.  Use Playwright's API to launch a browser and navigate to the web server's address.
    3.  Programmatically interact with the page: fill in the input fields of the configuration form with values for a simple test case.
    4.  Simulate a user clicking the "Start Run" button.
    5.  Assert that the UI updates to show a "Running" status message.
    6.  The test will then poll the filesystem, waiting for the output database file to be created by the backend process.
    7.  Once the file appears, the test can assert that the UI status has updated to "Completed".
    This E2E test provides the highest level of confidence, verifying that all layers of the system—from the user's browser interaction to the backend API, the pipeline orchestration, and the final file output—are working together correctly.