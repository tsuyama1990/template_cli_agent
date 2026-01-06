# SPEC.md: Cycle 1 - Core Framework and Basic Pipeline

## 1. Summary

This document provides the detailed technical specifications for Cycle 1 of the MLIP-AutoPipe project. The primary objective of this cycle is to deliver a Minimum Viable Product (MVP) that establishes the core architectural foundation and implements a complete, end-to-end data generation workflow. This foundational work will enable the creation of a basic, yet functional, dataset for a simple binary alloy system. The scope of this cycle is intentionally focused on creating a stable, testable, and robust framework upon which the more advanced features of Cycle 2 can be built.

Functionally, the system at the end of Cycle 1 will be a Command-Line Interface (CLI) tool. This tool will accept a single YAML configuration file as input, which specifies the parameters for a data generation run. The tool will then execute a four-stage pipeline: (1) Generation of initial atomic structures for a specified alloy, (2) Exploration of new structures using a simple Molecular Dynamics (MD) simulation with a computationally inexpensive potential, (3) Sampling of a subset of the generated structures using a straightforward random selection method, and (4) Storage of the final, curated structures into a standard ASE (Atomic Simulation Environment) database.

Technically, this cycle involves several key tasks. First is the scaffolding of the entire project structure, including the creation of all necessary directories and `__init__.py` files to define the Python package `mlip_autopipec`. A critical task is the implementation of a comprehensive suite of Pydantic models to define and validate all configuration parameters. This schema-first approach is paramount for ensuring robustness and preventing configuration errors. The core pipeline logic will be encapsulated within a `PipelineOrchestrator` class, which will manage the sequential execution of the four stages. For this cycle, concrete implementations for each stage will be developed: an `AlloyGenerator` for structure creation, a basic `ExplorationEngine` that uses ASE's built-in EMT potential for the MD simulation, a `RandomSampler`, and an `AseDBWrapper` for database interactions. The entire system will be driven by a `Typer`-based CLI. Crucially, this cycle also includes the development of a comprehensive test suite, covering unit tests for each module and a full end-to-end integration test to validate the entire workflow.

## 2. System Architecture

The architecture for Cycle 1 is focused on implementing the foundational components of the system. The file structure is designed to be modular and extensible, clearly separating concerns for each part of the pipeline.

**File Structure (Cycle 1 Focus):**

The files and directories to be created or modified in this cycle are marked in **bold**.

```
.
├── **pyproject.toml**
└── src/
    └── **mlip_autopipec/**
        ├── **__init__.py**
        ├── **cli/**
        │   ├── **__init__.py**
        │   └── **main.py**              # CLI application logic
        ├── **core/**
        │   ├── **__init__.py**
        │   ├── **orchestrator.py**      # PipelineOrchestrator class
        │   └── **models.py**            # Pydantic configuration models
        ├── **generators/**
        │   ├── **__init__.py**
        │   ├── **base.py**              # BaseStructureGenerator ABC
        │   └── **alloy.py**             # AlloyGenerator implementation
        ├── **exploration/**
        │   ├── **__init__.py**
        │   └── **engine.py**            # Basic ExplorationEngine
        ├── **sampling/**
        │   ├── **__init__.py**
        │   ├── **base.py**              # BaseSampler ABC
        │   └── **random_sampler.py**    # RandomSampler implementation
        └── **storage/**
            ├── **__init__.py**
            └── **database.py**          # AseDBWrapper class
```

**Code Blueprints:**

*   **`pyproject.toml`**: This file will be created to manage project dependencies and tooling. It must include `pydantic`, `typer`, `pyyaml`, `ase`, `pytest`, `pytest-cov`, and `pytest-mock` as dependencies. It will also be configured with the strict `ruff` and `mypy` settings as defined in the architectural overview.

*   **`cli/main.py`**:
    *   Will contain a `Typer` application.
    *   A main command, `run`, will be defined using `@app.command()`.
    *   This command will accept one argument: `config_path: Path`, which is the path to the YAML configuration file.
    *   Inside this function, it will first load the YAML file and parse it into the `FullConfig` Pydantic model. Error handling for file-not-found or parsing errors is required.
    *   It will then instantiate the `PipelineOrchestrator` with the loaded configuration.
    *   It will call the `orchestrator.run()` method to execute the pipeline.
    *   Rich progress indicators (e.g., `rich.progress`) should be used to provide user feedback during the run.

*   **`core/orchestrator.py`**:
    *   The `PipelineOrchestrator` class will be the central component.
    *   Its `__init__` method will accept a `FullConfig` object.
    *   It will have a public method `run()`.
    *   The `run()` method will instantiate and call the components for each of the four stages in sequence:
        1.  Instantiate `AlloyGenerator`, call `generate()`.
        2.  Instantiate `ExplorationEngine`, call `run_simulations()` with the generated structures.
        3.  Instantiate `RandomSampler`, call `sample()` with the simulation trajectories.
        4.  Instantiate `AseDBWrapper`, call `write_results()` with the final sampled structures.
    *   It will manage the data flow between these stages (e.g., passing the list of `ase.Atoms` objects).

*   **`generators/base.py` & `generators/alloy.py`**:
    *   `base.py` will define an abstract base class `BaseStructureGenerator` using `abc.ABC`. It will have one abstract method `generate(self) -> list[Atoms]`.
    *   `alloy.py` will contain `AlloyGenerator`, which inherits from `BaseStructureGenerator`. It will implement the `generate` method to create random alloy supercells using ASE's `ase.build.bulk` and `ase.build.make_supercell` functions, followed by random replacement of atoms to achieve the target composition. It will perform a simple check for minimum atomic distances to ensure physical plausibility.

*   **`exploration/engine.py`**:
    *   The `ExplorationEngine` class will be implemented.
    *   It will have a method `run_simulations(self, initial_structures: list[Atoms]) -> list[Path]`.
    *   This method will use Python's `concurrent.futures.ProcessPoolExecutor` to run multiple MD simulations in parallel.
    *   Each worker process will take one `ase.Atoms` object, attach an `ase.calculators.emt.EMT` calculator, and run a short NVT simulation using `ase.md.velocityscaling.VelocityVerlet`.
    *   The trajectory of each simulation will be saved to a file (e.g., `trajectory_0.traj`, `trajectory_1.traj`). The method will return a list of paths to these trajectory files.

*   **`sampling/base.py` & `sampling/random_sampler.py`**:
    *   `base.py` will define the `BaseSampler` ABC with an abstract method `sample(self, trajectories: list[Path]) -> list[Atoms]`.
    *   `random_sampler.py` will implement `RandomSampler`, which reads all frames from the trajectory files and randomly selects a pre-defined number of frames to be included in the final dataset.

*   **`storage/database.py`**:
    *   `AseDBWrapper` will be implemented.
    *   It will have a `__init__` method that takes the database path.
    *   It will have a method `write_results(self, structures: list[Atoms])`. This method will connect to the SQLite database using `ase.db.connect`, iterate through the list of `Atoms` objects, and write each one to the database, ensuring that properties like energy and forces are included.

## 3. Design Architecture

The design for Cycle 1 is centered around establishing a robust, schema-driven foundation using Pydantic. This ensures that the core logic is built upon a bedrock of validated, type-safe data structures, which is critical for the reliability of a scientific computing pipeline.

**Pydantic Schema Design (`core/models.py`):**

The Pydantic models are the canonical representation of the system's configuration. They are not just data containers; they are the primary mechanism for input validation and for defining the "contract" between different parts of the system.

*   **`FullConfig(BaseModel)`**: This is the top-level container that aggregates all other configuration models. It will act as the single source of truth for a pipeline run.
    *   `model_config = ConfigDict(extra="forbid")` will be set to prevent users from passing unknown configuration keys, catching typos and misconfigurations early.

*   **`SystemConfig(BaseModel)`**: Defines the intrinsic properties of the material being studied.
    *   `elements: list[str]`: A list of chemical symbols. A `@field_validator` will be used to ensure they are valid element symbols.
    *   `composition: dict[str, float]`: A mapping from element to its fractional composition. A `@model_validator` will ensure the values sum to 1.0.

*   **`GenerationConfig(BaseModel)`**: Parameters controlling the initial structure creation.
    *   `num_structures: int = Field(..., gt=0)`: The number of initial seed structures to generate, must be positive.
    *   `crystal_structure: Literal["fcc", "bcc", "hcp"]`: The lattice type.

*   **`ExplorationConfig(BaseModel)`**: Parameters for the MD simulation.
    *   `temperature_k: float = Field(..., gt=0)`: Simulation temperature in Kelvin.
    *   `md_steps: int = Field(..., ge=0)`: Number of MD steps.
    *   `calculator: Literal["emt"]`: In Cycle 1, only the EMT calculator is supported. This strict literal type prevents users from attempting to use unsupported calculators.

*   **`SamplingConfig(BaseModel)`**: Parameters for the sampling stage.
    *   `method: Literal["random"]`: Only "random" is supported in Cycle 1.
    *   `num_samples: int = Field(..., gt=0)`: The final number of structures to select.

*   **Consumers and Producers:**
    *   **Producer:** The user creates a `config.yml` file. The `cli.main` module is the **consumer** of this YAML file and the **producer** of the `FullConfig` Pydantic object.
    *   **Consumers:** The `PipelineOrchestrator` is the primary consumer of the `FullConfig` object. It then dispatches the relevant sub-models (`GenerationConfig`, `ExplorationConfig`, etc.) to the respective pipeline stage components (Generator, Explorer), which are the ultimate consumers of the validated configuration data.

*   **Versioning and Extensibility:** The design explicitly prepares for future extension. By using `Literal` types for parameters like `calculator` and `method`, we create a clear point for extension in Cycle 2. For example, `calculator: Literal["emt", "mace"]` will be a simple, backward-compatible change. The nested structure of `FullConfig` allows for new configuration sections (e.g., `WebServerConfig`) to be added without breaking existing YAML files.

## 4. Implementation Approach

The implementation will be approached in a logical, step-by-step manner, starting from the data structures and moving outwards to the application logic. This ensures that each layer is built on a tested and validated foundation.

1.  **Project Setup:**
    *   Create the directory structure as outlined in the System Architecture section.
    *   Create `__init__.py` files in all subdirectories to make them importable Python packages.
    *   Create the `pyproject.toml` file and add all the necessary dependencies. Use `uv sync` to create the virtual environment and install them.

2.  **Pydantic Models (`core/models.py`):**
    *   Implement all Pydantic models as specified in the Design Architecture.
    *   Write unit tests for these models first. Test validation logic, ensuring that both valid and invalid data are handled correctly. This is a Test-Driven Development (TDD) approach for the core data structures.

3.  **Component Implementation (Bottom-up):**
    *   Implement the `storage.AseDBWrapper` and its unit tests. Mock the `ase.db.connect` call to test the logic without file I/O.
    *   Implement the `generators.base.BaseStructureGenerator` and `generators.alloy.AlloyGenerator`. Write unit tests to verify that the generated structures have the correct properties.
    *   Implement the `sampling.base.BaseSampler` and `sampling.random_sampler.RandomSampler`. Unit test its logic with dummy trajectory files.
    *   Implement the `exploration.engine.ExplorationEngine`. This is more complex to test. The unit test should mock `ProcessPoolExecutor` and the ASE MD functions. The test should verify that the engine correctly prepares and submits the right number of jobs to the pool, not that the MD itself is correct.

4.  **Orchestrator and CLI (Top-down):**
    *   Implement the `core.orchestrator.PipelineOrchestrator`. Use `pytest-mock` to inject mock objects for the generator, explorer, sampler, and storage components. Write unit tests to verify that the `run()` method calls each component's methods in the correct sequence.
    *   Implement the `cli.main.py` Typer application. Use `typer.testing.CliRunner` to write tests. These will be integration-style tests that check the application's behavior from the command line.
        *   Test with a non-existent config file.
        *   Test with a malformed YAML file.
        *   Test with a YAML file that violates Pydantic validation rules.
        *   Test a successful run by mocking the `PipelineOrchestrator` to avoid running the full, slow pipeline.

5.  **Full Integration Test:**
    *   Create a final integration test that runs the actual CLI command on a minimal, valid `config.yml`.
    *   This test will not use mocks for the pipeline components. It will execute the entire, real workflow.
    *   It will run a very short simulation (e.g., 2 atoms, 5 MD steps).
    *   After the CLI command finishes, the test will use `ase.db.connect` to open the generated database file and assert that it contains the expected number of structures. This validates that all components work together correctly.

## 5. Test Strategy

The test strategy for Cycle 1 is focused on building a high-coverage test suite that validates the core functionality and ensures the framework is robust and reliable.

**Unit Testing Approach (Min 300 words):**

The unit testing approach for Cycle 1 is strictly component-focused, aiming to verify the correctness of each class and function in isolation. This is crucial for building a reliable system, as it ensures that the fundamental building blocks are sound before they are assembled. Mocks are used extensively to isolate the "unit under test" from its dependencies, leading to fast, deterministic, and targeted tests.

For the Pydantic models in `core/models.py`, the tests will not just check successful instantiation but will rigorously probe the validation logic. Using `pytest.raises(ValidationError)`, we will assert that the models correctly reject invalid data—for example, a `SystemConfig` where the composition does not sum to 1.0, or an `ExplorationConfig` with a negative temperature. This ensures the data integrity contract is enforced at the system's entry point.

The `generators.alloy.AlloyGenerator` will be tested by calling its `generate()` method and inspecting the returned list of `ase.Atoms` objects. Assertions will verify the number of atoms, the ratio of chemical species, and that the total number of structures matches the configuration. We will also mock any external I/O or computationally intensive parts to keep the tests fast.

Testing the `exploration.engine.ExplorationEngine` requires a more sophisticated use of mocking. We will not test the correctness of the ASE MD simulation itself (as that is the responsibility of the ASE library). Instead, we will mock the `ProcessPoolExecutor` to verify that the engine attempts to parallelize the simulations correctly. We will also mock the ASE `VelocityVerlet` dynamics object to assert that it is being instantiated with the correct parameters (e.g., temperature, timestep) derived from the Pydantic configuration model. This approach tests the *logic* of the engine, not the underlying physics simulation.

Finally, the `core.orchestrator.PipelineOrchestrator` will be tested using mock objects for each of the four pipeline stages. We will inject these mocks and then call the `run()` method. Assertions will then check that the methods on the mock objects (e.g., `mock_generator.generate()`, `mock_explorer.run_simulations()`) were called exactly once and in the correct order. This validates the orchestration logic without the overhead of running the actual pipeline.

**Integration Testing Approach (Min 300 words):**

The integration testing approach for Cycle 1 focuses on verifying that the independently tested units work together correctly to form a complete, functional pipeline. While unit tests ensure the components are correct in isolation, integration tests ensure that the "glue" between them is also correct. The primary goal is to simulate a real user scenario from the command line to the final database output.

The cornerstone of this strategy will be a single, comprehensive end-to-end test. This test will be located in the `tests/integration` directory. It will start by creating a temporary directory for the test run to isolate its outputs. Inside this directory, it will create a minimal but valid `config.yml` file. This configuration will be designed for speed: a very small system (e.g., 2 atoms of Cu), a small number of seed structures, and a very short MD simulation (e.g., 5-10 steps).

The test will then use the `typer.testing.CliRunner` to invoke the main CLI command (`ac-cdd run --config-path /path/to/temp/config.yml`). The test will execute the *real* code path, without mocks for the core pipeline components. This means it will actually generate structures, run a brief EMT simulation in parallel processes, sample the results, and write to a database file.

After the CLI process completes successfully (asserting a zero exit code), the test will proceed to the verification phase. It will locate the output ASE database file (`results.db`) in the temporary directory. Using `ase.db.connect`, it will open this database. The core assertions will then verify the integrity of the output:
1.  Check that the number of rows in the database matches the `num_samples` specified in the configuration.
2.  Read one of the rows (atomic configurations) from the database.
3.  Verify that it contains the expected data, such as `energy`, `forces`, and that the atoms are of the correct species.

This single test provides immense value by confirming that the data flows correctly through all four stages of the pipeline, that file I/O works as expected, and that the CLI can successfully orchestrate the entire process. It is the ultimate guarantee that the MVP is functional.
