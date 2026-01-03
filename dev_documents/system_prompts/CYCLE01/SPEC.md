# Specification: Cycle 1 - Core CLI Pipeline

## 1. Summary

This document provides the detailed technical specification for Cycle 1 of the MLIP-AutoPipe project. The primary objective of this cycle is to deliver a functional, end-to-end command-line interface (CLI) tool that automates the foundational pipeline for generating MLIP training data. This includes parsing a user-defined configuration file, generating initial atomic structures for a binary alloy, performing a simplified Molecular Dynamics (MD) simulation to explore new configurations, sampling structures from the resulting trajectory, and storing the final dataset into an ASE-compatible database.

This cycle establishes the architectural backbone of the entire application. We will focus on creating a robust, modular, and testable codebase. Key technical deliverables include the implementation of Pydantic models for strict configuration validation, a `click`-based CLI for user interaction, an abstract base class for structure generators with a concrete implementation for alloys, a basic MD exploration engine using a standard ASE potential, a simple random sampling module, and a database wrapper for data persistence. By the end of this cycle, a user will be able to define a simple binary alloy system in a YAML file and run a single command to produce a database of atomic structures, proving the viability of the core concept and providing a solid foundation for the advanced features planned in Cycle 2. The emphasis is on correctness, modularity, and establishing a clear workflow, rather than on computational performance or advanced scientific features, which will be the focus of the subsequent cycle.

## 2. System Architecture

The primary focus of Cycle 1 is to construct the foundational architecture of the MLIP-AutoPipe application. This involves establishing a robust, modular, and scalable file structure, implementing the core components for each stage of the pipeline, and ensuring they are correctly orchestrated. The architecture is designed with a strong emphasis on the separation of concerns, ensuring that each module has a single, well-defined responsibility. This approach not only facilitates parallel development and unit testing but also makes the system easier to maintain and extend in future cycles. All components will be created within the `src/mlip_autopipec/` package, adhering to standard Python packaging practices.

**File Structure (Cycle 1):**

The file structure for Cycle 1 is designed to logically group related functionalities. All files and directories listed below are to be created in this cycle.

```
src/mlip_autopipec/
├── __init__.py
├── **cli.py**                # The user-facing command-line interface entry point.
├── **config/**
│   ├── **__init__.py**
│   └── **models.py**         # Pydantic models for strict configuration validation.
├── **database/**
│   ├── **__init__.py**
│   └── **ase_db.py**         # A dedicated wrapper for all database interactions.
├── **generators/**
│   ├── **__init__.py**
│   ├── **base.py**           # Abstract base class defining the generator interface.
│   └── **alloy.py**          # Concrete implementation for generating alloy structures.
├── **explorers/**
│   ├── **__init__.py**
│   └── **md_engine.py**      # Logic for running basic Molecular Dynamics simulations.
├── **samplers/**
│   ├── **__init__.py**
│   ├── **base.py**           # Abstract base class defining the sampler interface.
│   └── **random_sampler.py** # Concrete implementation for random sampling.
├── **pipeline/**
│   ├── **__init__.py**
│   └── **orchestrator.py**   # The central orchestrator that connects all pipeline stages.
└── **shared/**
    ├── **__init__.py**
    └── **physics.py**        # Reusable functions for physics-based validation.
```

**Component Blueprints:**

This section provides a detailed blueprint for the classes and functions within each module, outlining their responsibilities, key methods, and interactions.

**`config/models.py`**: This module is the cornerstone of the application's reliability. It leverages Pydantic to define a strict, type-annotated schema for all user-configurable parameters. This schema-first approach ensures that any invalid configuration is caught at the earliest possible moment, preventing runtime errors deep within the computational pipeline. The `FullConfig` model will act as the root of the configuration tree, composing other, more specific models like `SystemConfig`, `ExplorationConfig`, and `SamplingConfig`. This nested structure improves readability and allows specific configuration objects to be passed to the components that need them, adhering to the principle of least knowledge. For instance, the `AlloyGenerator` will only receive the `SystemConfig` object, not the entire configuration.

**`database/ase_db.py`**: To decouple the application logic from the specifics of data persistence, all database operations will be encapsulated within an `AseDBWrapper` class. This class will act as a dedicated data access layer. It will manage the connection to the SQLite database via `ase.db.connect` and provide a clean, high-level API for the rest of the application. Its primary method, `write_atoms(atoms_list)`, will accept a list of `ase.Atoms` objects and iterate through them, writing each one to the database. This wrapper simplifies the logic in the main orchestrator and makes testing easier, as the entire database layer can be mocked by patching this single class.

**`generators/`**: This package introduces the factory pattern for creating initial atomic structures. The `base.py` module will define an abstract base class, `BaseStructureGenerator`, with a single abstract method: `generate() -> List[ase.Atoms]`. This enforces a common interface that all future generator classes must adhere to. For Cycle 1, the `alloy.py` module will provide the first concrete implementation, `AlloyGenerator`. This class will be initialized with a `SystemConfig` object. Its `generate` method will contain the logic for creating a disordered binary alloy. This process will involve several steps: starting with a primitive crystal structure (e.g., face-centered cubic) using `ase.build.bulk`, creating a sufficiently large supercell using `ase.build.make_supercell` to avoid periodic boundary issues, and then iterating through the atoms in the supercell to randomly assign their chemical symbols based on the specified composition. Crucially, after generation, each structure will be validated using functions from the `shared/physics.py` module to ensure it is physically plausible before it is returned.

**`explorers/md_engine.py`**: This module will house the computational core of the exploration stage. For Cycle 1, it will contain a function, `run_md(atoms, exploration_config)`, which performs a basic Molecular Dynamics simulation. This function will be responsible for setting up and running the simulation for a single input `ase.Atoms` object. The steps are as follows: first, an ASE-compatible calculator, such as `ase.calculators.emt.EMT`, is instantiated and attached to the `atoms` object. Second, initial velocities are assigned to the atoms from a Maxwell-Boltzmann distribution corresponding to the user-defined temperature. Third, an NVT thermostat, specifically the `ase.md.langevin.Langevin` dynamics algorithm, is attached to maintain the system's temperature. Finally, the simulation is executed for the specified number of steps, and the entire trajectory (a list of `ase.Atoms` objects, one for each timestep) is collected and returned.

**`samplers/`**: Similar to the generators, the samplers will follow a strategy pattern defined by an abstract base class, `BaseSampler`, in `base.py`. This ensures that different sampling algorithms can be used interchangeably by the orchestrator. In `random_sampler.py`, the `RandomSampler` class will provide the baseline implementation. It will be initialized with a `SamplingConfig` object and its `sample(trajectory)` method will simply use Python's built-in `random.sample` function to select a random subset of the structures from the input trajectory list. While simple, this provides a complete, end-to-end pipeline and serves as a benchmark against which the more advanced samplers of Cycle 2 can be compared.

**`pipeline/orchestrator.py`**: This module is the central nervous system of the application. The `WorkflowOrchestrator` class is responsible for executing the entire pipeline from start to finish. It is initialized with the validated `FullConfig` object. Its main public method, `run()`, will execute the four stages in a strict sequence:
1.  **Generation:** It will instantiate the `AlloyGenerator` with the `system` configuration and call its `generate` method to produce the initial seed structures.
2.  **Exploration:** It will iterate through each of the seed structures, calling the `run_md` function from the `md_engine` for each one. The resulting trajectories will be aggregated into a single, large list of all atomic configurations generated.
3.  **Sampling:** It will instantiate the `RandomSampler` with the `sampling` configuration and pass the aggregated trajectory to its `sample` method to obtain the final, curated list of structures.
4.  **Storage:** Finally, it will instantiate the `AseDBWrapper` with the `db_path` and call its `write_atoms` method to persist the final sampled structures to the database.
By encapsulating this control flow, the `WorkflowOrchestrator` decouples the individual components, which do not need to be aware of each other.

**`cli.py`**: This module will provide the user's primary entry point to the application. It will be implemented using the `click` library to create a clean and well-documented command-line interface. A main command, `run-pipeline`, will be defined. This command will use the `hydra` library to parse the user-provided YAML configuration file into the `FullConfig` Pydantic object. This powerful combination provides both a user-friendly CLI and robust, automatic configuration validation. Upon successfully parsing the configuration, this function will instantiate the `WorkflowOrchestrator` with the config object and invoke its `run` method to start the pipeline. It will also be responsible for catching any exceptions that bubble up from the pipeline and presenting them to the user as clear, informative error messages.

## 3. Design Architecture

This project's design is schema-first, driven by Pydantic. The configuration models are the canonical source of truth for all parameters used in the system, ensuring type safety and validation at the entry point of the application.

**Pydantic Schema Design (`config/models.py`):**

The configuration will be composed of several nested Pydantic models to ensure clarity and modularity.

```python
# config/models.py

from pydantic import BaseModel, Field
from typing import List, Dict

class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""
    elements: List[str] = Field(..., min_length=1, description="List of element symbols (e.g., ['Fe', 'Pt']).")
    composition: Dict[str, float] = Field(..., description="Composition of each element (e.g., {'Fe': 0.5, 'Pt': 0.5}).")
    num_structures: int = Field(10, gt=0, description="Number of initial seed structures to generate.")

class ExplorationConfig(BaseModel):
    """Configuration for the exploration (MD) stage."""
    temperature_k: float = Field(300.0, gt=0, description="Simulation temperature in Kelvin.")
    num_steps: int = Field(1000, gt=0, description="Number of MD steps to perform.")
    potential: str = Field("EMT", description="ASE-compatible potential to use (e.g., EMT).")

class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""
    method: str = Field("random", description="Sampling method. Only 'random' is supported in Cycle 1.")
    num_samples: int = Field(100, gt=0, description="Number of structures to sample from the trajectory.")

class FullConfig(BaseModel):
    """The root configuration model."""
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    db_path: str = Field("mlip_data.db", description="Path to the output ASE database file.")

```
**Key Invariants and Constraints:**
- The sum of values in `SystemConfig.composition` must equal 1.0. This will be enforced with a Pydantic root validator.
- `SystemConfig.elements` and the keys of `SystemConfig.composition` must contain the same set of elements.
- The `ExplorationConfig.potential` must be a calculator known to ASE. While not validated by Pydantic, the application logic will handle errors for unknown potentials.

**Consumers and Producers:**
- **Producer:** The user creates a YAML file that conforms to this schema.
- **Consumer:** The `cli.py` module, using Hydra, will be the primary consumer. It will parse the YAML file into a `FullConfig` object. This object will then be passed to the `WorkflowOrchestrator`, which will propagate the relevant sub-configs to the various pipeline components (e.g., `SystemConfig` is passed to the `AlloyGenerator`). This ensures that each component only has access to the configuration it needs.

This schema-driven approach prevents a large class of runtime errors by catching invalid configurations before any computationally expensive tasks are started.

## 4. Implementation Approach

The implementation will proceed in a logical order, building from the foundational components to the final orchestrator.

1.  **Pydantic Models:** First, implement all the Pydantic models as defined in the Design Architecture section in `src/mlip_autopipec/config/models.py`. This is the contract for the entire application.
2.  **Database Wrapper:** Implement the `AseDBWrapper` in `src/mlip_autopipec/database/ase_db.py`. It should have methods like `__init__(self, path)`, `connect(self)`, and `write_atoms(self, atoms_list: list)`. The implementation will be a thin wrapper around `ase.db.connect`.
3.  **Generator Interface and Implementation:**
    -   In `src/mlip_autopipec/generators/base.py`, define an abstract base class `BaseStructureGenerator` with an abstract method `generate(self) -> List[ase.Atoms]`.
    -   In `src/mlip_autopipec/generators/alloy.py`, create the `AlloyGenerator` class that inherits from `BaseStructureGenerator`. Its `__init__` will accept a `SystemConfig` object. The `generate` method will use `ase.build.bulk` to create a primitive cell, `ase.build.make_supercell` to expand it, and then randomly assign atomic symbols based on the composition.
4.  **Physics Validation:** In `src/mlip_autopipec/shared/physics.py`, implement a function `check_interatomic_distance(atoms: ase.Atoms, min_dist: float) -> bool` that returns `False` if any two atoms are closer than `min_dist`. This will be used by the generator.
5.  **Explorer Module:** In `src/mlip_autopipec/explorers/md_engine.py`, create a simple `run_md` function that takes an `ase.Atoms` object and an `ExplorationConfig` object. It will attach an ASE calculator (e.g., `ase.calculators.emt.EMT`), set up the MD dynamics (e.g., `ase.md.velocitydistribution.MaxwellBoltzmannDistribution`, `ase.md.langevin.Langevin`), and run the simulation, collecting the trajectory.
6.  **Sampler Interface and Implementation:**
    -   In `src/mlip_autopipec/samplers/base.py`, define a `BaseSampler` with an abstract method `sample(self, trajectory: List[ase.Atoms]) -> List[ase.Atoms]`.
    -   In `src/mlip_autopipec/samplers/random_sampler.py`, create a `RandomSampler` class. Its `__init__` will take a `SamplingConfig` object. The `sample` method will simply use Python's `random.sample` to select frames from the trajectory.
7.  **Orchestrator:** In `src/mlip_autopipec/pipeline/orchestrator.py`, create the `WorkflowOrchestrator`. Its `__init__` will take a `FullConfig` object. It will have a `run()` method that executes the full pipeline in sequence:
    -   Instantiate `AlloyGenerator` with `config.system`.
    -   Call `generator.generate()` to get seed structures.
    -   Loop through each seed structure and call `run_md()` from the explorer module.
    -   Aggregate all trajectories.
    -   Instantiate `RandomSampler` with `config.sampling`.
    -   Call `sampler.sample()` on the aggregated trajectory.
    -   Instantiate `AseDBWrapper` with `config.db_path`.
    -   Call `db_wrapper.write_atoms()` to save the final structures.
8.  **CLI Entry Point:** Finally, in `src/mlip_autopipec/cli.py`, create a `click` command. This command will use Hydra to load the YAML configuration file, instantiate the `WorkflowOrchestrator`, and call its `run()` method.

## 5. Test Strategy

Testing in Cycle 1 is crucial to validate the architecture and ensure each component behaves as expected.

**Unit Testing Approach:**
Each module will be tested in isolation using `pytest`. Mocks will be used extensively to isolate components from their dependencies.
-   **`config/models.py`:** Test that `FullConfig` correctly parses valid YAML files. Write specific tests that assert `pydantic.ValidationError` is raised for invalid inputs (e.g., `temperature_k = -100`, composition sum not equal to 1.0).
-   **`generators/alloy.py`:** Test the `AlloyGenerator`. Given a fixed `SystemConfig`, assert that the output is a list of `ase.Atoms` objects of the correct number, size, and chemical composition. Mock the `ase.build` functions to ensure reproducibility.
-   **`database/ase_db.py`:** Test the `AseDBWrapper` against a temporary database file created in the test directory. Write a list of atoms, then read them back and assert that the returned objects are identical to the originals. Ensure the test cleans up the temporary file.
-   **`samplers/random_sampler.py`:** Test the `RandomSampler`. Provide a list of 10 atoms and ask for 5 samples. Assert that the returned list has 5 elements and that all returned elements were present in the original list.
-   **`explorers/md_engine.py`:** This is the most complex unit to test. The `run_md` function will be tested by heavily mocking the ASE dynamics and calculator objects. The goal is not to test ASE's MD implementation, but to ensure our function correctly sets up the calculator, attaches it to the atoms object, and calls the `dyn.run()` method with the correct number of steps.

**Integration Testing Approach:**
A single, powerful integration test will be created to verify the end-to-end workflow of the CLI.
-   The test will live in `tests/integration/test_cli.py`.
-   It will use `click.testing.CliRunner` to invoke the main CLI command programmatically.
-   A minimal, valid `config.yaml` file will be created in a temporary directory managed by the test runner.
-   The most computationally expensive part, `explorers.md_engine.run_md`, will be mocked using `pytest-mock`. The mock will be configured to return a small, pre-defined trajectory of `ase.Atoms` objects immediately, without actually running any simulation.
-   The test will invoke the CLI runner, pointing to the temporary config file.
-   After the command finishes, the test will assert the following:
    1.  The command exited with a success code (0).
    2.  The output database file (e.g., `mlip_data.db`) was created at the expected location.
    3.  The test will connect to this database and verify that it contains the correct number of structures, as specified in the `sampling.num_samples` section of the config file.
    4.  It will verify that the structures in the database are a subset of the ones returned by the mocked `run_md` function.

This integration test ensures that all the components, from the CLI parsing to the final database write, are correctly wired together in the `WorkflowOrchestrator` and function as a cohesive whole.
