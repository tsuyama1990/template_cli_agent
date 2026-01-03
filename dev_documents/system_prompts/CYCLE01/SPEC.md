# CYCLE01 SPECIFICATION: Core Pipeline and CLI Foundation

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to build a robust, end-to-end, command-line-driven pipeline capable of generating a foundational dataset for a single class of materials: multi-component alloys. This cycle will focus on establishing the core architectural pillars of the application, including the main pipeline orchestrator, the data persistence layer using an ASE database, a modular structure for pipeline components, and a type-safe configuration system based on Pydantic and Hydra.

By the end of this cycle, a user will be able to define an alloy system in a configuration file and run a single command to execute the full, four-stage workflow:
1.  **Generation:** Create initial, physically valid alloy structures.
2.  **Exploration:** Run basic Molecular Dynamics (MD) simulations to generate diverse atomic configurations.
3.  **Sampling:** Select a subset of these configurations using a simple random sampling method.
4.  **Storage:** Save the final, selected structures into an ASE database.

The emphasis is on creating a solid, testable, and extensible foundation. The components built in this cycle, such as the `PipelineOrchestrator` and the `AseDBWrapper`, are designed to be generic. The exploration and sampling capabilities will be basic, serving as placeholders for the more advanced engines that will be developed in Cycle 2. The successful completion of this cycle will yield a functional, albeit simple, version of the tool that proves the viability of the overall architecture and provides a stable base for future development. The core deliverables include a working CLI, a well-defined project structure, and a comprehensive suite of unit and integration tests to ensure the reliability of the foundational code.

## 2. System Architecture

The architecture for Cycle 1 establishes the fundamental file and module structure for the entire project. The focus is on creating the essential components for the core pipeline. Bolded files and directories are the ones to be created or significantly modified in this cycle.

```
src/
└── mlip_autopipec/
    ├── __init__.py
    ├── cli/
    │   ├── __init__.py
    │   └── main.py              # Main CLI entrypoint (using Click)
    ├── core/
    │   ├── __init__.py
    │   ├── factories.py         # Factory classes for dispatching implementations
    │   ├── interfaces.py        # Abstract Base Classes defining component interfaces
    │   └── pipeline_orchestrator.py # The core PipelineRunner logic
    ├── domain/
    │   ├── __init__.py
    │   ├── configuration.py     # Pydantic models for Hydra configuration
    │   └── data_models.py       # Pydantic models for internal data structures (e.g., results)
    ├── services/
    │   ├── __init__.py
    │   ├── generation/
    │   │   ├── __init__.py
    │   │   ├── alloy.py         # **AlloyGenerator implementation**
    │   │   └── base.py          # **BaseGenerator abstract class**
    │   ├── exploration/
    │   │   ├── __init__.py
    │   │   └── md_engine.py     # **Basic MDEngine implementation**
    │   ├── sampling/
    │   │   ├── __init__.py
    │   │   └── random.py        # **RandomSampler implementation**
    │   └── storage/
    │       ├── __init__.py
    │       └── ase_db_wrapper.py # **ASE database wrapper**
    └── utils/
        ├── __init__.py
        └── physics.py           # Physics validation and utility functions

tests/
├── __init__.py
├── integration/
│   └── test_cli_pipeline.py # **End-to-end pipeline test**
└── unit/
    ├── __init__.py
    ├── services/
    │   ├── generation/
    │   │   └── test_alloy_generator.py # **Tests for AlloyGenerator**
    │   └── storage/
    │       └── test_ase_db_wrapper.py  # **Tests for AseDBWrapper**
    └── domain/
        └── test_configuration.py   # **Tests for Pydantic models**
```

This structure clearly separates concerns:
-   `domain`: Defines all data structures, ensuring a consistent data model throughout the application.
-   `core`: Contains the high-level business logic and abstract interfaces, decoupling the orchestrator from specific implementations.
-   `services`: Holds the concrete implementations of the interfaces for each pipeline stage.
-   `cli`: Provides the user-facing command-line interface.
-   `tests`: A parallel structure for unit and integration tests, ensuring code quality and correctness.

## 3. Design Architecture

The design of MLIP-AutoPipe is centered around a schema-first approach using Pydantic, which enforces strict data contracts between components. This ensures that data is always valid as it moves through the pipeline, reducing runtime errors and making the system more robust.

**Pydantic-Based Configuration (`domain/configuration.py`):**
The entire application will be configured through Pydantic models, which are then managed by Hydra. This provides a single source of truth for all parameters.

```python
# Example Pydantic models
from pydantic import BaseModel, Field
from typing import List, Dict

class AlloyGeneratorConfig(BaseModel):
    elements: Dict[str, float] = Field(..., description="Element symbols and their fractions, e.g., {'Cu': 0.5, 'Au': 0.5}")
    num_structures: int = Field(10, gt=0, description="Number of initial structures to generate")

class MDExplorerConfig(BaseModel):
    temperature_k: float = Field(300.0, gt=0, description="Simulation temperature in Kelvin")
    num_steps: int = Field(1000, gt=0, description="Number of MD steps")
    mlip_model_path: str = "path/to/model.pt"

class PipelineConfig(BaseModel):
    generator: AlloyGeneratorConfig
    explorer: MDExplorerConfig
    # ... other configs
```
**Key Invariants and Constraints:**
-   The sum of element fractions in `AlloyGeneratorConfig` must be validated to be close to 1.0.
-   All numerical parameters like `num_structures` and `temperature_k` must be strictly positive.
-   File paths for models will be validated to ensure they exist.

**Consumers and Producers:**
-   **Producer:** The user creates a `config.yaml` file on disk.
-   **Consumer:** Hydra, guided by the Pydantic models, loads and validates this configuration at the start of the `cli/main.py` execution. The resulting Pydantic object is then passed to the `PipelineOrchestrator`.

**Component Interfaces (`core/interfaces.py`):**
To ensure modularity and extensibility, the core pipeline components will be defined by abstract base classes (interfaces).

```python
from abc import ABC, abstractmethod
from typing import List
from ase import Atoms

class IStructureGenerator(ABC):
    @abstractmethod
    def generate(self) -> List[Atoms]:
        """Generates a list of initial atomic structures."""
        pass

class IExplorer(ABC):
    @abstractmethod
    def run_explorer(self, structures: List[Atoms]) -> None:
        """Runs the exploration process on the given structures."""
        pass

# ... other interfaces for Sampler, Storage
```
These interfaces act as a formal contract. The `PipelineOrchestrator` will be programmed against these interfaces, not the concrete implementations. This means we can easily add a new type of generator in the future (e.g., `IonicGenerator`) without changing the orchestrator's code, as long as the new generator implements the `IStructureGenerator` interface. This is the cornerstone of the system's extensibility.

## 4. Implementation Approach

The implementation will proceed in a logical, bottom-up sequence, starting with data structures and moving towards the high-level orchestration and user interface.

1.  **Project Scaffolding:** Create the directory structure as outlined in the System Architecture section. Initialise a `pyproject.toml` file with `uv`, and add initial dependencies like `click`, `pydantic`, `hydra-core`, `ase`, and `pytest`.

2.  **Domain Models:** Implement all configuration Pydantic models in `domain/configuration.py`. This includes creating models for each pipeline stage (Generator, Explorer, Sampler, Storage) and a main model that composes them. Write unit tests to verify the validation logic.

3.  **Storage Layer:** Implement the `AseDBWrapper` in `services/storage/ase_db_wrapper.py`. This class will encapsulate all interactions with the `ase.db` SQLite database. It will provide methods like `connect()`, `write_atoms()`, and `get_atoms_count()`. This abstraction is crucial for isolating the rest of the application from the specifics of the database implementation.

4.  **Component Interfaces:** Define the abstract base classes (`IStructureGenerator`, `IExplorer`, `ISampler`, `IStorage`) in `core/interfaces.py`.

5.  **Service Implementations (The Pipeline Stages):**
    *   **Generation:** Implement the `AlloyGenerator` in `services/generation/alloy.py`. It will take an `AlloyGeneratorConfig` Pydantic model. The logic will use ASE's `ase.build.bulk` and randomization functions to create alloy structures. It will also incorporate the physical validation utilities from `utils/physics.py`.
    *   **Exploration:** Implement the `MDEngine` in `services/exploration/md_engine.py`. It will use ASE's `Langevin` or other dynamics integrators. The MLIP calculator (e.g., MACE) will be attached to the `Atoms` objects to compute forces. The engine will be configured via the `MDExplorerConfig`.
    *   **Sampling:** Implement the `RandomSampler` in `services/sampling/random.py`. This service will read the trajectory files produced by the `MDEngine` and randomly select a specified number of frames.

6.  **Orchestration:** Implement the `PipelineOrchestrator` in `core/pipeline_orchestrator.py`. This is the heart of the application. Its `run()` method will execute the four stages in sequence. It will handle the I/O between stages, such as creating directories for intermediate files and ensuring that the output of one stage is correctly passed as the input to the next.

7.  **CLI Entrypoint:** Implement the main CLI function in `cli/main.py` using `click`. This function will be decorated with `@hydra.main` to handle configuration loading. It will instantiate the `PipelineOrchestrator` and the necessary service components (using the factories) and then call the orchestrator's `run()` method.

## 5. Test Strategy

Testing in Cycle 1 is critical for building a reliable foundation. We will use a combination of unit and integration tests.

**Unit Testing Approach (Min 300 words):**
The primary goal of unit testing in this cycle is to verify the correctness of each component in complete isolation. We will use the `pytest` framework and extensive mocking to achieve this. For the `AlloyGenerator`, tests will be written to confirm that it generates the correct number of structures, that the chemical composition of the generated structures matches the input configuration, and that basic physical constraints (like minimum atom distances) are respected. We will not need to run a real simulation; we only need to check the properties of the `Atoms` objects it returns.

For the `AseDBWrapper`, unit tests will use an in-memory SQLite database (`":memory:"`) to verify all CRUD (Create, Read, Update, Delete) operations. We will test that it correctly writes `Atoms` objects, reads them back, and handles edge cases like an empty database or corrupted entries gracefully.

The `PipelineOrchestrator` will be tested by providing it with mocked versions of the generator, explorer, sampler, and storage services. For example, we can use `unittest.mock.MagicMock` to create a fake generator object. The test will then assert that the orchestrator calls the `generate()` method on the mock object and that it correctly handles the list of `Atoms` objects that the mock returns. This allows us to test the orchestration logic (the "glue" code) without needing to run the actual, time-consuming pipeline stages. Similarly, Pydantic configuration models will be tested to ensure their validation logic works as expected, for instance by attempting to create a model with invalid data and asserting that a `ValidationError` is raised.

**Integration Testing Approach (Min 300 words):**
While unit tests verify components in isolation, integration tests are essential to ensure they work together correctly. The centerpiece of our integration testing strategy for Cycle 1 will be an end-to-end test of the CLI itself, located in `tests/integration/test_cli_pipeline.py`. This test will use `click.testing.CliRunner` to invoke the main CLI command, simulating a real user running the application.

To make the test fast and deterministic, we will use a simplified configuration. Instead of a computationally expensive MLIP like MACE, the test will configure the `MDEngine` to use ASE's built-in and very fast Effective Medium Theory (`EMT`) calculator. The test will run a complete, albeit miniature, pipeline:
1.  **Generate:** Create 1-2 small alloy structures.
2.  **Explore:** Run a very short MD simulation (e.g., 10 steps).
3.  **Sample:** Randomly select 5 frames.
4.  **Store:** Write the final structures to a temporary database file.

The test will be run within an isolated filesystem to avoid cluttering the project with test artifacts. After the CLI command finishes, the test will use the `AseDBWrapper` to connect to the temporary database and assert that the correct number of final structures have been saved. It will also perform basic checks on the stored structures, such as verifying their energy is a finite number. This single test provides high confidence that the main data flow, the interaction between components, the configuration loading, and the file I/O are all functioning correctly as a cohesive system.
