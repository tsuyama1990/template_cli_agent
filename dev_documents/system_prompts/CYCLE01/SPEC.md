# SPECIFICATION: Cycle 1 - Core CLI Pipeline and Foundational Components

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to build the foundational, non-graphical components of the application and deliver a fully functional, end-to-end command-line interface (CLI) pipeline. This cycle is fundamentally about establishing a robust and extensible architectural backbone for the entire system. It will focus on creating a clear and logical data flow, implementing a schema-driven configuration system for unparalleled robustness, and developing the core stages of the workflow: Generation, Exploration, Sampling, and Storage. By prioritizing the CLI, we cater first to the power-user persona—the computational scientist who will integrate this tool into larger, scripted workflows in high-performance computing (HPC) environments.

By the end of this cycle, a user will be able to define a complex physical system, such as a multi-component alloy or an ionic compound, within a structured YAML configuration file. They will then be able to execute the pipeline via a single, clear CLI command and, upon completion, receive a valid, ready-to-use ASE database containing a set of generated atomic structures. A key strategic decision for this cycle is to implement the exploration phase in a simplified manner. Instead of integrating a full, computationally expensive MLIP model, this initial version will use a placeholder calculator (like ASE's built-in EMT or Lennard-Jones potentials). This approach is critical for rapid, iterative development; it allows us to validate the entire pipeline's structural integrity, data flow, parallel processing logic, and checkpointing capabilities without the significant overhead and complexity of the final MLIP-based engine. This de-risks the project by ensuring the fundamental architecture is sound before we introduce the more complex scientific components in Cycle 2. Key deliverables for this cycle are therefore not just the final CLI tool, but also the critical architectural assets: the comprehensive Pydantic-based configuration system that guarantees type safety and user input validation; a set of abstract base classes for core components like generators and samplers, which ensures future extensibility; a suite of concrete implementations for common use-cases, including alloy and ionic structure generators; a basic random sampler; and the essential database writing module. This foundational work is paramount, providing a stable and well-tested platform upon which the advanced, intelligent features of the subsequent cycle can be confidently built.

## 2. System Architecture

The architecture for Cycle 1 is meticulously designed to establish the core modules and their interactions, laying the groundwork for all future development. The paramount focus is on creating a modular and decoupled system where each component has a single, well-defined responsibility, thereby adhering to sound software engineering principles. This separation of concerns is what will allow the system to be tested, maintained, and extended efficiently. The file structure is organized to reflect this modularity, with distinct directories for generators, explorers, samplers, and other primary components.

**File Structure (Cycle 1 Focus):**

The files and directories to be created or modified in this cycle are marked in **bold**. This structure clearly separates the application's internal logic (`mlip_autopipec`) from its entry points (`main_cli.py`) and common data structures (`common/`).

```
.
├── src/
│   └── mlip_autopipec/
│       ├── **__init__.py**
│       ├── **main_cli.py**         # Entry point for the command-line interface
│       ├── pipeline/
│       │   ├── **__init__.py**
│       │   └── **runner.py**       # Contains the main PipelineRunner orchestrator
│       ├── generators/
│       │   ├── **__init__.py**
│       │   ├── **base.py**         # Defines the BaseGenerator abstract class
│       │   ├── **factory.py**      # Implements the GeneratorFactory
│       │   ├── **alloy.py**        # Concrete implementation: AlloyGenerator
│       │   └── **ionic.py**        # Concrete implementation: IonicGenerator
│       ├── explorers/
│       │   ├── **__init__.py**
│       │   └── **md_engine.py**    # A simplified MD engine for this cycle
│       ├── samplers/
│       │   ├── **__init__.py**
│       │   ├── **base.py**         # Defines the BaseSampler abstract class
│       │   └── **random_sampler.py** # Concrete implementation: RandomSampler
│       ├── storage/
│       │   ├── **__init__.py**
│       │   └── **db_writer.py**    # The dedicated ASE DB writing component
│       └── common/
│           ├── **__init__.py**
│           ├── **config.py**       # All Pydantic models for configuration
│           └── **atoms_utils.py**  # Utility functions for ASE Atoms objects
└── tests/
    ├── **__init__.py**
    ├── **test_generators.py**
    ├── **test_pipeline.py**
    └── **test_config.py**
```

**Code Blueprints:**

The following blueprints provide a more detailed look at the core classes and their interactions, forming the technical specification for the implementation.

*   **`pipeline/runner.py` - `PipelineRunner`:** This class is the central nervous system of the application. It is initialized with the fully validated configuration object and is responsible for the entire lifecycle of a run. It will instantiate the necessary components using factories, execute the four pipeline stages in strict sequence, and manage the persistence of intermediate results (checkpointing) to ensure fault tolerance.
    ```python
    from mlip_autopipec.common.config import MainConfig
    from mlip_autopipec.generators.factory import GeneratorFactory
    # ... other imports

    class PipelineRunner:
        def __init__(self, config: MainConfig):
            self.config = config
            # Use the factory to decouple the runner from concrete generator classes
            self.generator = GeneratorFactory.create(config.system)
            self.explorer = MDExplorer(config.exploration)
            self.sampler = RandomSampler(config.sampling)
            self.db_writer = DBWriter(config.storage)
            # Add logging configuration here

        def run(self):
            # Stage 1: Generation
            # logger.info("Starting structure generation...")
            initial_structures = self.generator.generate()
            # self._save_checkpoint(1, initial_structures)

            # Stage 2: Exploration
            # logger.info("Starting exploration...")
            trajectories = self.explorer.run(initial_structures)
            # self._save_checkpoint(2, trajectories)

            # Stage 3: Sampling
            # logger.info("Starting sampling...")
            sampled_structures = self.sampler.sample(trajectories)
            # self._save_checkpoint(3, sampled_structures)

            # Stage 4: Storage
            # logger.info("Starting database storage...")
            self.db_writer.write(sampled_structures)
            # logger.info("Pipeline finished successfully.")
    ```

*   **`generators/base.py` - `BaseGenerator`:** This abstract base class defines the "contract" that all generator classes must follow. It ensures that every generator provides a `generate` method and has access to a shared, robust validation method.
    ```python
    from abc import ABC, abstractmethod
    from ase import Atoms
    from mlip_autopipec.common.config import SystemConfig

    class BaseGenerator(ABC):
        def __init__(self, config: SystemConfig):
            self.config = config

        @abstractmethod
        def generate(self) -> list[Atoms]:
            """Generates a list of initial seed structures."""
            raise NotImplementedError

        def _validate_structure(self, atoms: Atoms, min_dist_factor: float = 0.7) -> bool:
            """
            Performs physics-based validation on a structure.
            Checks for atomic overlaps based on covalent radii.
            Returns True if valid, False otherwise.
            """
            # Robust implementation for checking interatomic distances
            # using scaled covalent radii to prevent unphysical structures.
            pass
    ```

*   **`common/config.py` - Pydantic Schemas:** This is the schema-first heart of the application. It provides a single, unambiguous source of truth for all configuration parameters, complete with validation rules and documentation.
    ```python
    from pydantic import BaseModel, Field, validator
    from typing import List, Dict, Literal

    class SystemConfig(BaseModel):
        elements: List[str] = Field(..., min_length=1, description="List of elements, e.g., ['Cu', 'Au']")
        composition: Dict[str, float]
        generator_type: Literal["alloy", "ionic"]
        num_initial_structures: int = Field(..., gt=0, description="Number of seed structures to generate.")

        @validator('composition')
        def composition_must_sum_to_one(cls, v):
            if not abs(sum(v.values()) - 1.0) < 1e-6:
                raise ValueError("Composition fractions must sum to 1.0")
            return v
        # ... more validation logic

    class ExplorationConfig(BaseModel):
        md_steps: int = Field(1000, gt=0)
        temperature_k: float = Field(300.0, gt=0)
        # ...

    class SamplingConfig(BaseModel):
        num_samples: int = Field(..., gt=0)

    class StorageConfig(BaseModel):
        db_path: str = Field("final_dataset.db")

    class MainConfig(BaseModel):
        system: SystemConfig
        exploration: ExplorationConfig
        sampling: SamplingConfig
        storage: StorageConfig
    ```

## 3. Design Architecture

The design for Cycle 1 is fundamentally centered around a schema-first, Pydantic-driven architecture. This design choice is a pre-emptive measure against a vast category of common runtime errors and significantly enhances the robustness and user-friendliness of the system. All configuration that directs the pipeline's behavior—from the types of atoms to use, to the number of MD steps, to the final database filename—will be meticulously defined in strictly-typed Pydantic models. These models, located in the central `mlip_autopipec/common/config.py` module, act as the definitive "contract" for the application's configuration.

The primary producer of these configuration objects is the end-user, who will author a human-readable YAML file. This file is then consumed by the Hydra library, which serves as the parser and instantiator, transforming the YAML into a validated `MainConfig` Pydantic object. This powerful combination provides immense flexibility; users can easily override specific parameters from the command line (e.g., `mlip-autopipec exploration.temperature_k=500`) for quick experiments, without ever modifying the base configuration file.

The primary consumer of this configuration object is the `PipelineRunner` class, which acts as the application's central orchestrator. Upon initialization, it receives the top-level `MainConfig` object. It then follows the principle of delegation, passing the relevant sub-configuration objects to the components it constructs. For example, the `SystemConfig` is passed to the `GeneratorFactory`, the `ExplorationConfig` to the `MDExplorer`, and so on. This ensures that each component only has access to the configuration it needs, adhering to the principle of least privilege and making the components more modular and easier to test.

The data flow within the pipeline is designed to be a simple, linear, and sequential process, which makes it easy to reason about and debug. Each major component (Generator, Explorer, Sampler, Writer) is designed to be stateless and to act as a pure data-transformation step. It receives data from the previous stage, performs its specific task, and returns the transformed data to the `PipelineRunner`. For instance, the `Generator` takes no input and produces a `list[ase.Atoms]`. This list is then passed as input to the `Explorer`, which in turn produces a `list[list[Atoms]]` (a list of trajectories). This strict, unidirectional data flow prevents complex, state-dependent bugs. The `PipelineRunner` is the sole manager of state (e.g., saving intermediate results to disk for checkpointing), which cleanly separates the core scientific logic of each component from the operational concerns of the pipeline itself. For this cycle, versioning of the configuration schema is straightforward, but in the future, this Pydantic-based design will allow for graceful schema evolution and migration if needed.

## 4. Implementation Approach

The implementation will be executed in a structured, bottom-up fashion, starting with the foundational data structures and progressively building up to the high-level orchestration and user interface layers. This ensures that each layer is built upon a stable, tested foundation.

1.  **Project Scaffolding:** The first step is to create the complete directory structure as outlined in the System Architecture section. This includes creating all the necessary `__init__.py` files to define the Python packages. The `pyproject.toml` file will be populated with the project's core dependencies for this cycle: `ase` (for atomic structures), `pydantic` (for configuration), `hydra-core` (for config loading), `typer` (for the CLI), and `pytest` (for testing). A `uv sync` command will then be run to establish the virtual environment.

2.  **Configuration Schema (`common/config.py`):** The implementation will begin with the most critical file, defining all the Pydantic models: `SystemConfig`, `ExplorationConfig`, `SamplingConfig`, `StorageConfig`, and the top-level `MainConfig`. This schema-first approach is crucial, as it defines the data contracts for the rest of the application. Custom validators, such as ensuring that composition percentages sum to 1.0, will be implemented here.

3.  **Generator Modules (`generators/`):**
    *   The abstract base class `BaseGenerator` will be created in `base.py`. It will define the abstract `generate` method and include a concrete, shared implementation of the `_validate_structure` method for checking atomic overlaps.
    *   The `AlloyGenerator` will be implemented in `alloy.py`. It will leverage ASE's `ase.build.bulk` and random atom swapping capabilities to create initial alloy structures based on the specified composition.
    *   The `IonicGenerator` will be implemented in `ionic.py` for creating simple, charge-neutral ionic compounds like NaCl.
    *   Finally, the `GeneratorFactory` will be implemented in `factory.py`. It will have a static `create` method that takes the `SystemConfig` and returns the correct generator instance based on the `generator_type` field.

4.  **Placeholder Explorer (`explorers/md_engine.py`):** A simplified `MDExplorer` will be implemented. It will not use a real MLIP in this cycle. Instead, it will be configured to use a simple, fast calculator from ASE, such as `EMT` or `LennardJones`. Its `run` method will iterate through the provided initial structures, attach the calculator, and run a short `VelocityVerlet` MD simulation for each. The purpose of this component in Cycle 1 is not scientific accuracy but to ensure that the data pipeline (receiving structures, running a simulation in parallel, and returning trajectories) is working correctly and efficiently.

5.  **Sampling Module (`samplers/`):**
    *   An abstract `BaseSampler` class will be created in `base.py`.
    *   A `RandomSampler` will be implemented in `random_sampler.py`. It will take the list of trajectories from the explorer, flatten them into a single list of all generated structures, and then randomly select a specified number of `Atoms` objects from this list.

6.  **Storage Module (`storage/db_writer.py`):** The `DBWriter` class will be implemented. It will expose a single public method, `write`, which takes a list of `Atoms` objects. Internally, it will use `ase.db.connect` to open a connection to the SQLite database file specified in the configuration and then iterate through the structures, writing each one to the database along with its energy and forces (if available from the calculator).

7.  **Orchestration (`pipeline/runner.py`):** The `PipelineRunner` class will be implemented to tie everything together. Its `__init__` method will take the `MainConfig` object and instantiate all the necessary components. Its `run` method will implement the main application logic, calling each component in sequence and passing the data from one stage to the next, with logging at each step.

8.  **CLI (`main_cli.py`):** Finally, the command-line interface will be created using `Typer`. It will define a single main command that accepts the path to a Hydra configuration file. The function will use Hydra to load and validate the config, instantiate the `PipelineRunner` with the resulting `MainConfig` object, and then call the `run` method, wrapping the entire process in a try/except block for graceful error handling.

## 5. Test Strategy

The testing strategy for Cycle 1 is focused on rigorously verifying the correctness of the core pipeline logic and the individual components, ensuring a stable foundation for future development.

**Unit Testing Approach:**

Unit tests are essential for ensuring that each component behaves correctly in isolation. We will use the `pytest` framework and `pytest-mock` library to achieve this.

*   **Configuration (`test_config.py`):** We will write tests specifically for our Pydantic models. These tests will attempt to create model instances with both valid and invalid data. For example, we will assert that a `ValidationError` is raised if `num_initial_structures` is set to `-1` or if the composition fractions do not sum to 1.0. This verifies our custom validation logic.
*   **Generators (`test_generators.py`):**
    *   The `AlloyGenerator` will be tested to ensure it creates the correct number of structures with the correct number of atoms.
    *   We will assert that the chemical composition of the generated alloy structures is statistically consistent with the input configuration.
    *   The `_validate_structure` method will be tested directly by manually creating an `Atoms` object with two overlapping atoms and asserting that the method correctly returns `False`.
    *   The `GeneratorFactory` will be tested to ensure it returns an instance of `AlloyGenerator` when the config specifies `"alloy"`.
*   **Sampler (`test_sampler.py`):** The `RandomSampler` will be tested by providing it with a dummy trajectory (a list of lists of `Atoms` objects) and asserting that it returns a list containing the correct number of samples as specified in its configuration.
*   **PipelineRunner (`test_pipeline.py`):** The orchestration logic of the `PipelineRunner` is a critical test target. We will use `pytest-mock` to replace the real Generator, Explorer, Sampler, and DBWriter with mock objects (`mocker.patch`). We will then call the `runner.run()` method and assert that the `generate`, `run`, `sample`, and `write` methods on our mock objects were called in the correct sequence. This verifies the control flow without needing to perform any expensive computations.

**Integration Testing Approach:**

Integration tests are designed to verify that the individual components, once wired together, function correctly as a complete system. The focus will be on testing the application from its primary entry point: the CLI.

*   **CLI End-to-End Test (`test_pipeline.py`):** This will be the most important test of Cycle 1.
    *   We will use the `CliRunner` utility provided by `Typer` to invoke our CLI command programmatically within a test function.
    *   The test will be executed within an isolated temporary directory created by `pytest`'s `tmp_path` fixture. This ensures the test does not leave artifacts on the filesystem.
    *   A minimal, valid `config.yaml` file will be created within this temporary directory (e.g., generate 5 structures of a simple Cu-Au alloy with 10 MD steps, and sample 8 final structures).
    *   The `CliRunner` will execute the command: `mlip_autopipec --config config.yaml`.
    *   After the command execution, a series of assertions will validate the outcome:
        1.  The command must exit with a success code of 0.
        2.  The specified ASE database file (`final_dataset.db`) must exist in the temporary directory.
        3.  The test will connect to this database using `ase.db.connect` and assert that it contains exactly 8 rows, matching the `num_samples` in the configuration.
    *   This single test validates the entire application workflow, from configuration parsing by Hydra and Pydantic, through all pipeline stages, to the final database creation, ensuring all the pieces are correctly and robustly integrated.
