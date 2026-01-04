# CYCLE01 Specification: Core CLI Pipeline

## 1. Summary

This document provides the detailed technical specification for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 1 is to build and deliver a fully functional, command-line-driven pipeline capable of generating a foundational dataset for MLIP training. This cycle focuses on implementing the core components of the system: configuration management, initial structure generation, a basic exploration engine, a simple sampling method, and persistent storage. The scope is intentionally focused on a single, well-defined workflow: generating a dataset for a binary alloy system (e.g., FePt). This allows us to establish a robust and extensible software architecture that will serve as the backbone for all future enhancements.

By the end of this cycle, a user will be able to define a simple alloy in a YAML configuration file, execute the pipeline via a single command-line instruction, and receive a well-structured ASE database as output. This initial version will not include the advanced hybrid MD/MC simulations or the sophisticated Farthest Point Sampling; those features are deferred to Cycle 2. Instead, Cycle 1 implements a standard NVT molecular dynamics simulation and a straightforward random sampling technique. This approach ensures that the fundamental data flow and component integrations are correctly implemented and tested before adding more complex scientific logic. The engineering focus is on creating a modular, testable, and reliable system. Key architectural patterns, such as the use of Abstract Base Classes for generators and the factory pattern for their instantiation, will be implemented from the outset. Emphasis will also be placed on robust error handling and clear dependency management using `uv` and `pyproject.toml`. The successful completion of this cycle will yield a valuable tool for MLIP practitioners and validate the core architectural design of the MLIP-AutoPipe framework.

## 2. System Architecture

The architecture for Cycle 1 is a direct implementation of the foundational components outlined in the main `SYSTEM_ARCHITECTURE.md` document. The focus is on creating the essential file and class structures that form the pipeline's backbone.

**File Structure:**
The following tree shows the files to be created or modified in this cycle. Bold entries indicate new files to be created.

```
.
├── pyproject.toml
└── src/
    └── mlip_autopipec/
        ├── __init__.py
        ├── **cli.py**                  # Main entry point for the Command-Line Interface
        ├── pipeline/
        │   ├── __init__.py
        │   └── **runner.py**             # Contains the PipelineRunner orchestrator class
        ├── generators/
        │   ├── __init__.py
        │   ├── **base.py**               # Abstract Base Class for all generators
        │   ├── **alloy.py**              # Generator for alloy structures
        │   └── **factory.py**            # Factory to select the appropriate generator
        ├── explorers/
        │   ├── __init__.py
        │   └── **md_engine.py**          # Core MD simulation engine (basic version)
        ├── sampling/
        │   ├── __init__.py
        │   ├── **base.py**               # Abstract Base Class for samplers
        │   └── **random_sampler.py**     # Simple random sampling implementation
        ├── storage/
        │   ├── __init__.py
        │   └── **ase_db_writer.py**      # Module for writing to ASE DB
        └── common/
            ├── __init__.py
            ├── **atoms_validator.py**    # Physics-based validation logic
            └── **exceptions.py**         # Custom exception classes
```

**Code Blueprints:**

*   **`cli.py`**: This file will use the `typer` library to create the main command-line interface. A single command, `run`, will be implemented. This command will accept the path to a Hydra configuration file. Its primary responsibility is to parse the configuration and pass it to the `PipelineRunner`.
*   **`pipeline/runner.py`**:
    *   `PipelineRunner`: A class initialized with the Hydra config object. It will have a public method `run()`. Inside `run()`, it will orchestrate the entire workflow:
        1.  Instantiate the appropriate generator using `GeneratorFactory`.
        2.  Call the generator to create initial structures.
        3.  For each structure, instantiate and run the `MDEngine`.
        4.  Instantiate the `RandomSampler` and call it on the collected trajectory files.
        5.  Call the `save_to_db` function from `ase_db_writer` to store the final structures.
*   **`generators/base.py`**:
    *   `BaseStructureGenerator` (ABC): An abstract class defining the interface `generate(self, **kwargs) -> list[ase.Atoms]`. It will also contain concrete methods for common validation checks, leveraging `atoms_validator.py`.
*   **`generators/alloy.py`**:
    *   `AlloyGenerator(BaseStructureGenerator)`: A concrete implementation for generating alloy structures. It will take parameters like `elements`, `composition`, `crystal_structure`, and `supercell_size`. It will use `ase.build` functions to create the initial crystal lattice and then randomly substitute atoms to achieve the target composition.
*   **`generators/factory.py`**:
    *   `GeneratorFactory`: A class with a single static method `get_generator(config: DictConfig) -> BaseStructureGenerator`. This method will inspect the configuration (e.g., `config.system.type`) and return an instance of the corresponding generator (`AlloyGenerator` in this cycle).
*   **`explorers/md_engine.py`**:
    *   `MDEngine`: A class that performs a basic molecular dynamics simulation. Its `run(atoms: ase.Atoms)` method will:
        1.  Attach an MLIP calculator (e.g., ASE's EMT for testing, or a pre-trained MACE model).
        2.  Set up an NVT simulation using an ASE dynamics object (e.g., `Langevin`).
        3.  Run the dynamics for a specified number of steps.
        4.  Save the trajectory to a file.
*   **`sampling/random_sampler.py`**:
    *   `RandomSampler(BaseSampler)`: An implementation that reads a trajectory file and randomly selects a specified number of frames to be included in the final dataset.
*   **`storage/ase_db_writer.py`**:
    *   A single function `save_to_db(atoms_list: list[ase.Atoms], db_path: str)` that connects to an ASE database at `db_path` and writes each `Atoms` object to it, along with essential metadata like energy and forces.

## 3. Design Architecture

The design for Cycle 1 is centered around establishing clear data contracts and interactions between components using Pydantic models for validation and configuration. This ensures a robust, type-safe foundation.

*   **Configuration Schema (Conceptual Pydantic Models):**
    Although Hydra manages the configuration, we will internally validate the structure using Pydantic models to prevent runtime errors.

    ```python
    # In a new file, e.g., src/mlip_autopipec/config_models.py

    from pydantic import BaseModel, Field
    from typing import List, Dict

    class SystemConfig(BaseModel):
        type: str = "Alloy"
        elements: List[str]
        composition: Dict[str, float]
        crystal_structure: str
        num_initial_structures: int = Field(gt=0)

    class ExplorationConfig(BaseModel):
        temperature_k: float = Field(gt=0)
        num_steps: int = Field(gt=0)
        calculator: str # e.g., "MACE" or "EMT"

    class SamplingConfig(BaseModel):
        method: str = "Random"
        num_samples: int = Field(gt=0)

    class MainConfig(BaseModel):
        system: SystemConfig
        exploration: ExplorationConfig
        sampling: SamplingConfig
        db_path: str
    ```
    This `MainConfig` model will be used by the `PipelineRunner` to validate and access configuration data in a structured, predictable way. It serves as the primary data contract for the entire pipeline.

*   **Component Interactions and Data Flow:**
    1.  **Producer (`cli.py`):** The CLI's primary role is to produce a validated `MainConfig` object from the raw Hydra `DictConfig`.
    2.  **Consumer (`PipelineRunner`):** Consumes the `MainConfig`. It acts as a producer of `ase.Atoms` objects for the next stages.
    3.  **Producer (`AlloyGenerator`):** Consumes the `SystemConfig` portion of the main config. It produces a list of initial `ase.Atoms` objects. These objects are a critical data structure, containing atomic numbers, positions, and cell information.
    4.  **Consumer/Producer (`MDEngine`):** Consumes a single `ase.Atoms` object. It produces a trajectory file (e.g., `.xyz`), which is a textual representation of atomic coordinates over time.
    5.  **Consumer/Producer (`RandomSampler`):** Consumes the path to a trajectory file. It parses this file and produces a curated list of `ase.Atoms` objects, representing the final dataset.
    6.  **Consumer (`ase_db_writer`):** Consumes the final list of `ase.Atoms` objects and writes them to the database, producing no further data.

*   **Invariants and Constraints:**
    *   **Atom Validation:** A key invariant, enforced by `atoms_validator.py`, is that no two atoms can be closer than a specified threshold. The `AlloyGenerator` must ensure all its output structures satisfy this constraint.
    *   **Configuration Validity:** The `PipelineRunner` will immediately validate the entire configuration using the Pydantic model. If validation fails, the program will exit with a clear error message before any computation begins. This prevents difficult-to-debug errors deep within the pipeline.
    *   **Data Immutability (Conceptual):** While Python objects are mutable, the design encourages treating data as immutable between stages. The `MDEngine` receives an `Atoms` object and produces a new file; it does not modify the input object in place. This separation of concerns simplifies the logic and improves testability.

## 4. Implementation Approach

The implementation will proceed in a logical, bottom-up fashion, starting with the core components and culminating in the CLI that ties them together.

1.  **Setup `pyproject.toml`:** Define the project metadata and add initial dependencies: `typer`, `hydra-core`, `ase`, `pydantic`. Run `uv sync` to create the virtual environment.
2.  **Implement Core Data Structures:** Create the Pydantic configuration models in `src/mlip_autopipec/config_models.py`. This defines the "shape" of the data that will flow through the system.
3.  **Implement `common` Utilities:**
    *   Create `common/exceptions.py` with a custom `PhysicsViolationError`.
    *   Implement the `overlap_check` function in `common/atoms_validator.py`.
4.  **Implement `storage` Module:** Create the `ase_db_writer.py` file with the `save_to_db` function. This is a simple, isolated component that is easy to implement and test first.
5.  **Implement `generators` Module:**
    *   Define the `BaseStructureGenerator` ABC in `base.py`.
    *   Implement the `AlloyGenerator` in `alloy.py`. It will use `ase.build.bulk` to create a pristine supercell and then use `numpy.random.choice` to replace elements according to the specified composition. It will call the `overlap_check` validator before returning the structures.
    *   Implement the `GeneratorFactory` in `factory.py`.
6.  **Implement `sampling` Module:**
    *   Define the `BaseSampler` ABC in `base.py`.
    *   Implement the `RandomSampler` in `random_sampler.py`. It will use `ase.io.read` with the `index=':'` slice to read all frames from a trajectory, then use `random.sample` to select the desired number.
7.  **Implement `explorers` Module:**
    *   Implement the `MDEngine` class in `md_engine.py`. For this cycle, the calculator can be hard-coded to ASE's `EMT` for simplicity and speed. It will use `ase.md.velocitydistribution.MaxwellBoltzmannDistribution` to set initial velocities and `ase.md.langevin.Langevin` for the dynamics. An `ase.io.trajectory.Trajectory` object can be used to write the output.
8.  **Implement `pipeline` Orchestrator:**
    *   Implement the `PipelineRunner` class. This is where all the previously implemented components are brought together. The `run` method will execute the full sequence of operations. It should include logging (`print` statements are acceptable for this cycle) to indicate the current stage of the process.
9.  **Implement `cli.py`:**
    *   Create a `typer.Typer` app.
    *   Define a `run` command that takes the config path.
    *   Use Hydra's `initialize` and `compose` functions to load the YAML configuration.
    *   Instantiate and run the `PipelineRunner`.

This step-by-step process ensures that each component is built on a solid foundation, with dependencies being implemented before the components that rely on them.

## 5. Test Strategy

Testing in Cycle 1 is crucial for validating the core architecture and ensuring each component functions correctly in isolation and as part of the whole.

**Unit Testing Approach (Min 300 words):**
Unit tests will be written using `pytest`. The primary goal is to test each module's logic independently, using mocks to isolate it from its dependencies.
*   **`atoms_validator`**: We will create test `ase.Atoms` objects, some with overlapping atoms and some without, and assert that the `overlap_check` function returns the expected boolean values.
*   **`AlloyGenerator`**: We will test this by providing a simple configuration (e.g., a 2-atom SiC cell). We will assert that the output is a list of `ase.Atoms` objects, that the number of atoms is correct, and that the chemical symbols (`at.get_chemical_symbols()`) match the requested 1:1 composition. We don't need to check the exact positions, as they have a random component, but we can check that they satisfy the `overlap_check`.
*   **`RandomSampler`**: We will create a dummy trajectory file with a known number of frames (e.g., 10). We will run the sampler and assert that the returned list of `ase.Atoms` objects has the correct length (as specified in the sampling config) and that the selected atoms are indeed a subset of the original trajectory.
*   **`MDEngine`**: Testing the MD engine is more complex. We will not test the physics of the simulation itself, as that is the job of the ASE library. Instead, we will test the *orchestration* of the simulation. We will provide a simple `Atoms` object and run the engine for a very small number of steps (e.g., 5). We will use `mocker.patch` to mock the ASE dynamics object's `run()` method to avoid actual computation. The test will assert that a trajectory file is created and that it is not empty.
*   **`PipelineRunner`**: We will test the runner by mocking all the components it calls (`GeneratorFactory`, `MDEngine`, `RandomSampler`, `save_to_db`). The test will provide a sample configuration and assert that each mocked component's primary method (`generate`, `run`, `sample`, `save_to_db`) is called exactly once and in the correct sequence. This verifies the orchestration logic without performing any actual work.

**Integration Testing Approach (Min 300 words):**
A single, comprehensive integration test will validate the entire pipeline from the command line. This test ensures that the "plumbing" between the components is correct and that data flows through the system as expected.
*   The test will use `typer.testing.CliRunner`. The `runner.invoke()` method will be used to execute the `run` command defined in `cli.py`.
*   A temporary directory will be created for the test using `pytest`'s `tmp_path` fixture. A minimal Hydra configuration file will be created in this directory.
*   This configuration will be designed for speed and simplicity:
    *   `system.num_initial_structures`: 1
    *   The system will be a very small unit cell (e.g., 2 atoms of Si).
    *   `exploration.num_steps`: 5 (just enough to generate a short trajectory).
    *   `exploration.calculator`: `EMT` (ASE's fast built-in potential).
    *   `sampling.num_samples`: 2
    *   `db_path`: A path within the temporary directory.
*   The test will invoke the CLI and assert that it exits with a success code (`result.exit_code == 0`).
*   After the command completes, the test will inspect the temporary directory. It will assert that:
    1.  An intermediate trajectory file was created.
    2.  The final database file (e.g., `structures.db`) exists at the specified path.
*   Finally, the test will connect to the output database using `ase.db.connect()` and assert that it contains the correct number of entries (2, as specified in `num_samples`). This confirms that the data has successfully passed through every stage of the pipeline, from generation to storage.
