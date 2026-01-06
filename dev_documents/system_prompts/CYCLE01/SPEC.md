# Specification: CYCLE 01 - Core Engine and Structure Generation

## 1. Summary

This document outlines the detailed specifications for the first development cycle of the MLIP-AutoPipe project. The primary objective of Cycle 01 is to establish the foundational architecture and implement the core functionalities necessary for a minimal viable product (MVP). This cycle focuses on creating a robust command-line application that can generate physically valid initial atomic structures for two key material classes—alloys and ionic crystals—and subsequently perform basic molecular dynamics (MD) simulations to explore their potential energy surfaces. The core deliverable of this cycle is a functional, albeit simplified, end-to-end pipeline that can be executed from the command line, taking a user-defined configuration file and producing a database of explored atomic configurations.

The scope of this cycle is intentionally focused on building the essential components that will serve as the backbone for more advanced features in subsequent cycles. This includes setting up a clean, modular project structure, implementing a powerful and type-safe configuration system using Pydantic, and creating a flexible structure generation framework based on an abstract factory pattern. The generation process will be equipped with critical physics-based validation logic from the outset, ensuring that all generated structures are physically plausible by preventing atomic overlaps and adhering to geometric constraints.

Furthermore, this cycle will introduce a parallelized exploration engine. While it will be a simplified version of the final engine, it will leverage Python's `ProcessPoolExecutor` to run multiple MD simulations concurrently, demonstrating the system's capability for high-throughput computation. The simulations will be driven by pre-trained MLIP models like MACE, and the resulting trajectories will be saved. A simple, robust wrapper for the Atomic Simulation Environment (ASE) database will also be implemented, providing a standardized mechanism for storing and retrieving atomic structures. By the end of this cycle, the project will have a solid, testable foundation, including a command-line interface, a configuration parser, two distinct structure generators, a parallel simulation engine, and a database storage solution. This MVP will be a significant milestone, providing a tangible tool that can already be used for basic dataset generation tasks, paving the way for the integration of more sophisticated exploration and sampling techniques in Cycle 02. The emphasis is on building high-quality, well-tested, and modular components that ensure the long-term maintainability and extensibility of the MLIP-AutoPipe framework.

## 2. System Architecture

The architecture for Cycle 01 is designed to be a modular and scalable foundation for the entire MLIP-AutoPipe application. The system is orchestrated by a central `PipelineRunner` that is controlled via a command-line interface. The core logic is broken down into distinct services for configuration, generation, exploration, and storage, ensuring a clean separation of concerns.

**File Structure (Cycle 01):**

The following ASCII tree shows the files to be created or modified in this cycle. Files marked in **bold** are the primary targets for creation and implementation.

```
.
├── pyproject.toml
└── src/
    └── mlip_autopipec/
        ├── __init__.py
        ├── **cli.py**                # Main CLI entry point (Typer)
        ├── config/
        │   ├── __init__.py
        │   └── **models.py**         # Pydantic models for configuration
        ├── generators/
        │   ├── __init__.py
        │   ├── **base.py**           # Abstract BaseGenerator class
        │   ├── **alloy.py**          # Alloy structure generator
        │   └── **ionic.py**          # Ionic structure generator
        ├── explorers/
        │   ├── __init__.py
        │   └── **engine.py**         # Basic MD exploration engine
        ├── storage/
        │   ├── __init__.py
        │   └── **database.py**       # ASE Database wrapper
        └── **pipeline.py**           # Core PipelineRunner orchestrator
```

**Code Blueprints:**

This section provides the detailed blueprints for the key components to be implemented in Cycle 01.

**`src/mlip_autopipec/cli.py`**
This file will serve as the entry point for the user. It will use the `typer` library to create a simple, clean command-line interface.

```python
import typer
from pathlib import Path
import yaml

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.pipeline import PipelineRunner

app = typer.Typer()

@app.command()
def run(
    config_path: Path = typer.Option(..., "--config", "-c", help="Path to the configuration YAML file.", exists=True)
):
    """
    Runs the MLIP-AutoPipe pipeline from a configuration file.
    """
    try:
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)

        config = FullConfig(**config_dict)

        runner = PipelineRunner(config)
        runner.run()

    except Exception as e:
        typer.echo(f"An error occurred: {e}")
        raise typer.Exit(code=1)

if __name__ == "__main__":
    app()
```

**`src/mlip_autopipec/pipeline.py`**
The `PipelineRunner` is the central orchestrator. It will initialize all necessary services and execute the pipeline stages in the correct order.

```python
from .config.models import FullConfig
from .generators.base import BaseStructureGenerator
from .explorers.engine import ExplorerEngine
from .storage.database import AseDBWrapper

class PipelineRunner:
    def __init__(self, config: FullConfig):
        self.config = config
        self.db_wrapper = AseDBWrapper(path=self.config.db_path)

    def _get_generator(self) -> BaseStructureGenerator:
        # Factory logic will go here
        if self.config.system.type == "alloy":
            from .generators.alloy import AlloyGenerator
            return AlloyGenerator(self.config.system)
        elif self.config.system.type == "ionic":
            from .generators.ionic import IonicGenerator
            return IonicGenerator(self.config.system)
        else:
            raise NotImplementedError(f"System type '{self.config.system.type}' not supported.")

    def run(self):
        # Stage 1: Generation
        generator = self._get_generator()
        initial_structures = generator.generate()
        self.db_wrapper.write_atoms_list(initial_structures, "initial")

        # Stage 2: Exploration
        explorer = ExplorerEngine(self.config.exploration)
        explored_trajectories = explorer.run_parallel(initial_structures)

        # In Cycle 01, we just store the trajectories directly
        for i, trajectory in enumerate(explored_trajectories):
            self.db_wrapper.write_atoms_list(trajectory, f"explored_traj_{i}")
```

**`src/mlip_autopipec/generators/base.py`**
This file defines the abstract interface for all structure generators, ensuring that they all adhere to a common contract and share validation logic.

```python
from abc import ABC, abstractmethod
from typing import List
from ase import Atoms
import numpy as np

from mlip_autopipec.config.models import SystemConfig

class BaseStructureGenerator(ABC):
    def __init__(self, config: SystemConfig):
        self.config = config

    @abstractmethod
    def _create_structures(self) -> List[Atoms]:
        """Subclasses must implement this method to create the base structures."""
        pass

    def generate(self) -> List[Atoms]:
        """Generates and validates a list of atomic structures."""
        structures = self._create_structures()
        validated_structures = []
        for atoms in structures:
            if self._is_valid(atoms):
                validated_structures.append(atoms)
        return validated_structures

    def _is_valid(self, atoms: Atoms) -> bool:
        """Performs physics-based validation checks."""
        return self._passes_overlap_check(atoms)

    def _passes_overlap_check(self, atoms: Atoms, min_dist_factor: float = 0.7) -> bool:
        """Checks if any two atoms are too close to each other."""
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        # A simple check based on covalent radii could be implemented here
        # For now, we'll use a placeholder.
        return np.min(distances) > min_dist_factor
```
This detailed architectural blueprint provides a clear and actionable plan for the development team. It establishes a modular framework that is easy to understand, test, and extend. The use of an abstract base class for generators and a central pipeline orchestrator promotes a clean, decoupled design that will be critical for managing the complexity of the project as it grows in subsequent cycles. The immediate inclusion of a Pydantic-based configuration system ensures that from day one, all inputs are validated, leading to a more robust and user-friendly application.

## 3. Design Architecture

The design architecture for Cycle 01 is centered around creating a robust, type-safe, and configuration-driven system. The core of this design is the use of Pydantic for defining all configuration schemas. This approach, often referred to as "schema-first," ensures that all data flowing into the application is validated against a clear, explicit, and self-documenting contract. This prevents a wide range of runtime errors and makes the system easier to use and debug.

**Pydantic-based Schema Design:**

The entire configuration of the pipeline will be encapsulated within a single root Pydantic model, `FullConfig`. This model will be composed of several nested models, each representing a specific domain of configuration (e.g., system definition, exploration parameters).

File: `src/mlip_autopipec/config/models.py`

```python
from pydantic import BaseModel, Field
from typing import List, Dict, Literal

class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""
    type: Literal["alloy", "ionic"] = Field(..., description="The type of material system.")
    elements: List[str] = Field(..., min_length=1, description="List of element symbols.")
    composition: Dict[str, float] = Field(..., description="Composition of elements (fractions must sum to 1).")
    num_structures: int = Field(10, gt=0, description="Number of initial structures to generate.")
    # Add other relevant fields like lattice parameters, etc.

class ExplorationConfig(BaseModel):
    """Configuration for the exploration phase."""
    temperature_k: float = Field(300.0, gt=0, description="MD simulation temperature in Kelvin.")
    pressure_gpa: float = Field(0.0, ge=0, description="MD simulation pressure in GPa.")
    time_step_fs: float = Field(1.0, gt=0, description="Time step for MD integration in femtoseconds.")
    num_steps: int = Field(1000, gt=0, description="Total number of MD steps.")
    mlip_model: str = Field("MACE", description="Name of the MLIP model to use.")
    max_workers: int = Field(4, gt=0, description="Number of parallel processes for exploration.")

class FullConfig(BaseModel):
    """Root configuration model for the MLIP-AutoPipe pipeline."""
    db_path: str = Field("mlip_data.db", description="Path to the output ASE database.")
    system: SystemConfig
    exploration: ExplorationConfig
```

**Domain Concepts and Invariants:**
*   **`SystemConfig`**: This model represents the definition of the material to be simulated.
    *   **Invariants**: The sum of the values in the `composition` dictionary must be equal to 1.0. This can be enforced with a Pydantic `model_validator`. The `elements` list must contain only valid chemical symbols.
    *   **Producers**: This data is produced by the user, written in a YAML file.
    *   **Consumers**: The `PipelineRunner` consumes this to select the correct generator, and the specific `BaseStructureGenerator` implementation uses it to create the atomic structures.
*   **`ExplorationConfig`**: This model defines the parameters for the MD simulation.
    *   **Invariants**: All numerical values must be physically reasonable (e.g., temperature and time step must be greater than zero). Pydantic's `Field` constraints (`gt`, `ge`) enforce these rules at the time of parsing.
    *   **Producers**: User-defined in the YAML configuration.
    *   **Consumers**: The `ExplorerEngine` consumes this model to configure the parallel MD runs.

**Extensibility and Versioning:**
This Pydantic-based approach is highly extensible. Adding a new parameter is as simple as adding a new field to the model, with a default value for backward compatibility. If a future version of the pipeline requires a different configuration structure, a new version of the `FullConfig` model can be created, and the `PipelineRunner` can be updated to handle multiple configuration versions, ensuring a smooth transition for users. This schema-first design provides a solid foundation for building a complex application in a structured and maintainable way.

## 4. Implementation Approach

The implementation of Cycle 01 will proceed in a logical, step-by-step manner, starting with the foundational components and progressively building up to the fully integrated MVP. This approach ensures that each component can be tested independently before being integrated into the main pipeline.

**Step 1: Project Scaffolding and Dependency Management**
The first action is to set up the project structure as defined in the System Architecture section. The `pyproject.toml` file will be created and populated with the initial set of dependencies. These will include `typer` for the CLI, `pydantic` for configuration, `pyyaml` for parsing the config file, `ase` for atomic data structures, and `pymatgen` for crystallographic utilities. Using `uv` to manage the virtual environment and dependencies is recommended.

**Step 2: Implement the Configuration System**
The Pydantic models outlined in the Design Architecture (`SystemConfig`, `ExplorationConfig`, `FullConfig`) will be implemented in `src/mlip_autopipec/config/models.py`. This step is critical as it defines the data contract for the entire application. Unit tests will be written in parallel to verify that the models correctly parse valid YAML files and reject invalid ones.

**Step 3: Develop the CLI and PipelineRunner**
The `cli.py` file will be created with the main `run` command using Typer. This command will be responsible for loading the YAML file, instantiating the `FullConfig` model, and passing it to the `PipelineRunner`. The initial implementation of `PipelineRunner` in `pipeline.py` will be a skeleton class that takes the config object. This creates the main application entry point and the central orchestrator.

**Step 4: Implement the Structure Generators**
The abstract base class `BaseStructureGenerator` will be created in `generators/base.py`, defining the common interface and the shared validation logic (e.g., `_passes_overlap_check`). Subsequently, the two concrete implementations, `AlloyGenerator` and `IonicGenerator`, will be developed.
*   **`AlloyGenerator` (`generators/alloy.py`):** This will use `pymatgen` or ASE's tools to create a primitive cell and then build a supercell. It will then randomly replace atomic species according to the specified composition.
*   **`IonicGenerator` (`generators/ionic.py`):** This will be more complex, potentially using `pymatgen`'s `Structure.from_spacegroup` or similar methods to create charge-neutral structures based on known crystal prototypes.

**Step 5: Implement the Database Wrapper**
The `AseDBWrapper` class will be created in `storage/database.py`. This class will encapsulate all interactions with the ASE database. It will have a simple API with methods like `__init__(path)`, `write_atoms_list(atoms_list, group_name)`, and potentially `read_group(group_name)`. This isolates database logic from the main pipeline.

**Step 6: Integrate Generation and Storage**
The `PipelineRunner` will be updated to include the factory logic (`_get_generator`) that selects the appropriate generator based on the configuration. The `run` method will be implemented to execute the generation stage, call the generator to create the initial structures, and then use the `AseDBWrapper` to write these structures to the database.

**Step 7: Implement the Basic Exploration Engine**
The `ExplorerEngine` will be implemented in `explorers/engine.py`. Its primary responsibility will be to manage the parallel execution of MD simulations. It will have a `run_parallel` method that takes the list of initial structures. This method will use `ProcessPoolExecutor` to map a worker function to each structure. The worker function will take an `Atoms` object, attach an ASE calculator (e.g., a simple Lennard-Jones potential for initial testing, or a pre-trained MLIP), run a short MD simulation, and return the trajectory.

**Step 8: Final Integration**
The `PipelineRunner`'s `run` method will be extended to call the `ExplorerEngine` after the generation stage. The trajectories returned by the explorer will be written to the database using the `AseDBWrapper`. At this point, the full MVP pipeline for Cycle 01 will be complete and ready for end-to-end integration testing.

## 5. Test Strategy

The test strategy for Cycle 01 is designed to build confidence in the core components and their interactions. It will consist of a comprehensive suite of unit tests to validate individual modules in isolation and a set of integration tests to ensure they work together correctly.

**Unit Testing Approach (Min 300 words):**
Unit tests are the first line of defense against bugs and are crucial for ensuring the correctness of each component's logic. All unit tests will be placed in a `tests/unit` directory that mirrors the `src` directory structure.
*   **Configuration (`tests/unit/config/test_models.py`):** These tests will focus on the Pydantic models. We will create several sample YAML files, both valid and invalid. The tests will assert that `FullConfig.parse_file()` correctly loads the valid files and that `pydantic.ValidationError` is raised for invalid files. Specific tests will target edge cases, such as compositions that do not sum to 1.0 or negative values for physical parameters like temperature. This ensures the robustness of our input validation.
*   **Generators (`tests/unit/generators/test_alloy.py`, `test_ionic.py`):** Testing the generators is critical for ensuring physical validity. We will test the `AlloyGenerator` by providing it with a simple binary system configuration. The test will then assert that the generated `Atoms` objects have the correct number of atoms and that the ratio of the elements matches the requested composition. A key test will be to invoke the `_passes_overlap_check` method directly to ensure it correctly identifies and rejects structures where atoms are too close. Similarly, for the `IonicGenerator`, tests will verify charge neutrality in the final structures.
*   **Database (`tests/unit/storage/test_database.py`):** The `AseDBWrapper` will be tested using an in-memory SQLite database to avoid filesystem side effects. The test will create a list of sample `ase.Atoms` objects, use the wrapper to write them to the database, and then read them back. The core assertion will be that the `Atoms` objects read back are identical to the original ones, ensuring the integrity of our data persistence layer.

**Integration Testing Approach (Min 300 words):**
Integration tests will verify the interactions and data flow between the different components of the pipeline. These tests will be located in the `tests/integration` directory and will focus on complete or partial workflows.
*   **Generation Pipeline (`tests/integration/test_generation_pipeline.py`):** The first major integration test will verify the generation stage from end to end. It will involve creating a temporary YAML configuration file on the fly. The test will then invoke the `cli.py` runner with this config file. The primary assertion will be to check that a database file is created and that this database contains the correct number of initial structures as specified in the configuration. This test ensures that the CLI, configuration parsing, generator factory, and database writer all work together seamlessly.
*   **Full MVP Pipeline (`tests/integration/test_full_pipeline_cycle01.py`):** The cornerstone test for Cycle 01 will be an end-to-end test of the MVP. This test will use a minimal but complete configuration file that specifies the generation of one or two simple structures and a very short exploration run (e.g., 10 steps). The test will execute the full `run` command. After the run completes, the test will open the resulting ASE database and perform several assertions:
    1.  Verify that the "initial" group of structures exists.
    2.  Verify that the "explored_traj_X" groups exist.
    3.  Check that the number of atoms in the explored structures is consistent with the initial structures.
    4.  Optionally, check if the temperature of the system during the simulation was close to the requested temperature. This test validates the entire workflow, from user input to final output, providing high confidence in the integrated system.
