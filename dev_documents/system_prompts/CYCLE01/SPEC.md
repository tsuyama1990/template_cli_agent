# CYCLE01/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 1 of the MLIP-AutoPipe project. The primary objective of this cycle is to deliver a robust, command-line-driven application capable of generating physically valid seed structures for various material systems and storing them in a standardized database format. This foundational cycle focuses on establishing the core architecture, data structures, and workflow of the application. The scope of this cycle is strictly limited to the **Generation** and **Storage** stages of the pipeline. The more complex Exploration and Sampling stages are deferred to Cycle 2. At the end of this cycle, we will have a functional command-line interface (CLI) that accepts a user-defined configuration file, processes it through a pipeline, and outputs an ASE (Atomic Simulation Environment) database containing a set of high-quality initial atomic configurations.

The core components to be developed in this cycle include the Pydantic-based configuration schemas for defining material systems, the main `PipelineRunner` class that orchestrates the workflow, a factory for dynamically selecting the appropriate structure generator, and a comprehensive suite of `Generator` classes for different material types (e.g., alloys, ionic compounds). A significant emphasis is placed on ensuring the physical validity of all generated structures. Each generator will incorporate a series of rigorous validation checks, including minimum interatomic distance enforcement, automatic supercell expansion to satisfy periodic boundary conditions, and charge neutrality for ionic systems. The final output will be managed by an `AseDBWrapper` class, which will handle the serialization of generated structures and their associated metadata into a portable and widely-compatible SQLite database file. This initial deliverable will provide immediate value to researchers by automating the tedious and error-prone process of creating initial structures for computational simulations, while also laying a solid, well-tested groundwork for the advanced features planned in the subsequent cycle.

## 2. System Architecture

The architecture for Cycle 1 is a simplified version of the final system, focusing exclusively on the initial and final steps of the pipeline. The system will be implemented within the `src/mlip_autopipec` package.

**File Structure (Cycle 1 Focus):**

The files and directories to be created or modified in this cycle are marked in bold.

```
src/mlip_autopipec/
├── __init__.py
├── **cli.py**                  # CLI entry point using click
├── config/
│   ├── __init__.py
│   └── **schemas.py**          # Pydantic models for configuration
├── pipeline/
│   ├── __init__.py
│   └── **runner.py**           # PipelineRunner for Gen -> Store workflow
├── generators/
│   ├── __init__.py
│   ├── **base.py**             # Abstract BaseGenerator class
│   ├── **factory.py**          # GeneratorFactory implementation
│   ├── **alloy.py**            # AlloyGenerator implementation
│   ├── **ionic.py**            # IonicGenerator implementation
│   └── **...**                 # Other generator stubs
├── storage/
│   ├── __init__.py
│   └── **db_wrapper.py**       # ASE Database wrapper
└── utils/
    ├── __init__.py
    └── **physics.py**          # Physics validation functions
```

**Code Blueprints:**

-   **`cli.py`**:
    ```python
    import click
    from mlip_autopipec.config.schemas import FullConfig
    from mlip_autopipec.pipeline.runner import PipelineRunner
    # Hydra or a similar library will be used to load config from yaml

    @click.command()
    @click.argument('config_path', type=click.Path(exists=True))
    def main(config_path):
        """Runs the MLIP-AutoPipe structure generation pipeline."""
        # config = load_config_from_yaml(config_path, FullConfig) # Pseudocode
        config = FullConfig(...) # Load and validate using Pydantic
        runner = PipelineRunner(config)
        runner.run()
        click.echo("Pipeline completed successfully.")
    ```

-   **`pipeline/runner.py`**:
    ```python
    from mlip_autopipec.config.schemas import FullConfig
    from mlip_autopipec.generators.factory import GeneratorFactory
    from mlip_autopipec.storage.db_wrapper import AseDBWrapper

    class PipelineRunner:
        def __init__(self, config: FullConfig):
            self.config = config

        def run(self):
            # Stage 1: Generation
            generator = GeneratorFactory.create_generator(self.config.system)
            seed_structures = generator.generate()

            # Stage 4: Storage
            db_wrapper = AseDBWrapper(filepath=self.config.storage.db_path)
            db_wrapper.write_structures(seed_structures)
    ```

-   **`generators/base.py`**:
    ```python
    from abc import ABC, abstractmethod
    from typing import List
    from ase import Atoms
    from mlip_autopipec.config.schemas import SystemConfig
    from mlip_autopipec.utils.physics import validate_structure

    class BaseGenerator(ABC):
        def __init__(self, config: SystemConfig):
            self.config = config

        def generate(self) -> List[Atoms]:
            structures = self._generate_initial()
            validated_structures = []
            for atoms in structures:
                # Apply common transformations and validations
                # e.g., apply_strain, apply_rattle
                validated = self._validate(atoms)
                if validated:
                    validated_structures.append(atoms)
            return validated_structures

        @abstractmethod
        def _generate_initial(self) -> List[Atoms]:
            """Subclass-specific generation logic."""
            pass

        def _validate(self, atoms: Atoms) -> bool:
            """Perform rigorous physics-based validation."""
            return validate_structure(atoms, min_dist=self.config.validation.min_dist)

    ```

-   **`storage/db_wrapper.py`**:
    ```python
    from typing import List
    from ase import Atoms
    from ase.db import connect

    class AseDBWrapper:
        def __init__(self, filepath: str):
            self.filepath = filepath

        def write_structures(self, structures: List[Atoms]):
            with connect(self.filepath) as db:
                for atoms in structures:
                    # Extract metadata to be stored alongside the structure
                    metadata = {"energy": atoms.get_potential_energy()} # Example
                    db.write(atoms, key_value_pairs=metadata)
    ```

This architecture ensures a clear separation of responsibilities. The CLI handles user interaction, the `PipelineRunner` orchestrates the workflow, the `Generators` encapsulate the complex logic of structure creation, and the `AseDBWrapper` abstracts away the database interactions.

## 3. Design Architecture

This cycle establishes the foundational data schemas using Pydantic, ensuring that the entire application is built on a type-safe and validated data model. This schema-first approach is critical for the long-term stability and maintainability of the project.

**Configuration Schemas (`config/schemas.py`):**

The core of the design is a single, comprehensive `FullConfig` model that nests all other configuration objects. This provides a single source of truth for all parameters.

-   **`ValidationConfig`**: Defines parameters for physical validation.
    -   `min_dist` (float): The minimum allowed distance between any two atoms. Key constraint: must be a positive value.
    -   `min_cell_size` (float): The minimum supercell size, often tied to the MLIP's cutoff radius.

-   **`SystemConfig`**: The most critical schema, defining the material to be generated. It will likely be a discriminated union to handle different system types.
    -   `system_type` (Literal["alloy", "ionic", ...]): The discriminant field.
    -   `elements` (List[str]): List of chemical symbols.
    -   `num_structures` (int): Number of seed structures to generate. Constraint: Must be greater than 0.

-   **`AlloySystemConfig(SystemConfig)`**:
    -   `composition` (Dict[str, float]): A dictionary mapping elements to their fractional composition. Constraint: Values must sum to 1.0.
    -   `lattice_type` (str): E.g., 'fcc', 'bcc'.

-   **`IonicSystemConfig(SystemConfig)`**:
    -   `stoichiometry` (Dict[str, int]): A dictionary mapping elements to their integer counts in the formula unit.
    -   `charge_states` (Dict[str, float]): Expected charge state for each element. Constraint: The overall structure must be charge-neutral.

-   **`StorageConfig`**: Defines where and how to save the data.
    -   `db_path` (str): The file path for the output ASE database.

-   **`FullConfig(BaseModel)`**: The root model.
    -   `system: Union[AlloySystemConfig, IonicSystemConfig, ...]`
    -   `validation: ValidationConfig`
    -   `storage: StorageConfig`

**Data Consumers and Producers:**

-   **Producers**: The primary producer of data is the user, who creates a YAML configuration file. This file is parsed and validated by our Pydantic schemas at runtime. The `Generator` classes are internal producers, creating ASE `Atoms` objects.
-   **Consumers**: The `PipelineRunner` consumes the `FullConfig` object. The `Generator` classes consume the `SystemConfig` part of the configuration. The `AseDBWrapper` consumes the list of ASE `Atoms` objects produced by the generators.

**Versioning and Extensibility:**

This schema-driven design is inherently extensible. To add a new system type (e.g., "MolecularCrystal"), a developer would:
1.  Create a new `MolecularCrystalSystemConfig` Pydantic model.
2.  Add it to the `Union` in the `FullConfig` model.
3.  Implement a corresponding `MolecularCrystalGenerator` class.
4.  Update the `GeneratorFactory` to recognize the new `system_type`.

This ensures that adding new features is a predictable process with a low risk of introducing regressions. The schemas act as a clear contract between different parts of the application.

## 4. Implementation Approach

The implementation will proceed in a logical, step-by-step manner, building the system from the data layer up to the user interface layer.

1.  **Project Scaffolding**: Create the directory structure as outlined in the System Architecture section. Initialize `__init__.py` files in all sub-packages.
2.  **Configuration Schemas**: Implement all Pydantic models in `config/schemas.py`. This is the first and most critical step. We will start with `ValidationConfig` and `StorageConfig`, then define the base `SystemConfig` and the specific implementations for `AlloySystemConfig` and `IonicSystemConfig`. We will include custom validators, for example, to check that alloy compositions sum to 1.0.
3.  **Utility Functions**: Implement the core validation logic in `utils/physics.py`. The `validate_structure` function will take an ASE `Atoms` object and a minimum distance and return `True` or `False`. This function will be pure and have no dependencies on the rest of the application, making it easy to unit test.
4.  **Storage Wrapper**: Implement the `AseDBWrapper` in `storage/db_wrapper.py`. This class will be simple, with a constructor that takes the database path and a `write_structures` method. This component can be developed and tested independently.
5.  **Generator Base Class and Factory**: Implement the `BaseGenerator` abstract class in `generators/base.py`. This will establish the contract that all other generators must follow. Then, implement the `GeneratorFactory` in `generators/factory.py`. This factory will have a single static method, `create_generator`, that takes a `SystemConfig` object and returns an instance of the correct generator subclass based on the `system_type` field.
6.  **Concrete Generators**: Implement the `AlloyGenerator` and `IonicGenerator` classes. This is where the core domain logic resides. We will use libraries like `ase.build` and `pymatgen` to create the initial crystal structures and then apply the necessary modifications (e.g., creating a supercell, introducing random atom placements for alloys). Each generator will call the `_validate` method from its base class on every structure it produces.
7.  **Pipeline Runner**: Implement the `PipelineRunner` in `pipeline/runner.py`. This class will tie everything together. Its `run` method will instantiate the factory, create the generator, call its `generate` method, and then pass the results to the `AseDBWrapper`.
8.  **CLI**: Finally, implement the command-line interface in `cli.py` using `click`. This script will be responsible for loading a YAML configuration file, parsing it into the `FullConfig` Pydantic model, and then kicking off the `PipelineRunner`.

## 5. Test Strategy

Testing in Cycle 1 will be focused on ensuring the correctness of the data models and the generated structures.

**Unit Testing Approach (Min 300 words):**

Unit tests will be created in a parallel `tests/unit` directory. We will test each component in isolation.
-   **Schemas**: For the Pydantic schemas in `config/schemas.py`, we will write tests that verify both valid and invalid data. For example, for `AlloySystemConfig`, we will test that a configuration with compositions summing to 1.0 is accepted, while a configuration where they sum to 0.9 is rejected with a `ValidationError`. We will test edge cases, such as empty element lists or negative numbers where positive ones are expected. This ensures our configuration system is robust.
-   **Generators**: Each generator will be tested extensively. We will not use mocks for ASE, but rather create small, well-defined test cases. For `AlloyGenerator`, we will request a 4-atom cell of CuAu and assert that the resulting `Atoms` object contains 2 Cu and 2 Au atoms. We will also iterate through all pairs of atoms in the generated structure and assert that their distance is greater than the configured `min_dist`. This directly tests the core functionality and validation logic of each generator.
-   **Utilities**: The `validate_structure` function in `utils/physics.py` will be tested with hand-crafted `Atoms` objects. We will create one object where two atoms are deliberately placed too close and assert that the function returns `False`. We will create another object that is valid and assert it returns `True`. This ensures the core physical validation logic is correct.

**Integration Testing Approach (Min 300 words):**

Integration tests will be located in `tests/integration` and will test the full workflow from the CLI down to the database output. We will use `pytest` fixtures and `click.testing.CliRunner` to invoke our CLI from within the tests.
-   **End-to-End Run**: A typical test will involve creating a temporary directory and writing a small YAML configuration file into it (e.g., for generating 5 random alloy structures of AgPd). The test will then use `CliRunner.invoke` to run the `main` command, passing the path to this temporary config file. The test will assert that the command exits with a status code of 0.
-   **Database Verification**: After the CLI command completes, the test will check that the output ASE database file (`test.db`) was created in the temporary directory. It will then use ASE's `ase.db.connect` to open the database and inspect its contents. The test will assert that the database contains exactly 5 rows. It will then read the first row, reconstruct the `Atoms` object, and verify that its chemical formula matches the expected 'AgPd' and that its atoms have the correct types. This verifies that all components, from config parsing to generation to storage, are working together correctly. We will create multiple such tests for different system types (alloy, ionic) to ensure broad coverage.
