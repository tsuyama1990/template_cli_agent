# CYCLE 01 SPEC.md: Core Functionality - Foundational Scaffolding

## 1. Summary

This document provides the detailed technical specification for the inaugural development cycle of the MLIP-AutoPipe project. The paramount objective of Cycle 01 is to construct the foundational scaffolding of the application, establishing a robust and extensible framework upon which all future functionality will be built. This cycle is meticulously scoped to include three critical components: a user-friendly and robust Command-Line Interface (CLI), a powerful and strictly type-safe configuration system, and the initial service for generating physically plausible atomic structures. Upon successful completion of this cycle, a user will be able to author a simple YAML file to define an alloy system, invoke the application from their command line, and receive a correctly formatted XYZ file containing a set of atomic coordinates that are guaranteed to be physically valid. To manage complexity and ensure a high-quality initial delivery, this cycle deliberately omits the more advanced aspects of the pipeline, such as the molecular dynamics exploration, intelligent sampling, and database storage. These features will be the focus of the subsequent cycle.

The key deliverables for Cycle 01 are a functional CLI entry point powered by the `typer` library, a comprehensive configuration module built on `Pydantic` that provides clear, actionable validation and error messages, and an initial `AlloyGenerator` class capable of creating randomized alloy structures based on precise user specifications. The successful implementation of these components is not merely a feature delivery; it is a critical validation of the project's core architectural principles. This foundational work will establish the definitive project structure, the patterns for dependency injection, the conventions for service-based logic, and the comprehensive testing strategy that will govern the remainder of the project's development. By focusing on a narrow but deep slice of functionality, we ensure that the initial codebase is stable, well-tested, and perfectly prepared for the addition of more complex scientific features in Cycle 02. It is, in essence, the pouring of a solid foundation before the rest of the structure is erected. This phase prioritizes architectural integrity and developer experience to ensure the long-term health and maintainability of the software.

## 2. System Architecture

The work in this cycle will establish the initial, canonical project structure within the `src/mlip_autopipec` directory. This structure is not arbitrary; it is a deliberate design that embodies the principles of modularity and separation of concerns. The focus is on creating the specific modules required to deliver the core functionality of this cycle: the CLI, the configuration management, and the first structure generation service. This structured approach ensures that the codebase is logical, navigable, and easy for new developers to understand. Each component will have a clearly defined location and responsibility, which is crucial for maintaining a clean architecture as the project grows in complexity. The test suite will mirror this structure, providing a one-to-one mapping between application code and test code, which facilitates test discovery and maintenance.

**File Structure (Cycle 01 Additions/Modifications in bold):**

```
src/
└── mlip_autopipec/
    ├── **__init__.py**
    ├── **cli.py**                 # Typer-based Command-Line Interface (Infrastructure Layer)
    ├── **config.py**              # Pydantic models for configuration management (Application Layer)
    ├── **domain/**
    │   └── **__init__.py**
    ├── **services/**              # Core Business Logic Layer
    │   ├── **__init__.py**
    │   └── **generation.py**      # Structure generation logic and validation
    └── **orchestrators/**         # High-level workflow coordination
        ├── **__init__.py**
        └── **workflow.py**        # The main WorkflowOrchestrator class
tests/
└── unit/
    ├── **__init__.py**
    ├── **test_cli.py**
    ├── **test_config.py**
    └── **test_generation.py**
```

This file structure provides a clear and robust separation of concerns, which is a cornerstone of the project's architecture:
-   `cli.py`: This file is the outermost layer of the application, responsible exclusively for handling user interaction and command-line parsing. It contains no business logic.
-   `config.py`: This module is dedicated to defining the data contracts for the application's configuration. It validates and type-annotates all incoming settings.
-   `orchestrators/workflow.py`: This module acts as the central coordinator, managing the high-level process flow and directing the sequence of operations.
-   `services/generation.py`: This is part of the core business logic layer, containing the detailed scientific algorithms and validation rules for creating atomic structures.
-   `tests/`: This directory contains the corresponding unit tests, ensuring that each component of the application is thoroughly verified and behaves as expected. This parallel structure is vital for maintaining high code coverage and reliability.

## 3. Design Architecture

The design for Cycle 01 is fundamentally centered around a "schema-first" philosophy, with the Pydantic library as its cornerstone. This approach mandates that all data structures, especially the user-facing configuration, are rigorously defined, validated, and self-documented before they are used in the application's logic. This design choice is motivated by a desire for robustness and clarity. By validating the entire configuration at the entry point of the application, we eliminate a significant source of potential runtime errors and remove the need for repetitive, defensive error-checking code within the deeper service layers. The primary focus is therefore on the `config.py` module, which will serve as the single source of truth for all configurable parameters that drive the application's behavior. This module will define the contract that the user must fulfill, and in return, it provides the rest of the application with a guarantee of type safety and data integrity.

**`config.py` - Pydantic-based Schema Blueprint:**

The application's configuration will be defined by a series of hierarchically nested Pydantic models. This structure allows for a clean, readable, and organized representation in the user's YAML file, while providing the application with a powerful, strongly-typed object model internally.

```python
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Any

class GenerationConfig(BaseModel):
    """Configuration for the initial structure generation phase."""
    generator_type: str = Field("Alloy", description="The type of generator to use (e.g., 'Alloy').")
    elements: List[str] = Field(..., min_items=1, description="A list of all element symbols to be included.")
    composition: Dict[str, float] = Field(..., description="A dictionary mapping element symbols to their fractional composition.")
    lattice_constant: float = Field(..., gt=0, description="The lattice constant of the primitive cell in Angstroms.")
    num_structures: int = Field(10, gt=0, description="The total number of unique, valid structures to generate.")
    supercell_size: List[int] = Field([3, 3, 3], min_items=3, max_items=3, description="The 3x3 matrix defining the supercell expansion.")

    @validator('composition')
    def composition_must_sum_to_one(cls, v: Dict[str, float]) -> Dict[str, float]:
        if not abs(sum(v.values()) - 1.0) < 1e-6:
            raise ValueError("The sum of composition values must be equal to 1.0")
        return v

    @validator('composition')
    def elements_must_match_composition(cls, v: Dict[str, float], values: Dict[str, Any]) -> Dict[str, float]:
        if 'elements' in values and set(v.keys()) != set(values['elements']):
            raise ValueError("The keys in 'composition' must exactly match the symbols in 'elements'")
        return v

class FullConfig(BaseModel):
    """The root configuration model for the entire MLIP-AutoPipe pipeline."""
    generation: GenerationConfig

# Utility function to be placed in this module
def load_config_from_yaml(path: str) -> FullConfig:
    """Loads a YAML configuration file from the given path and validates it against the schema."""
    # The implementation will use the PyYAML library to load the raw data
    # and then Pydantic's model_validate method to parse and validate it.
    pass
```

**Key Invariants, Constraints, and Validation Rules:**
-   The `composition` dictionary's floating-point values must sum to approximately 1.0. A custom Pydantic validator will be implemented to enforce this with a small tolerance.
-   The set of keys in the `composition` dictionary must be identical to the set of elements in the `elements` list. A cross-field validator in Pydantic will enforce this critical consistency check.
-   Numeric fields like `lattice_constant` and `num_structures` must be strictly positive, which will be enforced by Pydantic's `gt=0` field constraint.
-   The `supercell_size` must always be a list containing exactly three integers.

**Consumers and Producers of the Configuration:**
-   **Producer:** The end-user, who authors a `config.yaml` file on their local filesystem according to the defined schema.
-   **Consumer:** The `cli.py` module is the primary internal consumer. It will invoke the `load_config_from_yaml` utility function to read, parse, and validate the YAML file. Upon successful validation, it produces a `FullConfig` instance, which is then passed to the `WorkflowOrchestrator` to begin the execution of the pipeline.

## 4. Implementation Approach

The implementation for Cycle 01 will proceed in a logical, dependency-first sequence. We will start with the most fundamental component—the configuration—and build outwards towards the user-facing CLI. This ensures that each layer of the application is built upon a solid, tested foundation.

1.  **Implement `config.py` (The Foundation):**
    -   First, define the `GenerationConfig` and `FullConfig` Pydantic models exactly as specified in the Design Architecture blueprint. This includes all fields, type hints, and descriptions.
    -   Implement the custom validators within the `GenerationConfig` model. This is a critical step to ensure the integrity of the configuration. The first validator will check that the composition fractions sum to 1.0. The second, more complex validator will ensure that the elements listed are perfectly consistent with the keys of the composition dictionary.
    -   Implement the `load_config_from_yaml` utility function. This function will encapsulate the logic for reading the YAML file from disk using the `pyyaml` library and then safely parsing the raw dictionary into the `FullConfig` Pydantic model using `model_validate`. It must include robust `try...except` blocks to handle potential `FileNotFoundError` and `yaml.YAMLError`, providing user-friendly error messages for each case.

2.  **Implement `services/generation.py` (The Core Logic):**
    -   Define a Python `abc.ABC` (Abstract Base Class) named `BaseStructureGenerator`. It will declare a single abstract method, `generate() -> List[ase.Atoms]`, thereby defining the contract that all future generators must adhere to.
    -   Implement the `AlloyGenerator` class, ensuring it inherits from `BaseStructureGenerator`.
    -   The `__init__` method of `AlloyGenerator` will accept a fully validated `GenerationConfig` object as its sole argument, storing it for later use.
    -   The public `generate` method will orchestrate the structure creation. It will contain the high-level logic to first create a primitive cell (using `ase.build.bulk`), then expand it into a supercell (using `ase.build.make_supercell`), and finally, randomize the atomic placements according to the specified `composition`.
    -   A private helper method, `_randomize_composition`, will be implemented to handle the logic of randomly assigning atomic symbols to the sites in the supercell lattice.
    -   A second private method, `_is_physically_valid`, will implement the critical `overlap_check` to ensure no two atoms are unrealistically close. This function will take an `ase.Atoms` object and return a boolean. The main `generate` loop will use this method to discard any invalid structures it creates, continuing until it has generated the required `num_structures` of valid configurations.

3.  **Implement `orchestrators/workflow.py` (The Conductor):**
    -   Create the `WorkflowOrchestrator` class.
    -   Its `__init__` method will accept a `FullConfig` object.
    -   Implement a public `run` method. This method will contain the logic to:
        a.  Implement a simple factory pattern. It will inspect the `config.generation.generator_type` string and, based on its value, will instantiate the appropriate generator class (for this cycle, it will only handle the "Alloy" type).
        b.  It will then call the generator instance's `generate` method to receive the list of `ase.Atoms` objects.
        c.  Finally, it will use `ase.io.write` to serialize the list of atoms to an output file named `initial_structures.xyz` in the current working directory.

4.  **Implement `cli.py` (The Entry Point):**
    -   Create a `typer` application instance.
    -   Define a single command, `run`, that accepts a required `config_path` argument of type `typer.Argument`.
    -   The function implementing this command will be the final piece that connects everything. It will:
        a.  Call `config.load_config_from_yaml` with the path provided by the user.
        b.  Instantiate the `WorkflowOrchestrator` with the resulting `FullConfig` object.
        c.  Call the orchestrator's `run` method to execute the entire process.
        d.  Upon successful completion, it will use the `rich` library to print a formatted, user-friendly success message to the console.

## 5. Test Strategy

The testing strategy for Cycle 01 is of paramount importance as it establishes the quality standards and testing patterns for the entire project. All tests will be implemented using the `pytest` framework and will be co-located in the `tests/` directory, mirroring the application's package structure.

**Unit Testing Approach:**
The unit testing philosophy is to verify each component in complete isolation, using mocks to sever dependencies on other components or the filesystem. This ensures that tests are fast, deterministic, and precisely target the logic under review.
-   **`test_config.py`:** This will be the most comprehensive and detailed test suite in this cycle. It will rigorously validate the `load_config` function and the Pydantic models. We will create a suite of test YAML files in a temporary directory.
    -   **Success Path:** A test will load a perfectly valid `config.yaml` and assert that the returned object is a `FullConfig` instance with attributes that exactly match the file's content.
    -   **Failure Paths:** A series of tests will be written to cover all conceivable invalid configurations. This includes YAML files with missing required fields (e.g., no `elements` list), fields with incorrect data types (e.g., `lattice_constant` as a string), and values that violate the custom validation rules (e.g., compositions that do not sum to 1.0, or a mismatch between `elements` and `composition` keys). For each of these invalid cases, the test will use the `pytest.raises(ValidationError)` context manager to assert that Pydantic's validation is correctly triggered and that the error messages are informative. This suite acts as a crucial quality gate for all application inputs.
-   **`test_generation.py`:** This suite will focus on the `AlloyGenerator`. To ensure isolation, it will not rely on file I/O or the config loader. Instead, a valid `GenerationConfig` object will be instantiated directly within the test code and passed to the `AlloyGenerator`.
    -   The `generate` method will be called, and detailed assertions will be performed on the returned list of `ase.Atoms` objects. We will verify: (1) that the list contains the correct number of structures (`num_structures`), (2) that each individual structure has the correct total number of atoms corresponding to the supercell size, (3) that the chemical composition of each generated structure (e.g., the ratio of Cu to Au atoms) is statistically consistent with the requested `composition`, and (4) that the physical validation logic is working. This last point will be explicitly tested by creating a configuration with a very large atomic radius, forcing overlaps, and asserting that the generator correctly discards these invalid structures and still manages to produce the required number of valid ones.
-   **`test_cli.py`:** This suite will use the `typer.testing.CliRunner` to perform black-box tests on the command-line interface.
    -   It will test the success case by invoking the `run` command with a path to a valid temporary configuration file. The test will assert that the command exits with a success code (0), that a success message is printed, and that the expected output file (`initial_structures.xyz`) is created.
    -   It will also test the CLI's robustness by running the command with a path to a non-existent file and asserting that it exits with a non-zero exit code and prints a clear "File not found" error message.

**Integration Testing Approach:**
While unit tests verify the parts, integration testing verifies the whole. For Cycle 01, the goal is to ensure that all the newly created components—from the CLI to the file output—work together seamlessly.
A single, comprehensive integration test will be created to simulate a real user's workflow from beginning to end.
1.  **Setup:** The test function will use `pytest`'s `tmp_path` fixture to create a temporary, isolated directory. Inside this directory, it will programmatically create a `config.yaml` file with a complete and valid configuration for generating a simple binary alloy (e.g., a 2x2x2 supercell of CuAu).
2.  **Execution:** The test will instantiate the `typer.testing.CliRunner` and invoke the main CLI command (`run`), passing the path to the temporary `config.yaml`. Critically, no part of the application's internal logic will be mocked. The CLI will call the real `load_config`, which will instantiate the real `WorkflowOrchestrator`, which in turn will instantiate and run the real `AlloyGenerator`.
3.  **Verification:** After the command has finished executing, the test will perform a series of assertions to verify the complete outcome.
    -   It will first assert that the CLI command exited with a status code of 0.
    -   It will then check for the existence of the `initial_structures.xyz` file within the temporary directory.
    -   Finally, and most importantly, it will read the generated `initial_structures.xyz` file back into memory using `ase.io.read(..., index=':')`. This crucial step validates that the output is in a correct, standard, and readable format. Once loaded, it will perform detailed assertions on the resulting list of `ase.Atoms` objects, verifying that the number of structures and the physical properties of each structure (e.g., number of atoms, cell dimensions, chemical symbols) are all fully consistent with the parameters defined in the original input configuration file. This end-to-end test provides high confidence in the correctness of the entire application stack for this cycle.
