# Specification: MLIP-AutoPipe - Cycle 1

## 1. Summary

This document provides the detailed technical specification for the first implementation cycle of the MLIP-AutoPipe project. The primary goal of Cycle 1 is to establish the foundational architecture and a minimal, end-to-end runnable pipeline. This cycle is fundamentally about building the "chassis" of the application: a solid, testable skeleton upon which the complex scientific "engine" of Cycle 2 can be mounted. While this initial version will not contain the advanced scientific logic for exploration or sampling, its successful completion is a critical milestone that de-risks the project by proving the viability of the core software architecture. The core deliverables for this cycle are a functional command-line interface (CLI) for user interaction, a robust, self-validating configuration system based on Pydantic for defining workflows, and the basic scaffolding for the main application pipeline. This includes a simplified structure generator, a placeholder (mock) for the computationally intensive labeling engine, and a functional database wrapper for storing results. The successful completion of this cycle will result in a command-line tool that can be executed with a simple configuration, run without errors, and produce a correctly formatted, albeit scientifically simplistic, output database.

Specifically, the implementation will begin by defining the project's dependencies and structure in `pyproject.toml`. The centerpiece of this cycle's design is the schema-first approach. We will meticulously define a set of Pydantic models in `src/mlip_autopipec/config.py` that represent the entire user-facing configuration schema. This ensures that from day one, all user input is rigorously validated, preventing a large class of potential runtime errors. The main application logic will be encapsulated in a `WorkflowOrchestrator` class. This class will be responsible for interpreting the validated configuration object and sequentially invoking the different stages of the simplified pipeline: Generation, Labeling, and Storage. For the "Generation" stage, a basic `AlloyGenerator` will be implemented. This component will be capable of producing simple, randomized alloy structures based on user-specified elemental compositions. It will also incorporate essential physical validation checks, such as ensuring atoms do not overlap, establishing a pattern of physical correctness that will be followed by all future generators. The "Labeling" stage, which in Cycle 2 will involve complex and time-consuming DFT calculations, will be implemented here as a simple placeholder or "mock." This `MockLabelingEngine` will fulfill the required software interface by accepting an atomic structure and returning a data structure with dummy energy and force values. This is a crucial strategic decision, as it allows us to develop and test the entire data flow of the pipeline without the dependency on external simulation software. Finally, the "Storage" stage will be implemented with an `AseDBWrapper`, a dedicated class for handling all interactions with the ASE-compatible SQLite database, ensuring that data is written correctly and safely.

## 2. System Architecture

The focus of Cycle 1 is to construct the foundational file and class structure of the application. The components implemented in this cycle are the essential building blocks that define the software's architecture and enable future expansion. The architecture is designed to be modular, testable, and to enforce a clear separation of concerns, which is critical for the long-term maintainability of a scientific software project. Every file and module has a distinct and well-defined purpose. The directory structure is not arbitrary; it is a direct reflection of the layered design of the application, separating the user interface, business logic, and data access layers. This meticulous organization ensures that developers can easily locate relevant code, understand its function, and modify it with minimal risk of unintended side effects. For instance, all database-related code is confined to `ase_db_wrapper.py`, meaning any future changes to the storage mechanism will be localized to a single file. Similarly, the logic for creating atomic structures is isolated within `generators.py`, allowing for the addition of new, complex generators in Cycle 2 without impacting the core workflow. This disciplined approach is paramount for building a robust and extensible system.

**File Structure (Cycle 1 Focus):**
The files and directories to be created or modified in this cycle are marked in **bold**. This file tree represents the exact blueprint for the initial codebase.

```
.
├── **pyproject.toml**              # Project metadata, dependencies, and tool configuration.
└── src
    └── mlip_autopipec
        ├── **__init__.py**         # Makes the directory a Python package.
        ├── **cli.py**              # **UI Layer**: Handles all command-line interaction using Typer.
        ├── **config.py**           # **Schema Layer**: Defines the Pydantic models for configuration validation.
        ├── **main.py**             # **Application Entrypoint**: Contains the main function to run the orchestrator.
        ├── common
        │   ├── **__init__.py**
        │   ├── **pydantic_models.py** # **Domain Layer**: Defines core internal data structures like DFTResult.
        │   └── **exceptions.py**     # Defines custom exception classes for the application.
        └── core
            ├── **__init__.py**
            ├── **ase_db_wrapper.py** # **Infrastructure Layer**: Manages all interaction with the ASE SQLite DB.
            ├── **generators.py**     # **Service Layer**: Contains logic for generating atomic structures.
            ├── **labeling.py**       # **Service Layer**: Contains the placeholder engine for labeling structures.
            └── **workflow.py**       # **Orchestration Layer**: The core class that manages the pipeline execution.
```

This structure clearly separates concerns into distinct layers:
-   **`cli.py` (UI Layer):** This module's sole responsibility is to interact with the user. It parses command-line arguments, provides help text, and handles the initial loading of the configuration file. It contains no scientific or business logic. Its role is to be a thin wrapper that translates user commands into calls to the orchestration layer.
-   **`config.py` (Schema Layer):** This module is the public contract of the application. It defines, through Pydantic models, exactly what constitutes a valid configuration file. This is the single source of truth for all application settings.
-   **`workflow.py` (Orchestration Layer):** This module contains the high-level business logic. The `WorkflowOrchestrator` class within it directs the overall flow of the pipeline, calling the various service-layer components in the correct order and managing the flow of data between them.
-   **`core/*.py` (Service and Infrastructure Layers):** These modules are the workhorses of the application. Each module provides a specific service: `generators.py` creates things, `labeling.py` calculates things, and `ase_db_wrapper.py` saves things. They are designed to be independent of each other and of the UI layer, which makes them highly reusable and easy to unit-test in isolation. This modular design is what will allow for the seamless addition of complex features in Cycle 2.

## 3. Design Architecture

This cycle is fundamentally about establishing a robust, schema-driven foundation. The "Schema-First" principle is the most important aspect of the design architecture. All data, from the user's configuration file to the data objects passed between internal components, will be defined by Pydantic models. This approach provides numerous benefits, including automatic and transparent validation, clear and explicit data contracts, and a significant reduction in common data-related bugs. It forces a level of design discipline that is often absent in scientific scripts, leading to a much more professional and maintainable codebase. By defining the structure of our data before implementing the logic that acts upon it, we create a stable and predictable environment. This means that any function or method can be confident that the data it receives is already validated and conforms to the expected structure, allowing the developer to focus purely on the core logic of that function. This upfront investment in data modeling pays significant dividends in terms of reduced debugging time and increased code clarity. The Pydantic schemas also serve as a form of living documentation; by reading the models in `config.py`, a new user can immediately understand all the available settings and their constraints.

**Pydantic Schema Design (`config.py` and `common/pydantic_models.py`):**
The following code blocks are not just examples; they are the exact, complete blueprints for the Pydantic models to be implemented in Cycle 1.

**In `src/mlip_autopipec/config.py`:**
This file will define the user-facing configuration schema. This is the contract with the user.

```python
# config.py
from pydantic import BaseModel, Field, validator
from typing import List, Dict

class SystemConfig(BaseModel):
    elements: List[str] = Field(..., min_items=1)
    composition: Dict[str, float]
    lattice_constant: float = Field(4.0, gt=0)
    crystal_structure: str = 'fcc'

    @validator('composition')
    def check_composition_sum(cls, v):
        if not abs(sum(v.values()) - 1.0) < 1e-6:
            raise ValueError('Composition probabilities must sum to 1.0')
        return v

class GenerationConfig(BaseModel):
    num_structures: int = Field(..., gt=0, description="The number of structures to generate.")
    supercell_size: int = Field(3, ge=1)
    apply_rattle: bool = True
    rattle_strength: float = Field(0.01, ge=0)

class LabelingConfig(BaseModel):
    # Placeholder for Cycle 1. In Cycle 2, this will contain DFT settings.
    engine: str = Field('mock', description="The labeling engine to use.")

class StorageConfig(BaseModel):
    db_path: str = "output_database.db"

class FullConfig(BaseModel):
    """ The root model for the entire configuration file. """
    system: SystemConfig
    generation: GenerationConfig
    labeling: LabelingConfig
    storage: StorageConfig
```

**In `src/mlip_autopipec/common/pydantic_models.py`:**
This file will define the internal data models that are passed between components.

```python
# common/pydantic_models.py
from pydantic import BaseModel, Field
from typing import List, Optional
import numpy as np

class DFTResult(BaseModel):
    """
    Represents the result of a single-point calculation.
    This model acts as a clear data contract between the LabelingEngine
    and the rest of the application.
    """
    class Config:
        arbitrary_types_allowed = True

    energy: float
    forces: np.ndarray
    stress: Optional[np.ndarray] = None

    @validator('forces', 'stress', pre=True)
    def to_numpy_array(cls, v):
        if isinstance(v, list):
            return np.array(v)
        return v
```
The design decisions here are deliberate. In `SystemConfig`, a custom validator is included to ensure the composition percentages sum to 1.0, catching a common user error. In `DFTResult`, `arbitrary_types_allowed = True` is necessary to allow for NumPy arrays, a common requirement in scientific computing. These models establish the clear data contracts that underpin the entire application. The `WorkflowOrchestrator` will be initialized with a `FullConfig` object. The `AlloyGenerator`'s `generate` method will return a `List[ase.Atoms]`. The `MockLabelingEngine`'s `run` method will return a `DFTResult` object. And finally, the `AseDBWrapper`'s `write_atoms` method will accept an `ase.Atoms` object and a `DFTResult` object. This strict, type-hinted, and validated flow of data is the core of this cycle's architectural design.

## 4. Implementation Approach

The implementation for Cycle 1 will proceed in a logical, bottom-up fashion, starting with the foundational data models and dependencies and culminating in the user-facing command-line interface. This ensures that each component is built upon a solid and already-tested foundation.

1.  **Project Setup (`pyproject.toml`):** The very first step is to establish the project's environment. The `pyproject.toml` file will be created and populated with all the necessary dependencies for this cycle. This includes `pydantic` for data modeling, `typer` for the CLI, `pyyaml` for parsing the configuration file, `ase` and `numpy` for the scientific data structures, and `pytest` for testing. An entry point for the CLI, `mlip-autopipec = "mlip_autopipec.cli:app"`, will also be defined in the `[project.scripts]` section.

2.  **Schema Implementation (`config.py`, `pydantic_models.py`):** Following the schema-first principle, the next step is to translate the Pydantic schema designs from the section above into actual code. The classes `SystemConfig`, `GenerationConfig`, `LabelingConfig`, `StorageConfig`, and `FullConfig` will be implemented in `src/mlip_autopipec/config.py`. The `DFTResult` model will be implemented in `src/mlip_autopipec/common/pydantic_models.py`.

3.  **Core Component Implementation (`core/`):** With the data contracts in place, we can now build the core service and infrastructure components.
    *   **`AseDBWrapper` (`ase_db_wrapper.py`):** This class will be implemented to handle all database operations. Its constructor will accept the database file path as a string. The primary method will be `write_atoms(self, atoms: ase.Atoms, dft_result: DFTResult)`. This method will use `ase.db.connect` to open the database. Crucially, it will use a context manager (`with db.managed_connection():`) to ensure the database connection is automatically and safely closed, even if errors occur. Inside the `with` block, it will call `db.write(atoms, data=dft_result.model_dump())`. The use of `model_dump()` serializes the Pydantic model into a dictionary suitable for storage in the ASE DB's `data` column.
    *   **`generators.py`:** First, an abstract base class `BaseGenerator` will be defined using `abc.ABC`. It will have one abstract method, `generate(self) -> List[ase.Atoms]`. Then, the `AlloyGenerator` will be implemented, inheriting from `BaseGenerator`. Its constructor will accept the `SystemConfig` and `GenerationConfig` objects. The `generate` method will contain the logic: it will use `ase.build.bulk` to create a primitive cell, then use `ase.build.make_supercell` to expand it. It will then iterate through the atoms in the supercell, randomly assigning their chemical symbols based on the probabilities defined in the `composition` dictionary. A final validation step within this method will check for atomic overlaps. The method will return a list containing the requested `num_structures` of `ase.Atoms` objects.
    *   **`labeling.py`:** This module will contain the placeholder `MockLabelingEngine`. The class will have a single method, `run(self, atoms: ase.Atoms) -> DFTResult`. For Cycle 1, this method will perform no calculations. It will simply instantiate and return a `DFTResult` object with default, non-physical values (e.g., `energy=0.0`, `forces=np.zeros((len(atoms), 3))`). This mock is essential for decoupling the pipeline logic from the yet-to-be-implemented DFT interface.

4.  **Orchestration (`workflow.py`):** Now we implement the `WorkflowOrchestrator`. Its constructor will accept a `FullConfig` object. Its main public method, `run_workflow()`, will execute the simplified pipeline logic:
    a. Log the start of the workflow.
    b. Instantiate the `AlloyGenerator` using the `system` and `generation` sections of its config.
    c. Call `generator.generate()` to create the list of initial structures.
    d. Instantiate the `MockLabelingEngine` and the `AseDBWrapper`.
    e. Begin a loop through the generated structures. Inside the loop, for each structure:
        i. Call the `labeling_engine.run()` to get the dummy `DFTResult`.
        ii. Call the `db_wrapper.write_atoms()` to save the structure and its dummy result.
    f. Log the successful completion of the workflow.

5.  **CLI and Entrypoint (`cli.py`, `main.py`):** Finally, we create the user-facing entry point.
    *   In `cli.py`, a `typer.Typer` app will be created. A main function, `run`, will be decorated with `@app.command()`. This function will accept a single argument: `config_path: Path = typer.Option(..., help="Path to the YAML configuration file.")`.
    *   Inside `run`, a helper function will be called to load and parse the YAML file. This function will use `pyyaml` to load the file and then `FullConfig.model_validate(data)` to parse it into the Pydantic model. This `model_validate` call is where the automatic validation happens.
    *   The `main.py` file will instantiate the `WorkflowOrchestrator` with the validated config object and call its `run_workflow()` method. The CLI `run` function will call this main application logic, wrapped in a `try...except` block to catch any potential errors and print a user-friendly message.

## 5. Test Strategy

Testing in Cycle 1 is arguably the most critical testing phase of the project, as it validates the core architecture and ensures that the foundation is stable, reliable, and correct. A robust test suite at this stage will prevent foundational flaws from propagating into the more complex features of Cycle 2. The strategy is two-pronged: exhaustive unit tests to verify each component in isolation, and a comprehensive integration test to verify that the components work together as a cohesive whole.

**Unit Testing Approach (Min 300 words):**
The primary goal of unit testing in this cycle is to achieve a high degree of confidence in each individual module's correctness, independent of its collaborators. We will use `pytest` as the test runner and `unittest.mock` for creating test doubles where necessary. All unit tests will reside in the `tests/unit/` directory.

*   **`test_config.py`:** This test suite is paramount for ensuring the robustness of our user-facing API. It will not test the logic of the application, but rather the Pydantic data models themselves. We will create a series of test cases, each representing a different user input scenario. For example:
    *   `test_valid_config_loads_successfully`: This will use a string containing a perfectly valid YAML configuration and assert that `FullConfig.model_validate(yaml.safe_load(config_string))` executes without raising any exceptions.
    *   `test_missing_required_field_raises_validation_error`: This test will use a YAML string where a critical field, like `elements` in the `system` block, is missing. It will use `pytest.raises(ValidationError)` to assert that this specific exception is thrown.
    *   `test_incorrect_data_type_raises_validation_error`: Here, a field will have the wrong type, such as `num_structures: "five"`. Again, we will assert that a `ValidationError` is raised, and we can even inspect the exception message to ensure it is informative.
    *   `test_custom_validator_for_composition`: This will specifically test our custom validator. A test case will provide a composition like `{'Fe': 0.6, 'Pt': 0.5}` (which sums to 1.1) and assert that our custom `ValueError` is raised.

*   **`test_generators.py`:** To test the `AlloyGenerator`, we need to control the randomness for reproducibility. Each test will start with `np.random.seed(42)`. The main test will instantiate the generator with a fixed configuration (e.g., an 8-atom cell of 50/50 SiGe). The assertions will be precise: (1) The length of the returned list of `ase.Atoms` objects must be exactly `num_structures`. (2) For each `Atoms` object, `atoms.get_chemical_symbols().count('Si')` must equal `4`, and the same for `Ge`. (3) A loop will check all pairwise atomic distances to ensure none are smaller than a physically realistic cutoff, proving the overlap check is functional.

*   **`test_ase_db_wrapper.py`:** This test must not leave any `test.db` files on the system. The use of `pytest`'s `tmp_path` fixture is essential. A test function `test_write_and_read_back` will be created. It will instantiate the `AseDBWrapper` with a path to a database file within `tmp_path`. It will then create a sample `ase.Atoms` object and a dummy `DFTResult` object. After calling `write_atoms`, the test will immediately use `ase.db.connect` to open the same temporary database file. It will read the first (and only) row using `db.get(1)`. It will then perform detailed assertions: `assert row.energy == 0.0`, `assert np.array_equal(row.forces, np.zeros((len(atoms), 3)))`, and `assert row.symbols == atoms.get_chemical_symbols()`. This confirms that the data serialization and database writing process is working perfectly.

**Integration Testing Approach (Min 300 words):**
While unit tests are essential, they cannot guarantee that the components will work together correctly. The integration test is designed to fill this gap by simulating a real user running the application from the command line and verifying the final, concrete output. A single, powerful integration test will be created in `tests/integration/test_cli_integration.py`.

*   **Test Setup:** The test function will use the `typer.testing.CliRunner` to invoke the CLI without needing a separate subprocess. It will also use the `tmp_path` fixture to create a temporary, isolated working directory for the test run. This is crucial for ensuring that the test does not depend on or interfere with the state of the developer's file system.

*   **Arrange Phase:** Inside the test function, the first step is to "arrange" the test conditions. This involves programmatically creating a `config.yaml` file inside the `tmp_path` directory. The content of this YAML file will be carefully crafted for the test. It will specify a small number of structures (e.g., 5) to keep the test fast. Crucially, the `db_path` within this config file will be set to an absolute path also inside the `tmp_path` (e.g., `f"{tmp_path}/integration_test.db"`).

*   **Act Phase:** The next step is to "act" by running the application. The `CliRunner` is used for this: `result = runner.invoke(app, ["run", "--config", str(config_file_path)])`. This line executes the `mlip-autopipec run` command in-process and captures all of its outputs (exit code, stdout, stderr, and any exceptions).

*   **Assert Phase:** The final and most important phase is to "assert" that the outcome was correct. The assertions will be multi-faceted:
    1.  **CLI Success:** `assert result.exit_code == 0` and `assert result.exception is None`. This is the first check; if the application crashed, the test fails here.
    2.  **File System Side Effects:** The test will then check the file system. `db_path = tmp_path / "integration_test.db"; assert db_path.exists()`. This confirms that the application produced the expected artifact.
    3.  **Database Content Verification:** This is the deepest level of verification. The test will use `ase.db.connect(db_path)` to open the database that the application just created. It will then perform assertions against the database content itself: `db = ase.db.connect(db_path); assert len(db) == 5`. This confirms that the correct number of structures were generated and saved. It can even go further and read a row to ensure the dummy energy value is correct.

This end-to-end test, from simulating a user's command to introspecting the final database artifact, provides extremely high confidence that the core pipeline is wired correctly and is fully functional.
