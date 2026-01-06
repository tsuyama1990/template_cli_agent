# Specification: MLIP-AutoPipe Cycle 1

## 1. Summary

This document provides the exhaustive technical specification for the inaugural implementation cycle of the MLIP-AutoPipe project. The paramount objective of Cycle 1 is to construct the foundational software framework, culminating in a functional, end-to-end command-line tool capable of executing a simplified, yet complete, data generation workflow. This initial cycle is strategically focused on establishing the core architecture, defining the primary data structures, and implementing the essential components of the pipeline. By creating a solid, well-tested, and reliable base, we pave the way for the seamless integration of more advanced and scientifically complex features in subsequent cycles. The scope is deliberately and carefully limited to the most critical-path functionalities to ensure that the core pipeline is demonstrably robust, thoroughly tested, and dependable before we introduce additional layers of complexity. This approach mitigates development risk and provides a valuable, functional tool at the earliest possible stage.

The principal deliverable for this cycle is a command-line interface (CLI) application, professionally implemented using the Typer library for its robustness and user-friendly nature. This CLI will serve as the primary user entry point, accepting a single argument: a path to a detailed configuration file in YAML format, which declaratively defines a specific data generation task. Upon invocation, the CLI will orchestrate a sequential, four-stage pipeline: Generation, Exploration, Sampling, and Storage.

The **Generation** stage will be implemented with a single, highly robust `AlloyGenerator`. This component will be responsible for creating random binary alloy structures within a supercell, based on user-defined parameters such as chemical elements, precise compositions, and the underlying crystal structure (e.g., FCC, BCC). A critical feature of this generator is the incorporation of rigorous physical validation checks. These checks, such as ensuring no two atoms are closer than a physically realistic bonding distance, are essential to prevent impossible or unstable starting configurations, thereby ensuring the quality of the entire downstream simulation.

The **Exploration** stage, for this cycle, will utilize a pure Molecular Dynamics (MD) engine. It will consciously omit the more complex hybrid Monte Carlo features planned for Cycle 2, focusing instead on the core task of running stable MD simulations. To facilitate rapid development and testing of the parallel simulation infrastructure, this stage will use a simple, fast, and reliable classical potential, specifically the Effective Medium Theory (EMT) potential which is a standard component of the ASE library. This choice allows us to fully build out and debug the parallel execution logic without the significant computational overhead and added complexity of full-scale machine learning models.

The **Sampling** stage will be implemented with a straightforward `RandomSampler`. This component will be responsible for selecting a specified number of frames (atomic configurations) at random from the MD trajectories generated in the previous stage. While less sophisticated than the descriptor-based sampling methods planned for the future, this provides the essential functionality to reduce the vast amount of trajectory data to a manageable and useful size for the final dataset.

Finally, the **Storage** stage will be fully implemented, providing a durable and efficient `AseDBWrapper`. This service will be responsible for saving the final sampled atomic structures and their associated metadata (potential energy, forces) into a structured, persistent, and queryable SQLite-based ASE database. This ensures that the final output of the pipeline is a high-quality, ready-to-use data artifact, perfectly formatted for subsequent MLIP training tasks.

## 2. System Architecture

The architecture for Cycle 1 is meticulously designed to establish the fundamental file and module structure that will serve the project for its entire lifecycle. All source code will be located within the `src/mlip_autopipec` directory, a standard practice that establishes it as a modern, installable Python package. The design philosophy is rooted in the principle of separation of concerns, logically dividing the codebase into distinct areas of responsibility: configuration, domain models, services (which encapsulate the core scientific and business logic), and the user interface (the CLI). This structured approach is paramount for ensuring the project is maintainable, scalable, and testable.

The primary user entry point will be `cli.py`, which will define the Typer application and its commands. This file will be responsible for parsing command-line arguments, loading the user's YAML configuration file, and validating it against the defined Pydantic schema. The core workflow orchestration, however, will be delegated to the `PipelineRunner` class, located in `services/orchestration.py`. This orchestrator is the heart of the application, and it will sequentially invoke the specialized services for each of the four pipeline stages. The contracts between all components are rigorously defined by the Pydantic models in `config/models.py`. These models serve as the single source of truth for the application's configuration schema, ensuring that all inputs are validated at runtime before any computation begins.

A key aspect of this architecture is the decoupling of the core logic from the user interface. The `PipelineRunner` and its subsidiary services are designed to be completely independent of the CLI. They operate on a strongly-typed Pydantic `FullConfig` object. This means that the core logic is agnostic as to whether that configuration object was created from a YAML file by the CLI, generated by a future Web UI, or constructed programmatically within a test suite. This modularity is a critical enabler for the project's long-term extensibility and testability.

The `shared` directory is another important architectural element. It will contain utility functions that are broadly applicable across multiple services, such as functions for applying physical constraints (`physics.py`) or common helper functions for interacting with ASE `Atoms` objects (`ase_utils.py`). This promotes code reuse and helps to avoid circular dependencies between the service modules, leading to a cleaner and more maintainable codebase.

**File Structure and Code Blueprints (Cycle 1):**

The files and directories to be created or modified in this cycle are marked in **bold**. This blueprint provides a clear guide for the development process.

```
.
├── pyproject.toml
└── src
    └── mlip_autopipec
        ├── __init__.py
        ├── **cli.py**                # Typer-based CLI. Defines `run` command.
        │   └── `def run(config_file: Path):`
        │       # 1. Load and parse YAML into FullConfig Pydantic model.
        │       # 2. Instantiate PipelineRunner.
        │       # 3. Call runner.execute().
        │       # 4. Handle exceptions and print user feedback.
        ├── config
        │   ├── __init__.py
        │   └── **models.py**         # Pydantic models for the entire configuration schema.
        │       └── `class FullConfig(BaseModel):`
        │           # Contains nested models for system, generator, etc.
        ├── domain
        │   ├── __init__.py
        │   └── models.py         # Placeholder for future non-config domain models.
        ├── services
        │   ├── __init__.py
        │   ├── **orchestration.py**    # Contains the main pipeline orchestrator.
        │   │   └── `class PipelineRunner:`
        │   │       `def execute(self, config: FullConfig):`
        │   │           # Calls Generator, Explorer, Sampler, Storage in sequence.
        │   ├── **generators.py**       # BaseGenerator ABC and AlloyGenerator.
        │   │   └── `class AlloyGenerator(BaseStructureGenerator):`
        │   │       `def generate(self, ...) -> list[ase.Atoms]:`
        │   │           # Uses ase.build.bulk and random substitution.
        │   ├── **exploration.py**      # Explorer and MDRunner logic.
        │   │   └── `class Explorer:`
        │   │       `def run(self, atoms_list: list[ase.Atoms], ...):`
        │   │           # Uses ProcessPoolExecutor to run _run_single_md in parallel.
        │   ├── **sampling.py**         # BaseSampler ABC and RandomSampler.
        │   │   └── `class RandomSampler(BaseSampler):`
        │   │       `def sample(self, trajectories: list[Path], ...):`
        │   │           # Reads trajectory files and randomly selects frames.
        │   └── **storage.py**          # AseDBWrapper for database interaction.
        │       └── `class AseDBWrapper:`
        │           `def write(self, atoms_list: list[ase.Atoms]):`
        │               # Connects to SQLite DB and writes atoms.
        └── shared
            ├── __init__.py
            ├── **ase_utils.py**        # Helper functions for working with ASE.
            └── **physics.py**          # Physical constraint validation functions.
                └── `def check_atomic_overlap(atoms: ase.Atoms, min_dist: float):`
```

This detailed architectural blueprint not only defines the structure but also provides a high-level view of the key classes and their responsibilities, ensuring a coherent and organized development process.

## 3. Design Architecture

The design of the MLIP-AutoPipe framework is fundamentally centered around a schema-first, service-oriented approach. This design philosophy leverages the power of Pydantic for robust, declarative data validation and for establishing clear, statically-typed contracts between all components of the system. This meticulous approach is critical for ensuring that the data flowing through the complex pipeline is well-defined, correct, and consistent at every step, which dramatically reduces the likelihood of runtime errors and simplifies debugging.

**Pydantic-based Schema Design:**

The absolute foundation of the system's configuration is the `FullConfig` Pydantic model, which will be centrally defined in `config/models.py`. This model serves as the single, unambiguous source of truth for all configurable parameters within the application. It is designed as a composite object, constructed from several smaller, more focused models, each representing a specific logical part of the pipeline. This compositional approach makes the configuration schema modular and easy to understand.

*   `SystemConfig(BaseModel)`: This model defines the core physical system to be simulated.
    *   `elements`: `list[str]` - A list of chemical symbols (e.g., `['Cu', 'Au']`). This will be validated to ensure they are valid chemical symbols.
    *   `composition`: `dict[str, float]` - A mapping from a chemical symbol to its fractional composition. A custom `@model_validator` will be implemented to ensure the values of this dictionary sum to approximately 1.0.
    *   `generator_type`: `Literal['alloy']` - For Cycle 1, this provides a strict enumeration, allowing only 'alloy'. This prevents users from attempting to use features that are not yet implemented.
    *   `num_initial_structures`: `int` - The number of distinct seed structures to generate. This will use `Field(gt=0)` to ensure it is a positive integer.
*   `GeneratorConfig(BaseModel)`: This model holds parameters specific to the structure generation process.
    *   `supercell_size`: `int` - The size of the supercell to create (e.g., 3 for a 3x3x3 supercell). Will also be validated to be a positive integer.
*   `ExplorationConfig(BaseModel)`: This model contains all parameters for the Molecular Dynamics simulation stage.
    *   `temperature_k`: `float` - The target simulation temperature in Kelvin. Validation `Field(gt=0)` ensures physical correctness.
    *   `num_steps`: `int` - The total number of MD steps to perform for each trajectory, validated to be positive.
*   `SamplingConfig(BaseModel)`: Parameters for the data sampling stage.
    *   `sampler_type`: `Literal['random']` - Strictly allows only 'random' for this cycle.
    *   `num_samples`: `int` - The number of structures to select from each trajectory, validated to be positive.
*   `FullConfig(BaseModel)`: The top-level model that aggregates all other configuration models.
    *   `system: SystemConfig`
    *   `generator: GeneratorConfig`
    *   `exploration: ExplorationConfig`
    *   `sampling: SamplingConfig`
    *   `db_path: Path` - The file path for the output ASE database. The path will be resolved relative to the configuration file location.

These Pydantic models are the primary data carriers for all configuration data. The key consumer of the `FullConfig` object is the `PipelineRunner`, which will unpack it and pass the relevant sub-models (e.g., `ExplorationConfig`) to each specialized service. The primary producer is the user, who authors a YAML file that must conform to this rigidly defined schema. The use of Pydantic provides automatic, declarative, and user-friendly validation out-of-the-box. If a user provides invalid input, the application will fail fast with a clear error message before any computationally expensive tasks are initiated. This design prevents the system from ever entering an invalid state and significantly enhances the user experience by reducing debugging time. For future extensibility, this design is ideal. Adding a new generator type is as simple as adding its name to the `Literal` type and creating a new configuration model for its parameters, without breaking any existing functionality.

## 4. Implementation Approach

The implementation of Cycle 1 will proceed in a logical, bottom-up sequence, building the system from the foundational data layer up to the user-facing interface. This methodical approach ensures that each component can be thoroughly unit-tested as it is developed, before it is integrated into the larger pipeline. This strategy minimizes integration issues and simplifies debugging.

**Step 1: Pydantic Configuration Models (`config/models.py`)**
The first and most crucial step is the implementation of all Pydantic models as described in the Design Architecture section. This will take place in `src/mlip_autopipec/config/models.py`. This involves defining all fields with their strict types, default values where appropriate, and declarative validators using `pydantic.Field` (e.g., `Field(gt=0)`). A custom `@model_validator` will be implemented for the `SystemConfig` model to enforce the constraint that the composition dictionary values must sum to 1.0. This initial step establishes the immutable "contract" for all configuration data that will flow through the system.

**Step 2: Storage Service (`services/storage.py`)**
With the data contracts defined, the next step is to implement the final output stage. The `AseDBWrapper` class will be implemented in `src/mlip_autopipec/services/storage.py`. This class will be designed as a context manager, using `__enter__` and `__exit__` to correctly handle the opening and closing of the connection to the ASE database. It will provide a simple, high-level public method, `write_atoms(atoms_list: list[ase.Atoms])`, which will iterate through a list of `ase.Atoms` objects and write each one to the database. This component is developed early to provide a concrete target for the pipeline and to allow subsequent components to save their results for testing and inspection.

**Step 3: Generator Service (`services/generators.py` and `shared/physics.py`)**
Next, the initial data creation stage will be built. The `BaseStructureGenerator` abstract base class will be defined in `generators.py` to establish the interface. Then, the concrete `AlloyGenerator` will be implemented. Its logic will use `ase.build.bulk` to create a primitive unit cell, followed by `ase.build.make_supercell` to expand it to the desired size. It will then iterate through the atoms and randomly replace their chemical species according to the specified composition. Concurrently, the physical validation logic, such as the `check_atomic_overlap` function, will be implemented in `shared/physics.py`. The `AlloyGenerator` will call this validation function on every structure it creates before returning the final list.

**Step 4: Exploration Service (`services/exploration.py`)**
The core computational component, the `Explorer` service, will be implemented in `exploration.py`. It will contain the main public method, `run(atoms_list: list[ase.Atoms], ...)` that takes the list of seed structures. This method will use `concurrent.futures.ProcessPoolExecutor` to map each input structure to a private worker function, `_run_single_md`, achieving parallel execution. This worker function is where the simulation logic resides. It will first initialize the velocities on the `Atoms` object using `ase.md.velocitydistribution.MaxwellBoltzmannInit`. Then, it will attach an ASE `EMT` calculator. Finally, it will run a `Langevin` dynamics simulation for the specified number of steps, writing the trajectory to a temporary XYZ file using `ase.io.Trajectory`.

**Step 5: Sampling Service (`services/sampling.py`)**
With the ability to generate trajectories, the `RandomSampler` will be implemented in `sampling.py`. Its public `sample` method will accept a list of file paths to the trajectory files produced by the explorer. For each trajectory, it will read all the frames into memory, and then use Python's `random.sample` function to select `num_samples` frames. It will aggregate the results from all trajectories and return a single flat list of `ase.Atoms` objects.

**Step 6: Orchestration Service (`services/orchestration.py`)**
With all the individual pipeline stages now complete, the `PipelineRunner` will be implemented in `orchestration.py`. Its primary `execute` method will be the glue that holds the application together. It will instantiate each of the required services (Generator, Explorer, Sampler, Storage). It will then call them in the correct sequence, carefully managing the data flow between them: the list of atoms from the generator is passed to the explorer; the list of trajectory file paths from the explorer is passed to the sampler; and the final list of sampled atoms is passed to the storage service.

**Step 7: CLI Implementation (`cli.py`)**
Finally, with the entire backend logic complete, the user interface will be implemented in `cli.py` using Typer. A single command, `run`, will be created. It will accept one required option: `--config-file`, which will be a `pathlib.Path` object. The command's logic will be responsible for loading the specified YAML file, and crucially, parsing it into the `FullConfig` Pydantic model. This parsing step implicitly triggers Pydantic's validation. If validation passes, the command will then instantiate and call the `PipelineRunner` with the validated config object. Rich feedback will be printed to the console to inform the user of the pipeline's progress.

## 5. Test Strategy

The testing strategy for Cycle 1 is designed to be rigorous and multi-faceted, focusing on ensuring the correctness, robustness, and reliability of the core components and their seamless integration. We will employ a combination of unit testing to validate components in isolation and integration testing to validate the end-to-end workflow.

**Unit Testing Approach:**
The unit testing strategy is designed to validate each individual component of the pipeline in complete isolation, ensuring that its specific responsibilities are met correctly. We will use the `pytest` framework as our test runner and the `unittest.mock` library for creating mock objects to isolate the components under test from their dependencies.

For the **Generator Service**, we will write a comprehensive suite of tests for the `AlloyGenerator`. These tests will verify several key properties of the generated structures. One test case will assert that the generated `ase.Atoms` objects have the exact atomic composition that was requested in the configuration (e.g., a 50/50 mix of Cu and Au). Another test will verify that the total number of atoms in the structure and the dimensions of the simulation cell are correct based on the primitive cell type and the supercell expansion factor. A particularly critical set of tests will focus on the physical validation logic located in `shared/physics.py`. We will programmatically create invalid `ase.Atoms` objects with deliberately overlapping atoms and assert that the `check_atomic_overlap` function correctly raises a `ValueError`, ensuring our physical validation is working as expected.

For the **Storage Service**, the `AseDBWrapper` will be tested against a temporary database file created on-the-fly using `pytest`'s `tmp_path` fixture. These tests will confirm that the database file is correctly created upon first write, that single and multiple `Atoms` objects can be written to it without error, and, most importantly, that the data read back from the database is identical to the original data. This includes verifying positions, cell parameters, atomic numbers, and attached metadata like potential energy.

The **Configuration Models** in `config/models.py` are a cornerstone of the application's robustness. They will be tested by attempting to instantiate them with various types of invalid data. For example, we will test with a composition that does not sum to 1.0, a negative simulation temperature, or an invalid chemical symbol in the `elements` list. In each case, we will assert that a `pydantic.ValidationError` is raised, and we will inspect the error message to ensure it is clear and informative for the user.

The **CLI** itself will be tested using `typer.testing.CliRunner`. We will simulate command-line invocations with different arguments. We'll test the success path by providing a valid configuration file and asserting that the application exits with a code of 0. We will also test failure paths, such as providing a path to a non-existent configuration file or a file with invalid YAML syntax, and assert that the application exits with a non-zero error code and prints a helpful, user-friendly error message to the console. The `PipelineRunner` will be mocked in these specific tests to ensure we are only testing the CLI's logic (argument parsing and runner invocation) and not the entire pipeline.

**Integration Testing Approach:**
The integration testing approach for Cycle 1 is focused on verifying that the four distinct stages of the pipeline (Generation, Exploration, Sampling, Storage) function together correctly as a cohesive, end-to-end workflow. We will create a single, comprehensive integration test that exercises the entire application through its primary public interface: the command line. This test will be located in `tests/integration/test_pipeline.py`.

The setup for this test will be crucial. It will use `pytest`'s `tmp_path` fixture to create a temporary, isolated directory for the test run. Inside this directory, the test will programmatically create a minimal, valid YAML configuration file. This configuration will define a simple, fast-to-compute system to ensure the test completes quickly within a CI environment. For example, it might specify a 2x2x2 supercell of a Cu-Au alloy, to be simulated using the fast EMT potential for a very short duration (e.g., only 10 MD steps), and from which we will sample 5 frames.

The core of the test will then use `typer.testing.CliRunner` to invoke the `mlip-autopipec run` command, passing the path to this dynamically created configuration file. The assertions will be performed on the results and side effects of this command. First, the test will check that the command completes with an exit code of 0, which indicates a successful run. Second, it will verify that the specified output database file (e.g., `output.db`) has been created in the temporary directory. The most important assertion will be on the content of this database. The test will use the ASE database library to connect to the output database and will assert that it contains the correct number of entries, which in our example case should be `num_initial_structures * num_samples`. We can also add a more detailed assertion to check that the `key_value_pairs` stored in the database for each entry contain the correct metadata. This single, powerful test validates the entire application flow: configuration file parsing, CLI invocation, the sequential execution of all four pipeline services, and the final creation of the expected data artifact, providing a high degree of confidence in the overall system integrity.
