# CYCLE 01: SPECIFICATION

## 1. Summary

This document provides the detailed technical specification for the foundational first cycle of the MLIP-AutoPipe project. The overarching goal of Cycle 01 is to deliver a functional, end-to-end, non-interactive pipeline that serves as the architectural backbone for the entire system. This cycle is about building the "engine" of the car before adding the "self-driving" capabilities. By the conclusion of this cycle, the project will have a command-line executable capable of taking a minimal user input—essentially just the chemical composition of a material—and autonomously executing the entire workflow of initial structure generation, robust DFT calculation, data persistence, and baseline model training. This will result in the creation of a "Generation 0" Machine Learning Interatomic Potential (MLIP). The core achievement of this cycle is the successful automation of the most labor-intensive and error-prone parts of the traditional MLIP creation process.

The scope of this cycle is deliberately focused on creating a solid, reliable, and well-engineered foundation. We will establish the complete project structure using modern Python tooling (`uv` and `pyproject.toml`). A cornerstone of this cycle is the implementation of a schema-first design, where all data structures for configuration and scientific data are rigorously defined and validated using Pydantic. This ensures type safety and contract-driven development from day one. The central `WorkflowOrchestrator` will be implemented to manage the linear, sequential flow of this initial pipeline. The key scientific modules to be fully implemented are the `StructureGenerator` (Module A), the `LabelingEngine` (Module C) specifically for Quantum Espresso, and the `TrainingEngine` (Module D) for the ACE model. A crucial piece of "embedded intelligence" in this cycle is the `ConfigExpander`, the component responsible for translating the user's simple request into a comprehensive and physically sound set of parameters for the DFT calculations and other components, effectively encapsulating a significant amount of expert knowledge.

It is critical to note what is *not* included in this cycle. The advanced features that give the system its full autonomy—namely the surrogate-based exploration with MACE (Module B) and the on-the-fly active learning loop (Module E)—are explicitly deferred to Cycle 02. The pipeline in this cycle is a "one-shot" process; it generates a static dataset and trains a single model. However, the successful completion of this cycle is non-negotiable for the project's success. It will provide the first tangible demonstration of the core value proposition: the robust automation of complex, multi-step computational materials science workflows. It also allows us to validate our architectural choices, the integrity of our data flow, the reliability of our database interactions, and the correctness of our interface with the external DFT code. This foundational work will de-risk the subsequent development and provide the stable platform upon which the more advanced, intelligent features will be built.

## 2. System Architecture

The system architecture for Cycle 01 is designed to be a robust, modular, and scalable foundation for the entire project. The blueprint focuses on implementing the essential components required for a linear, end-to-end workflow. The file structure is explicitly laid out to enforce a clear separation of concerns, isolating configuration, data models, interfaces, and concrete implementations. All files and directories that will be created or significantly modified during this cycle are highlighted in **bold**. This structure is not temporary; it is the final intended structure, with certain components being implemented as stubs to be filled in during the next cycle. This approach ensures that the architectural vision is established from the very beginning.

**File Structure Blueprint:**

```
mlip_autopipec/
├── **pyproject.toml**  # Project metadata, dependencies, and build configuration
├── uv.lock
└── src/
    └── mlip_autopipec/
        ├── **__init__.py**
        ├── **cli.py**           # User-facing command-line interface entry point
        ├── **config/**          # Package for configuration models and loading
        │   ├── **__init__.py**
        │   ├── **models.py**      # **Defines MinimalConfig and FullConfig Pydantic models**
        │   └── **loader.py**      # **Implements the ConfigExpander heuristic engine**
        ├── **data/**            # Package for scientific data models
        │   ├── **__init__.py**
        │   └── **models.py**      # **Defines AtomsModel and DFTResult Pydantic models**
        ├── **database/**        # Package for abstracting database interactions
        │   ├── **__init__.py**
        │   └── **ase_db_wrapper.py** # **Concrete implementation of the database wrapper**
        ├── **workflows/**       # Package for high-level orchestration
        │   ├── **__init__.py**
        │   └── **orchestrator.py**  # **Manages the linear pipeline logic for Cycle 01**
        ├── **interfaces/**      # Package defining the abstract contracts for components
        │   ├── **__init__.py**
        │   └── **engines.py**     # **Defines ABCs: IStructureGenerator, ILabelingEngine, etc.**
        └── **engines/**         # Package for concrete implementations of the engine interfaces
            ├── **__init__.py**
            ├── **structure_generator.py** # **Module A: SQS, NMS, etc. implementation**
            ├── explorer_sampler.py        # Stub file for Module B, not implemented in C01
            ├── **labeling_engine.py**     # **Module C: Quantum Espresso implementation**
            ├── **training_engine.py**     # **Module D: ACE model training implementation**
            └── simulation_engine.py       # Stub file for Module E, not implemented in C01
```

**Architectural Data Flow for Cycle 01:**

The workflow is strictly sequential and orchestrated from a single point of control, ensuring a clear and traceable process.

1.  **Invocation**: The process starts when the user runs the application from the command line, invoking a function within `cli.py`.
2.  **Configuration Expansion**: The CLI instantiates the `ConfigExpander` from `config/loader.py`. This component reads the user-provided `input.yaml`, parses it into the `MinimalConfig` Pydantic model, and then applies a series of complex, physics-based heuristics to generate the comprehensive, fully-parameterized, and validated `FullConfig` object.
3.  **Orchestrator Initialization**: The `WorkflowOrchestrator` is instantiated, receiving the `FullConfig` object and concrete instances of the required engines (`StructureGenerator`, `QuantumEspressoEngine`, `TrainingEngine`) as dependencies. This use of dependency injection is a key architectural choice for ensuring testability.
4.  **Initial Structure Generation**: The orchestrator makes its first call to `structure_generator.generate_structures()`. This engine (Module A) uses the information in `FullConfig` to determine the material type and generates a list of `AtomsModel` objects.
5.  **Data Persistence (Write)**: The orchestrator iterates through the newly generated list of `AtomsModel` objects and passes each one to the `AseDBWrapper`, which saves them to the SQLite database, marking them with a status of "uncalculated".
6.  **Data Retrieval**: The orchestrator then queries the `AseDBWrapper` to retrieve the batch of "uncalculated" structures. This step ensures that the pipeline is data-driven and can potentially be resumed if interrupted.
7.  **DFT Labeling**: The orchestrator iterates through the retrieved structures, passing each one to the `QuantumEspressoEngine.run_calculation()` method (Module C). This is the most computationally intensive step. The engine encapsulates all the logic for creating input files, running the `pw.x` subprocess, parsing the potentially complex output, and robustly handling any errors. Upon success, it returns a validated `DFTResult` object.
8.  **Data Persistence (Update)**: For each successfully calculated structure, the orchestrator calls the `AseDBWrapper` again to update the corresponding database record with the new `DFTResult` data and change its status to "calculated".
9.  **Dataset Aggregation**: Once all structures in the batch have been labeled, the orchestrator queries the database one last time to retrieve the complete, labeled dataset (a list of `AtomsModel` and their corresponding `DFTResult` objects).
10. **Model Training**: This complete dataset is then passed to the `TrainingEngine.train_model()` method (Module D). This engine handles the final step of fitting the ACE model to the data and saving the final potential file (e.g., `final_model.ace`) to the filesystem. The path to this file is the final output of the pipeline.

## 3. Design Architecture

The design of the system's components for Cycle 01 is rigorously defined by a schema-first philosophy, where Pydantic models serve as the unambiguous contracts for all data that flows between different parts of the application. This approach is critical for ensuring data integrity, providing runtime validation, and creating a system that is largely self-documenting. The architecture is designed to be highly decoupled, relying on abstract interfaces to minimize dependencies between the high-level workflow logic and the low-level implementation details of the scientific engines.

**`config.models.py` - The Configuration Contract:**
This module is the bedrock of the "expert in a box" concept. It will define the two-tier configuration schema that drives the entire pipeline.
-   `MinimalConfig`: This Pydantic `BaseModel` defines the user-facing API. It will be extremely sparse, containing only a nested `System` model with two required fields: `elements: list[str]` and `composition: str`. Every other parameter will be optional (`Optional[...]`). This design minimizes the barrier to entry for the user.
-   `FullConfig`: This is the canonical, internal configuration object. It will be a large, deeply nested Pydantic model that contains a validated, type-safe parameter for every conceivable option in the workflow. It will be composed of several sub-models, such as `DFTCompute`, `TrainingParams`, and `StructureGenerationParams`. For example, `DFTCompute` will contain fields like `ecutwfc: float = Field(..., gt=0)`, `smearing: Literal["gaussian", "mv"]`, and a nested `Magnetism` model. Crucially, the `FullConfig` model will be configured with `frozen=True` (or the Pydantic v2 equivalent), making it immutable after creation. This prevents components from modifying the configuration state mid-workflow, ensuring consistency and reproducibility. The `ConfigExpander` is the designated factory for this object.

**`data.models.py` - The Scientific Data Contract:**
This module defines the canonical representations for the scientific data that is the currency of the pipeline. This decouples the internal logic from the specific object implementations of external libraries like ASE.
-   `AtomsModel`: This model will represent an atomic structure. It will contain fields such as `symbols: list[str]`, `positions: list[list[float]]`, `cell: list[list[float]]`, and `pbc: tuple[bool, bool, bool]`. It will employ Pydantic's powerful validation features. For instance, a `@model_validator` will be used to ensure that the length of the `positions` list is identical to the length of the `symbols` list, and that the shape of the `cell` matrix is 3x3. This prevents malformed structural data from ever entering the system.
-   `DFTResult`: This model will represent the output of a quantum mechanical calculation. It will contain fields like `energy: float`, `forces: list[list[float]]`, and `stress: list[float]`. It will also have its own validators. A key validator will ensure that the shape of the `forces` array is consistent with the `AtomsModel` to which the result belongs. Another will enforce that the `stress` is always a 6-element Voigt vector, as required by many simulation engines.

**`interfaces.engines.py` - The Abstract Engine Contracts:**
This module is the key to the system's modularity and testability. It defines the abstract "roles" that concrete classes will fulfill, using Python's `abc` module.
-   `IStructureGenerator(ABC)`: Will define the contract for all structure generation methods. It will have one abstract method: `@abstractmethod def generate_structures(self, config: FullConfig) -> list[AtomsModel]: pass`.
-   `ILabelingEngine(ABC)`: Will define the contract for any DFT calculation engine. Its core method will be `@abstractmethod def run_calculation(self, atoms: AtomsModel, config: FullConfig) -> DFTResult: pass`.
-   `ITrainingEngine(ABC)`: Will define the contract for model training. Its method will be `@abstractmethod def train_model(self, dataset: list[tuple[AtomsModel, DFTResult]], config: FullConfig) -> Path: pass`.

**Producers and Consumers of Data:**
This strict, schema-driven design creates a clear and verifiable network of data producers and consumers. The `ConfigExpander` is the sole producer of the `FullConfig` object; all other components are consumers. The `StructureGenerator` is the initial producer of `AtomsModel` objects. The `LabelingEngine` is a consumer of `AtomsModel` and a producer of `DFTResult`. The `AseDBWrapper` acts as both a consumer and producer of these models to the persistence layer. The `TrainingEngine` is a consumer of the complete dataset. This clear flow, enforced by Pydantic's runtime validation, is essential for building a reliable and predictable system.

## 4. Implementation Approach

The implementation of Cycle 01 will be a methodical, step-by-step process that builds the system from the ground up, starting with the foundational data contracts and moving progressively towards the high-level orchestration logic. This approach ensures that each new component is built upon a stable and well-tested foundation.

1.  **Project Initialization**: The first step is to establish the development environment. We will run `uv init` to create the `pyproject.toml` file. This file will be immediately populated with the project's metadata, command-line entry points, and the list of core dependencies for this cycle, which includes `pydantic`, `ase`, `typer`, `numpy`, and the chosen ACE training library. The basic directory structure (`src/mlip_autopipec/...`) will also be created.
2.  **Schema-First Development**: The first coding task is to implement the Pydantic models that define the system's data contracts. We will create `config/models.py` and `data/models.py`. This step involves translating the design architecture into concrete, typed, and validated `BaseModel` classes for `MinimalConfig`, `FullConfig`, `AtomsModel`, and `DFTResult`. This is done first because these models are the glue that holds the entire system together.
3.  **Configuration Logic**: With the schemas defined, the `ConfigExpander` will be implemented in `config/loader.py`. This is a critical piece of logic that involves researching and codifying the heuristics for DFT calculations. For example, it will contain a dictionary mapping elements to their recommended pseudopotentials and energy cutoffs from the SSSP library.
4.  **Database Abstraction**: The `AseDBWrapper` class will be implemented in `database/ase_db_wrapper.py`. This class will encapsulate all interactions with the `ase.db` module. Methods will be created for `initialize_db`, `add_structures`, `update_structure_with_result`, `get_structures_by_status`, and `get_full_dataset`. Each method will handle the translation between our internal Pydantic models (`AtomsModel`, `DFTResult`) and the format required by the ASE database.
5.  **Interface Definition**: The abstract base classes (`IStructureGenerator`, `ILabelingEngine`, `ITrainingEngine`) will be formally defined in `interfaces/engines.py`. This step solidifies the API contracts that the concrete engines must adhere to.
6.  **Core Engine Implementation - Labeling**: The implementation of the `QuantumEspressoEngine` is the most complex task in this cycle and will be tackled with care. It will involve:
    *   A method for dynamically generating Quantum Espresso input file strings, populating templates with parameters from the `FullConfig` object.
    *   A method that uses Python's `subprocess` module to execute `pw.x`. This method must be robust, correctly handling command-line arguments, capturing stdout/stderr, and managing timeouts.
    *   A detailed and resilient parser for the Quantum Espresso output format, using regular expressions to find and extract the key physical quantities.
    *   The implementation of the multi-stage error recovery state machine, which will be triggered if the output parser detects a convergence failure.
7.  **Core Engine Implementation - Others**: Following the labeling engine, the `StructureGenerator` and `TrainingEngine` will be implemented. The `StructureGenerator` will contain the logic for the different generation methods (SQS, etc.), dispatching to the correct one based on the configuration. The `TrainingEngine` will be a wrapper around the chosen ACE library, responsible for formatting the data correctly and calling the library's training function.
8.  **Orchestration Logic**: With all the components built, the `WorkflowOrchestrator` will be implemented in `workflows/orchestrator.py`. Its `run_linear_pipeline` method will contain the high-level logic that calls each engine in the correct sequence, as detailed in the architectural flow. It will manage the state transitions of the data in the database (e.g., from "uncalculated" to "calculated").
9.  **User Interface**: Finally, the `cli.py` file will be created. Using the Typer library, a simple CLI application will be built. It will have a single command, `run`, which will take the path to an `input.yaml` file as an argument. The command's function will be to instantiate all the necessary objects (`ConfigExpander`, `WorkflowOrchestrator`, etc.) and start the pipeline, printing informative status messages to the console as the process unfolds.

## 5. Test Strategy

The testing strategy for Cycle 01 is designed to build a pyramid of confidence. It starts with a wide base of granular unit tests to ensure each individual component is correct, followed by a layer of integration tests to ensure the components communicate and transfer data correctly through the core pipeline.

**Unit Testing Approach (Min 300 words):**
The foundation of our testing strategy is comprehensive unit tests for every component, ensuring they are correct in isolation. This will be achieved by extensively using mocking to replace external dependencies and other system components.
-   **Configuration & Data Models**: The `ConfigExpander` is a primary candidate for unit testing. We will create a test suite with a wide variety of sample `input.yaml` files as input. These tests will cover different material types (alloys, molecules), different user specifications (e.g., requesting magnetism), and edge cases. For each input, we will run the expander and perform detailed assertions on the resulting `FullConfig` object, verifying that the heuristics have produced the expected, physically correct parameters. The Pydantic models themselves will also be tested. We will write tests that attempt to instantiate them with malformed data (e.g., a `forces` array of the wrong shape for a given `AtomsModel`) and assert that a `pydantic.ValidationError` is raised, confirming that our data contracts are being enforced.
-   **Engines**: Each engine will be tested in complete isolation. The `QuantumEspressoEngine` will receive the most thorough testing. Its dependency on the external `pw.x` executable will be mocked using `unittest.mock.patch`. We will create a test data directory containing a collection of saved, real-world Quantum Espresso output files. This collection will include examples of successful calculations, different convergence failures, crashes, and other error states. The tests will feed the content of these files to the engine's output parser. We will then assert that for successful runs, the physical quantities are parsed with high precision, and for failed runs, the correct exception type or error code is returned. We will also write specific tests to verify that the multi-stage error recovery logic is triggered in the correct sequence when a failure is detected.
-   **Database and Orchestrator**: The `AseDBWrapper` will be tested by having it connect to a temporary, in-memory SQLite database, ensuring that tests are fast and don't have side effects. We will test every method, verifying that data is written and read back correctly, and that our Pydantic models are correctly serialized and deserialized. The `WorkflowOrchestrator`'s logic will be tested by providing it with mock instances of all the engines. We can then verify the correctness of its high-level logic by asserting that the mock engines' methods are called in the correct order and with the expected arguments.

**Integration Testing Approach (Min 300 words):**
While unit tests verify the components, integration tests verify the connections between them. For Cycle 01, the focus is on ensuring the integrity of the linear data pipeline from end to end.
-   **Core Pipeline Data Flow**: The most critical integration test will validate the main success path. This test will use a *real* `ConfigExpander`, a *real* `StructureGenerator`, a *real* `AseDBWrapper` (connected to a temporary file-based database), and a *real* `WorkflowOrchestrator`. However, the `LabelingEngine` will be *mocked*. This is a crucial strategic decision, as it allows us to test the entire data flow without the extreme cost and potential flakiness of running real DFT calculations. The mock `LabelingEngine` will be programmed to return deterministic, predictable `DFTResult` objects. The test will then assert the following sequence:
    1. The CLI is invoked with a test configuration.
    2. The `StructureGenerator` creates `AtomsModel` objects, which are then written to the database. We will query the DB to confirm their existence.
    3. The `Orchestrator` calls the mock `LabelingEngine`. We will verify it was called with the correct `AtomsModel`.
    4. The predictable `DFTResult` from the mock is correctly used to update the structures in the database. We will query the DB again to confirm the update.
    5. The final, complete dataset is correctly retrieved and passed to the `TrainingEngine` (which can remain a lightweight mock for this test).
-   **Configuration Propagation Test**: A second important integration test will verify that the parameters from the `FullConfig` object are correctly propagated through the orchestrator and into the engines. For example, we will create a configuration with a specific, non-default energy cutoff. The test will use a real orchestrator but a mock `LabelingEngine`. We will then use mocking tools to "spy" on the `run_calculation` method of the mock engine and assert that the `config` object it receives contains the exact, correct energy cutoff value that was specified in the initial input. This confirms that the configuration is being correctly passed down through the layers of the application. Together, these tests will provide strong confidence in the robustness and correctness of the entire foundational pipeline.
