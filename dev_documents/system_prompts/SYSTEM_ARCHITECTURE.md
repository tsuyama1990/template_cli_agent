# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) is a sophisticated, command-line-driven software framework meticulously engineered to automate the generation of high-quality, physically-valid, and diverse training datasets for modern Machine Learning Interatomic Potentials (MLIPs), such as MACE, NequIP, and SevenNet. The core philosophy underpinning this project is the strategic removal of the human expert from the traditionally laborious and often intuition-driven data generation loop. It aims to replace this manual process with a robust, reproducible, and intelligent automated workflow. This system directly addresses a critical bottleneck in the lifecycle of MLIP development: the creation of extensive and representative datasets. Such datasets must span a wide spectrum of atomic configurations, encompassing not only stable, low-energy equilibrium states but also high-energy transition states, grain boundaries, defect structures, and other configurations that are prone to causing simulation failures. By systematically and intelligently exploring the thermodynamic phase space using a combination of simulation techniques, the tool is designed to produce datasets that cultivate MLIPs with superior accuracy, enhanced stability, and broader transferability across different physical conditions.

The pipeline is architected as a modular, sequential four-stage process: Generation, Exploration, Sampling, and Storage. Each stage is a distinct, self-contained component, ensuring maintainability and extensibility. The **Generation** stage is responsible for creating the initial seed structures that form the starting point for the simulations. It is designed to be highly flexible, supporting a wide variety of physical systems, from simple binary alloys and ionic crystals to complex multi-layered interfaces and molecular adsorption systems. A crucial aspect of this stage is the enforcement of fundamental physical constraints, such as ensuring minimum interatomic distances and appropriate simulation cell sizes, which guarantees that the initial configurations are physically plausible and ready for simulation.

The **Exploration** stage constitutes the scientific core of the system. It takes the seed structures and evolves them through time using a powerful hybrid simulation engine that combines molecular dynamics (MD) with Monte Carlo (MC) methods. This hybrid approach allows the system to not only explore the potential energy surface via the time-evolution of MD but also to perform larger, more disruptive configurational changes (like atomic swaps) via MC, which is essential for overcoming energy barriers and discovering novel structures. This stage is engineered for extreme robustness, incorporating advanced features like the automatic switching between thermodynamic ensembles (NVT vs. NPT) based on the system's topology and the dynamic integration of classical potentials (like ZBL) to prevent unphysical atomic overlaps during high-temperature simulations.

Following exploration, the **Sampling** stage intelligently curates the vast amount of data generated. Recognizing that not all simulation frames are equally valuable, this stage moves beyond simplistic random selection. It implements advanced, descriptor-based techniques like Farthest Point Sampling (FPS), which uses atomic environment descriptors (e.g., SOAP) to select a subset of configurations that is maximally diverse and non-redundant. This ensures the final dataset is information-rich and computationally efficient for training. Finally, the **Storage** stage systematically archives the curated atomic structures and their associated metadata (energy, forces, stress tensors) into a structured, queryable ASE (Atomic Simulation Environment) database. The entire workflow is governed by a flexible configuration system, enabling users to precisely define and customize complex data generation campaigns. The project will deliver both a powerful command-line interface (CLI) for experts and a user-friendly web interface (Web UI) for broader accessibility, with the ultimate ambition of accelerating the materials discovery and design cycle by democratizing the development of bespoke, high-performance MLIPs.

## 2. System Design Objectives

The primary objective of the MLIP-AutoPipe project is to deliver a fully automated, scientifically robust, and computationally efficient pipeline for the generation of MLIP training data. The design of the system is guided by a comprehensive set of goals, constraints, and success criteria to ensure it meets the needs of the computational science community.

**Goals:**
*   **Automation:** The foremost goal is to achieve end-to-end automation, minimizing the need for manual intervention. The ideal user workflow involves creating a single, declarative configuration file and executing a single command. The system must handle all subsequent steps, including parallel process management, error handling, and data organization, without user input. This reduces the cognitive load on the user and ensures reproducibility.
*   **Physical Validity:** Every structure generated by the pipeline, at every stage, must be physically meaningful. This is a non-negotiable principle. The system will enforce this through a series of validation checks: minimum interatomic distances to prevent atomic overlap, charge neutrality checks for ionic systems, and the use of simulation techniques that correctly sample thermodynamic ensembles. The goal is to produce "clean" datasets that do not require manual post-processing.
*   **Diversity and Representativeness:** The system must generate datasets that are not just large but also diverse and representative of the material's potential energy surface. The design will actively promote diversity through several mechanisms: applying strain and random perturbations in the generation phase, using high temperatures and hybrid MD/MC methods in the exploration phase, and employing advanced descriptor-based sampling techniques (FPS) in the sampling phase. This ensures the MLIP is trained on a wide variety of local environments, leading to a more robust model.
*   **Robustness and Fault Tolerance:** Scientific simulations, especially at high temperatures, can be unstable. The exploration phase, which involves many parallel simulations, must be highly resilient to failures. The system architecture will isolate individual simulation runs in separate processes. If one simulation crashes (e.g., due to a "Coulomb explosion"), it must be caught, logged for later analysis, and discarded without terminating the entire parent workflow. This ensures that a single problematic trajectory does not jeopardize a long-running data generation campaign.
*   **Modularity and Extensibility:** The system must be designed as a collection of loosely-coupled, high-cohesion modules. This is crucial for long-term maintainability and community-driven extension. It should be straightforward for a developer to add a new `StructureGenerator` for a novel material class, a new `Explorer` that uses a different simulation algorithm, or a new `Sampler` with a different selection strategy, all by inheriting from a common base class and without modifying the core pipeline orchestration logic.
*   **Usability and Accessibility:** The tool must be accessible to a wide range of users. For computational experts and for integration into larger automated workflows, a powerful and scriptable Command Line Interface (CLI) is essential. For students, experimentalists, and those less familiar with the command line, an intuitive Web User Interface (Web UI) will be provided to lower the barrier to entry. The UI will allow for interactive configuration, job monitoring, and visualization of results.

**Constraints:**
*   **Computational Resources:** The primary target environment is a multi-core workstation or a single node on a high-performance computing (HPC) cluster. The design must efficiently utilize multi-core CPUs for parallelizing the simulations. While single-GPU acceleration for MLIP models is a critical design consideration, the system will not initially be designed for multi-GPU or multi-node MPI-based parallelism, though the modular design should not preclude this in the future. Careful management of GPU memory within the multi-process environment is a key constraint to prevent resource conflicts.
*   **Dependencies:** The project will build upon the rich ecosystem of open-source Python libraries for computational science. Key dependencies will include ASE (for atomic structure representation and simulation), Pymatgen (for advanced structure generation), Hydra (for configuration management), Typer (for the CLI), and PyTorch (as the backend for MLIP models). The selection of these dependencies imposes a constraint to work within their respective APIs and data structures, but the benefit is a vast amount of pre-existing, well-tested functionality.
*   **Python Ecosystem:** The application will be developed exclusively in modern Python (3.11+), leveraging its latest features. We will enforce strict static typing using `mypy` and adhere to rigorous code quality standards using `ruff`. This constrains developers to write clean, maintainable, and well-documented code.

**Success Criteria:**
*   **Functional:** The system can successfully generate a database containing at least 1,000 unique, physically valid structures for a standard benchmark system (e.g., a Ni-Au binary alloy) from a single, complete configuration file, with the run completing in a reasonable timeframe on a standard workstation.
*   **Performance:** The Farthest Point Sampling implementation demonstrates at least a 20% higher structural diversity score (based on a defined SOAP-space metric) compared to a random sampling run on the same set of trajectories.
*   **Robustness:** The pipeline can complete a full run, successfully generating a final database, even when one of the child MD simulation processes is artificially terminated to simulate a crash.
*   **Extensibility:** A developer can add a new `StructureGenerator` for a simple cubic crystal by creating a new class that inherits from `BaseStructureGenerator` and implements the `generate` method, with no more than 5 lines of code change required in the core factory/orchestration logic.
*   **Usability:** The CLI will be successfully packaged and installable via `pip`. Its command structure (`mlip-autopipec run --config ...`) will be clear and provide useful, real-time progress feedback and logging. The Web UI will be successfully used to configure and launch a run that produces a valid, verifiable database.

## 3. System Architecture

The MLIP-AutoPipe is architected as a layered, modular system that rigorously separates concerns into three distinct layers: a user interface layer, an application/orchestration layer, and a core logic/domain services layer. This layered design is a cornerstone of the project, ensuring that the fundamental scientific logic is completely decoupled from the mechanisms through which a user interacts with the system. This promotes flexibility, maintainability, and testability.

At the highest level, the **User Interface Layer** provides the entry points for the user. This includes the Typer-based Command Line Interface (CLI) and the Streamlit-based Web User Interface (Web UI). The sole responsibility of this layer is to gather configuration parameters from the user and present the results and logs back to them. It contains no business logic.

The **Application Layer** acts as the central nervous system of the application. Its primary component is the `PipelineRunner`, a powerful orchestrator class. The `PipelineRunner`'s job is to manage the end-to-end execution of the data generation workflow. It receives a validated configuration object from the UI layer and then invokes the various services from the core logic layer in the correct sequence. It is also responsible for managing the flow of data between the stages (e.g., passing the list of generated structures to the explorer) and for handling high-level state, such as checkpointing and resuming the pipeline.

The **Core Logic Layer** contains the specialized services that perform the actual work. Each service is a self-contained module with a specific responsibility, corresponding to a stage in the pipeline:
1.  **Generation Service:** This service contains a `GeneratorFactory` which, based on the user's configuration, dynamically selects the appropriate `BaseStructureGenerator` (e.g., `AlloyGenerator`, `IonicGenerator`). The chosen generator then creates a set of initial `ase.Atoms` objects. These objects are then passed through a series of `PhysicsConstraint` validators to ensure their physical plausibility before being handed off to the next stage.
2.  **Exploration Service:** This service's `Explorer` component is the computational heart of the pipeline. It receives the list of seed structures and manages their simulation. It uses a `ProcessPoolExecutor` to distribute multiple `MDRunner` instances across available CPU cores for parallel execution. Each `MDRunner` is an independent worker responsible for a single simulation trajectory. A critical architectural decision here is the "late binding" of the MLIP calculator. The calculator object, which can be memory-intensive, is instantiated within each worker process rather than in the main process. This avoids computationally expensive serialization of large objects and provides a clean way to manage GPU resources in a multi-process environment. The `MDRunner` also encapsulates the complex logic for hybrid MD/MC and automatic ensemble switching.
3.  **Sampling Service:** The `Sampler` component receives the raw trajectory data (as file paths or in-memory objects) from the exploration stage. It employs a factory pattern to select the desired sampling strategy (e.g., `RandomSampler`, `FPSSampler`) and applies it to select a representative subset of frames.
4.  **Storage Service:** The final `Storage` service, implemented as an `AseDBWrapper`, takes the curated list of `ase.Atoms` objects from the sampler and commits them to a persistent SQLite database, along with all relevant metadata.

Data flows unidirectionally through this pipeline, with each stage producing artifacts that are consumed by the next. State between these long-running stages is managed through the filesystem (e.g., temporary XYZ trajectory files), which makes the process robust against interruptions. The final, curated state is persisted in the central ASE database.

```mermaid
graph TD
    subgraph User Interfaces
        CLI_Interface[CLI (Typer)]
        Web_UI[Web UI (Streamlit)]
    end

    subgraph Application Layer
        A[PipelineRunner Orchestrator]
    end

    subgraph Core Logic / Domain Services
        subgraph Generation
            B[GeneratorFactory]
            C[AlloyGenerator]
            D[IonicGenerator]
            E[...]
        end
        subgraph Exploration
            F[Explorer]
            G[MDRunner Worker]
            H[Hybrid MD/MC Logic]
        end
        subgraph Sampling
            I[SamplerFactory]
            J[RandomSampler]
            K[FPSSampler]
        end
        subgraph Storage
            L[AseDBWrapper]
        end
    end

    subgraph Data & State
        M[Configuration (Pydantic/YAML)]
        N[ASE Database (SQLite)]
        O[Intermediate Trajectories (.xyz)]
    end

    CLI_Interface -- Provides Config --> A
    Web_UI -- Provides Config --> A

    A -- Invokes & Passes Data --> B
    A -- Invokes & Passes Data --> F
    A -- Invokes & Passes Data --> I
    A -- Invokes & Passes Data --> L

    B --> C
    B --> D
    B --> E

    F -- Manages Parallel Runs --> G
    G -- Contains --> H

    I --> J
    I --> K

    L -- Writes To --> N
    C -- Produces --> O
    D -- Produces --> O
    E -- Produces --> O
    G -- Produces --> O
    I -- Consumes --> O


    M -- Parsed By --> CLI_Interface
    M -- Parsed By --> Web_UI
```
This architecture promotes high cohesion (components with a single, well-defined purpose) and low coupling (components that are independent of each other's internal implementation). This makes the system easier to understand, test, and extend. The strict use of Pydantic models as the data transfer objects between layers enforces clear, statically-typed contracts throughout the application.

## 4. Design Architecture

The project will be structured as a standard, modern Python package named `mlip_autopipec`. All source code will be located within the `src/` directory, making it compliant with `pip` and other standard build tools. The internal structure of the package is designed to reflect the logical layers of the system architecture, promoting clarity and ease of navigation.

**File Structure:**

```
.
├── pyproject.toml              # Project metadata, dependencies, and tool configuration
├── src
│   └── mlip_autopipec
│       ├── __init__.py
│       ├── cli.py                # Typer-based CLI application entry point
│       ├── web_ui.py             # Streamlit-based Web UI application entry point
│       ├── config
│       │   ├── __init__.py
│       │   └── models.py         # Centralized Pydantic models for all configuration
│       ├── domain
│       │   ├── __init__.py
│       │   └── models.py         # Pydantic models for core domain data (e.g., results)
│       ├── services
│       │   ├── __init__.py
│       │   ├── orchestration.py    # Contains the PipelineRunner class
│       │   ├── generators.py       # BaseGenerator ABC and all concrete implementations
│       │   ├── exploration.py      # Explorer class and MDRunner worker logic
│       │   ├── sampling.py         # BaseSampler ABC and all concrete implementations
│       │   └── storage.py          # AseDBWrapper for database interactions
│       └── shared
│           ├── __init__.py
│           ├── ase_utils.py        # Helper functions for ASE interoperability
│           └── physics.py          # Physics-based validation and analysis functions
└── tests
    ├── unit                      # Unit tests for individual components
    │   ├── test_generators.py
    │   └── ...
    └── integration               # End-to-end tests for the full pipeline
        └── test_pipeline.py
```

**Class/Function Definitions Overview:**

*   **`mlip_autopipec.cli.app` (Typer App):** The main entry point for the command-line interface. It will define commands such as `run` to execute a pipeline from a configuration file, and `ui` to launch the web interface. It will be responsible for loading YAML files and parsing them into Pydantic models before handing them off to the application layer.
*   **`mlip_autopipec.config.models.FullConfig` (Pydantic Model):** This will be the master Pydantic model that encapsulates every possible configuration option for the pipeline. It will be composed of smaller, nested models (e.g., `SystemConfig`, `ExplorationConfig`, `MCConfig`), providing a single, validated source of truth for all system parameters. This enables static analysis, autocompletion in IDEs, and robust validation of user input.
*   **`mlip_autopipec.services.orchestration.PipelineRunner`:** This is the central orchestrator class. Its primary public method, `run(config: FullConfig)`, will execute the four pipeline stages in the correct sequence. It will manage the in-memory data flow (e.g., lists of `ase.Atoms` objects) and filesystem paths between the services.
*   **`mlip_autopipec.services.generators.BaseStructureGenerator` (ABC):** An abstract base class defining the interface for all structure generators. It will enforce a single method signature, `generate(config: SystemConfig) -> list[ase.Atoms]`, ensuring that any new generator can be seamlessly integrated into the factory and the pipeline.
*   **`mlip_autopipec.services.exploration.Explorer`:** This class manages the parallel execution of simulations. Its main method will accept a list of `ase.Atoms` objects and use a `concurrent.futures.ProcessPoolExecutor` to distribute the simulation tasks to worker functions, collecting the paths to the resulting trajectory files.
*   **`mlip_autopipec.services.sampling.BaseSampler` (ABC):** This abstract base class defines the interface for all sampling strategies, with a single public method `sample(trajectories: list[Path]) -> list[ase.Atoms]`.
*   **`mlip_autopipec.services.storage.AseDBWrapper`:** This will be a context-managed class that encapsulates all interactions with the ASE SQLite database. It will provide clean, high-level methods like `write_atoms(atoms_list: list[ase.Atoms])` and abstract away the underlying database connection and transaction details.

**Data Models:**
The system's design is fundamentally schema-first, driven by Pydantic. This approach ensures data integrity and provides clear contracts between components.
*   **Configuration Models (`config/models.py`):** These models will map directly to the structure of the user-provided YAML configuration files. They will make extensive use of Pydantic's validators to enforce constraints at the point of data entry (e.g., `temperature > 0`, composition fractions must sum to 1.0). This fails fast and provides clear error messages to the user.
*   **Domain Models (`domain/models.py`):** While ASE's `Atoms` object is the primary data carrier for atomic structures, Pydantic models will be used to wrap results and metadata for structured data transfer. For instance, a `SimulationResult` model could be defined to encapsulate not just the final `Atoms` object but also metadata like the potential energy, the simulation time, and any error codes. This ensures that data passed between the services is well-defined, type-checked, and validated. This schema-first design is critical for reducing integration errors, improving developer productivity, and making the system more maintainable in the long run.

## 5. Implementation Plan

The project will be developed over two distinct, sequential cycles. This iterative approach allows us to deliver a functional core product quickly and then build upon it with more advanced capabilities.

**Cycle 1: Core Pipeline and CLI (Foundation)**

This foundational cycle is focused on building the essential backbone of the application. The primary goal is to produce a functional, end-to-end command-line tool that can perform a simplified but complete data generation workflow. The scope will be intentionally constrained to the most critical components to ensure a solid, well-tested foundation. We will implement the main `PipelineRunner` and its four constituent stages with simplified, baseline logic. The **Generation** stage will be represented by a robust `AlloyGenerator`, capable of creating random binary alloy supercells with user-specified compositions and crystal lattices. This generator will already include essential physical validation checks. The **Exploration** stage will be implemented as a pure Molecular Dynamics (MD) simulation engine, deliberately excluding the more complex Monte Carlo features for this cycle. To ensure speed and reliability for testing the core pipeline, it will use a standard, fast classical potential like EMT, which is readily available in ASE, rather than a full MLIP. The **Sampling** stage will be limited to a simple but effective `RandomSampler`, which will select a specified number of frames from the generated MD trajectories. Finally, the **Storage** component will be fully implemented, providing a robust wrapper around the ASE database to ensure that the final, curated data can be correctly and persistently saved. The sole user interface in this cycle will be a Typer-based CLI, which will accept a YAML configuration file to drive the entire pipeline. The comprehensive Pydantic configuration models will be fully defined in this cycle to establish a stable and extensible configuration schema for the entire project's lifecycle. All code will be fully type-hinted, documented with docstrings, and accompanied by a comprehensive suite of unit tests for each component. A crucial deliverable is a full end-to-end integration test that verifies the entire pipeline execution for a simple case.

**Cycle 2: Advanced Features and Web UI**

Building directly on the stable foundation of Cycle 1, this cycle will introduce the advanced, scientifically sophisticated features that provide the core value proposition of MLIP-AutoPipe, alongside a user-friendly Web UI to enhance accessibility. The **Generator** service will be significantly expanded to support a wider range of physical systems, with the implementation of an `IonicGenerator` (enforcing charge neutrality) and an `InterfaceGenerator` (for creating heterostructures and surfaces). The **Explorer** service will be the focus of the most intensive development. It will be upgraded to a full hybrid MD/MC engine, incorporating Monte Carlo moves such as atom swaps and vacancy hops to enable more efficient exploration of complex energy landscapes. A major milestone will be the integration of support for external MLIP models (e.g., MACE), which will include the critical feature of mixing the MLIP with the ZBL classical potential to ensure physical realism at short interatomic distances and prevent simulation failures. The explorer will also be enhanced with "intelligent" features like automatic ensemble switching. The **Sampling** stage will be upgraded with the implementation of the `FPSSampler`, which requires integrating a third-party library for computing SOAP structural descriptors and implementing the farthest-point sampling algorithm itself. The final major deliverable will be a Web UI, likely developed using Streamlit, which will provide an interactive and intuitive way for users to build configuration files, launch pipeline runs (by invoking the CLI/core services in a subprocess), monitor their progress, and visualize the resulting structures directly from the output database. The CLI will also be updated to expose all the new configuration parameters available in this cycle. Testing will be significantly more complex, requiring mock MLIP models and dedicated tests to verify the statistical correctness of the MC moves and the diversity metrics of the FPS algorithm.

## 6. Test Strategy

Our testing strategy is designed to be comprehensive, ensuring the reliability, numerical correctness, and physical validity of the generated data. The strategy is multi-layered, combining unit, integration, and user acceptance testing, and is designed to be fully automated within a Continuous Integration (CI) framework.

**Cycle 1 Test Strategy:**
*   **Unit Testing:** Each component will be rigorously tested in isolation to verify its specific functionality.
    *   **Generators:** Tests for the `AlloyGenerator` will programmatically verify that the output structures have the correct atomic composition, lattice type, and total number of atoms for a given supercell size. Critically, the physical constraint validators located in the `shared` module will be tested with deliberately "illegal" structures (e.g., `ase.Atoms` objects with overlapping atoms) to ensure they raise the expected exceptions. Filesystem interactions will be mocked to keep the tests fast and isolated.
    *   **Storage:** The `AseDBWrapper` will be tested against a temporary database file created using `pytest`'s `tmp_path` fixture. These tests will confirm that the wrapper can correctly create a new database, connect to it, write single and multiple `Atoms` objects, and that the data read back from the database is bit-for-bit identical to the original data.
    *   **Configuration:** The Pydantic configuration models are a critical part of the system's robustness. They will be tested by attempting to instantiate them with various forms of invalid data (e.g., negative temperatures, compositions that do not sum to 1.0, invalid element symbols) and asserting that a `pydantic.ValidationError` is raised with a user-friendly error message.
    *   **CLI:** The Typer-based CLI will be tested in isolation using the `typer.testing.CliRunner`. These tests will simulate command-line invocations, checking that the `run` command correctly invokes the `PipelineRunner` with a configuration object parsed from a test YAML file. Error conditions, such as a path to a non-existent file, will also be tested to ensure the CLI exits gracefully with a non-zero status code and a helpful error message.
*   **Integration Testing:** A single, crucial end-to-end integration test will be created to verify the entire pipeline. This test will use `typer.testing.CliRunner` to invoke the `cli.py run` command, pointing it to a minimal but complete configuration file for a simple binary alloy (e.g., CuAu). To ensure the test runs quickly, it will use the fast EMT potential and specify a very small number of simulation steps. The test will assert that the pipeline completes successfully (exit code 0) and that the final ASE database file is created and contains the exact expected number of atomic structures (`num_initial_structures * num_samples`). This test provides confidence that all the components are correctly wired together.

**Cycle 2 Test Strategy:**
*   **Unit Testing:**
    *   **Advanced Generators:** The new generators will have dedicated unit tests. For the `IonicGenerator`, tests will assert that all generated structures are charge-neutral. For the `InterfaceGenerator`, tests will verify the correct creation of slabs and interface geometries.
    *   **Explorer:** Testing the hybrid MD/MC logic is of paramount importance. We will design tests that run very short simulations on simple, well-defined systems (e.g., a 2-atom box). One test will attempt an atom swap move and assert that the atomic numbers of the selected atoms are correctly exchanged and that the move is accepted or rejected according to the Metropolis criterion for a pre-calculated energy change. Mocking the underlying MLIP calculator will be essential to isolate and test the logic of the `Explorer`'s new features.
    *   **Sampling:** The `FPSSampler` will be tested by providing it with a small, known set of structures with clear similarity differences. The test will assert that the sampler correctly and deterministically selects the most diverse structures from the set, validating the correctness of the FPS algorithm's implementation.
*   **Integration Testing:** The integration test suite will be expanded to cover the new, more complex workflows. A new test will be added that runs the full pipeline with a mock MLIP model and the hybrid MD/MC explorer enabled. The assertions for this test will be more sophisticated; for example, it will analyze the output database to confirm a wider range of compositional order parameters, demonstrating that the swap moves were effective. Another new integration test will specifically target the FPS sampling. It will run two full pipeline executions on the same input, one with random sampling and one with FPS. The test will then analyze both output databases, calculate a diversity metric for each, and assert that the score for the FPS-generated database is significantly higher.
*   **User Acceptance Testing (UAT):**
    *   **Web UI:** The Web UI requires its own dedicated testing. This will be done using an automated end-to-end testing framework like Playwright. The test script will automate the process of a user opening the UI, filling in the configuration form, clicking the "Run" button, and then verifying that the results displayed in the UI (e.g., the number of structures in a table) match the contents of the generated database.
    *   **Jupyter Notebooks:** For both cycles, UAT will be facilitated by providing users with well-documented Jupyter Notebooks. These notebooks will serve as both a test artifact and a tutorial, guiding the user on how to programmatically load the generated database, visualize the structures, and plot property distributions (e.g., energy histograms, radial distribution functions). This provides a transparent and verifiable way for users to confirm the quality of the generated data.
