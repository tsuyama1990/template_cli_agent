# SYSTEM ARCHITECTURE: MLIP-AutoPipe

## 1. Summary

MLIP-AutoPipe is a highly automated and robust framework designed to generate high-quality, physically valid, and diverse datasets for training modern Machine Learning Interatomic Potentials (MLIPs), such as MACE and SevenNet. The core philosophy of the project is to "remove the human expert from the loop" by creating a system that intelligently explores the thermodynamic phase space of a given material system to produce training data. Traditional methods of data generation often rely on manual, time-consuming, and error-prone processes, requiring deep domain expertise to select appropriate atomic configurations. MLIP-AutoPipe aims to replace this with a systematic, reproducible, and efficient pipeline.

The system is engineered to handle a wide variety of physical systems, including multi-component alloys, ionic crystals, covalent materials, interfaces between different materials, and surface adsorption phenomena. It also includes a knowledge-based approach that can generate initial structures from simple chemical formulas by leveraging crystallographic databases and symmetry principles. The primary output is a curated database of atomic structures, complete with associated metadata, ready to be consumed by MLIP training codes.

The architecture is built as a sequential pipeline, comprising four distinct stages: Generation, Exploration, Sampling, and Storage. The **Generation** stage creates initial seed structures based on user-defined parameters. Crucially, it incorporates physics-based validation checks to ensure that these initial structures are physically plausible, for instance by preventing atoms from being too close to one another. The **Exploration** stage, which is the most complex part of the system, uses Molecular Dynamics (MD) and hybrid Monte Carlo (MC) simulations to evolve these seed structures. This allows the system to discover a wide range of atomic configurations, including high-energy states and transition states, which are often critical for training robust MLIPs. The **Sampling** stage then intelligently selects a subset of the most informative structures from the vast number generated during exploration, using techniques like Farthest Point Sampling (FPS) to maximize diversity. Finally, the **Storage** stage saves the selected structures and their properties into a structured, queryable ASE database.

A key design principle is robustness and fault tolerance. The system is designed for long-running simulations, with features like progressive saving to prevent data loss in case of a crash. It also incorporates sophisticated safety mechanisms, such as preventing "Coulomb explosions" in ionic systems and automatically handling structures that violate physical constraints. Furthermore, the system is highly configurable through the Hydra framework, allowing for flexible customisation of every pipeline stage. It is intended to be used primarily as a Command-Line Interface (CLI) for batch processing and integration into larger automated workflows, but will also feature a secondary Web UI for interactive use, configuration, and visualisation.

## 2. System Design Objectives

The design of MLIP-AutoPipe is guided by a set of clear objectives aimed at creating a powerful, flexible, and user-friendly tool for materials scientists and MLIP developers.

**Primary Objectives:**
- **Automation:** The foremost objective is to completely automate the process of generating MLIP training data. The system should require minimal user intervention, taking high-level inputs (e.g., chemical composition) and producing a complete, analysis-ready dataset.
- **Data Quality and Physical Validity:** The generated structures must be physically realistic. The system must enforce physical constraints at every stage, such as minimum interatomic distances and correct charge balancing in ionic systems. The goal is to produce data that accurately represents the potential energy surface of the material.
- **Diversity of Structures:** A good training set must be diverse, covering a wide range of configurations, temperatures, and pressures. The system is designed to explore the potential energy landscape broadly, capturing not just stable, low-energy states but also the less common, high-energy configurations that are often the "difficult" cases for an MLIP to learn. The integration of hybrid MD/MC methods and advanced sampling algorithms like FPS is central to achieving this objective.
- **Robustness and Reliability:** The exploration phase can involve thousands of independent, long-running simulations. The system must be robust to failures. This includes graceful error handling, automatic recovery, prevention of common simulation artifacts (like Coulomb explosions), and checkpointing to ensure that progress is not lost.

**Technical Objectives:**
- **Modularity and Extensibility:** The architecture is designed to be modular. Each stage of the pipeline (Generation, Exploration, Sampling) and each component within those stages (e.g., different structure generators, different simulation engines) are implemented as interchangeable modules. This makes the system easy to extend with new generators, simulation techniques, or sampling methods in the future.
- **Performance and Scalability:** The framework must be able to efficiently leverage modern computing hardware. It is designed for parallel execution, using Python's multiprocessing capabilities to run many simulations simultaneously. It also includes optimisations to handle large ML models efficiently, such as a "late binding" approach for the MLIP calculators to avoid costly serialisation overhead.
- **Usability and Configuration:** While the primary interface is a powerful CLI for experts and automated workflows, the system must also be accessible. Configuration is managed via the Hydra framework, providing a clear and flexible way to define and override parameters. A secondary Web UI will further lower the barrier to entry for interactive exploration and visualisation.

**Success Criteria:**
- The system will be considered successful if it can reliably generate a training dataset that, when used to train a reputable MLIP model (e.g., MACE), results in a potential with predictive accuracy comparable to or exceeding that of a model trained on a manually curated dataset.
- The system should demonstrate the ability to generate valid structures for at least three distinct classes of materials (e.g., alloy, ionic, and adsorption systems) without requiring code changes, only configuration updates.
- The pipeline should successfully complete end-to-end runs for systems containing at least 100 atoms, demonstrating its stability and performance for non-trivial simulations.
- The modular design will be validated by the successful addition of a new structure generator or sampling algorithm with minimal changes to the core framework.

## 3. System Architecture

The MLIP-AutoPipe framework is designed as a modular, four-stage pipeline orchestrated by a central `PipelineRunner`. Data flows sequentially from one stage to the next, with outputs from each stage being persisted to disk and a central database to ensure state isolation and fault tolerance.

```mermaid
graph TD
    A[Start] --> B{1. Generation};
    B --> C{2. Exploration};
    C --> D{3. Sampling};
    D --> E{4. Storage};
    E --> F[End];

    subgraph "User Input"
        G[Hydra Config]
    end

    subgraph "Pipeline Stages"
        B; C; D; E;
    end

    subgraph "Data Persistence"
        H[(ASE Database)];
        I[/file/system/xyz];
    end

    G --> B;
    B -- Seed Structures --> I;
    B -- Records --> H;
    I -- Input for Exploration --> C;
    C -- Trajectories --> I;
    I -- Input for Sampling --> D;
    D -- Selected Structures --> I;
    I -- Final Structures --> E;
    E -- Final Records & Metadata --> H;

```

**Component Breakdown:**

1.  **Configuration (Hydra):** The entire workflow is controlled by a set of YAML configuration files managed by Hydra. This allows users to define the physical system, select the desired components (e.g., `AlloyGenerator`, `MDEngine`), and tune all necessary parameters (temperature, pressure, etc.) without modifying the source code.

2.  **`PipelineRunner`:** This is the central orchestrator. It reads the configuration, initialises the components for each stage, and executes them in the correct order. It manages the overall state of the workflow, including handling restarts by checking for the existence of intermediate files.

3.  **Stage 1: Generation:**
    *   A `GeneratorFactory` reads the configuration and instantiates the appropriate generator class (e.g., `AlloyGenerator`, `IonicGenerator`) which inherits from a common `BaseGenerator`.
    *   The selected generator creates a set of initial "seed" structures.
    *   These structures are passed through a series of `PhysicsConstraints` validators (e.g., `overlap_check`, `ensure_supercell_size`) to guarantee physical plausibility.
    *   The validated structures are saved to disk as `.xyz` files and their creation is logged in the ASE database.

4.  **Stage 2: Exploration:**
    *   The `ExplorerEngine` (e.g., `MDEngine`) reads the seed structures from the previous stage.
    *   It uses a `ProcessPoolExecutor` to run multiple simulations in parallel, distributing the structures among worker processes.
    *   Each worker process performs a simulation (e.g., MD or hybrid MD/MC) on its assigned structure. A key feature is the "late binding" of the MLIP calculator, which is instantiated inside the worker process to avoid serialisation issues and improve efficiency.
    *   The engine includes advanced logic for automatic ensemble switching (NVT vs. NPT) based on vacuum detection, and can mix the MLIP with a classical potential like ZBL for improved stability at high temperatures.
    *   Trajectories (the sequence of atomic positions over time) are progressively saved to disk.

5.  **Stage 3: Sampling:**
    *   The `Sampler` component processes the raw trajectory files generated during exploration.
    *   Based on the user's configuration, it applies a sampling strategy. This can be simple `RandomSampling` or the more sophisticated `FarthestPointSampling` (FPS), which uses SOAP descriptors to select a structurally diverse subset of frames.
    *   The output is a smaller, curated set of `.xyz` files representing the most informative structures.

6.  **Stage 4: Storage:**
    *   The `Storage` component takes the final sampled structures.
    *   It iterates through them, calculates final properties if needed, and writes each one as a new row in the central ASE (SQLite) database.
    *   Crucially, it stores rich metadata alongside each structure, such as its origin, the simulation parameters used, and calculated properties like energy and forces.

**Data Flow:**
The primary data objects are ASE `Atoms` objects, which represent the atomic structures. These are passed between stages primarily via the file system (as `.xyz` or trajectory files). The ASE database acts as a central ledger and checkpointing mechanism, tracking the state and history of all structures throughout the pipeline. This file-based handoff between major stages ensures that they are decoupled and that the entire process can be stopped and resumed.

## 4. Design Architecture

The software design emphasizes modularity, separation of concerns, and adherence to object-oriented principles. The file structure reflects this clear separation.

**File Structure (High-Level):**
```
src/mlip_autopipec/
├── cli/
│   ├── __init__.py
│   └── main.py              # Main CLI entrypoint (using Click)
├──- core/
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
│   ├── generation/          # Structure generation modules
│   │   ├── __init__.py
│   │   ├── alloy.py
│   │   └── base.py
│   ├── exploration/         # Exploration (MD/MC) modules
│   │   ├── __init__.py
│   │   └── md_engine.py
│   ├── sampling/            # Sampling modules
│   │   ├── __init__.py
│   │   ├── fps.py
│   │   └── random.py
│   └── storage/             # Storage (database) modules
│       ├── __init__.py
│       └── ase_db_wrapper.py
├── utils/
│   ├── __init__.py
│   └── physics.py           # Physics validation and utility functions
└── gui/
    └── main_gui.py          # Web UI entrypoint
```

**Key Design Patterns:**

*   **Dependency Inversion & Interfaces:** The `core` modules depend on abstract interfaces defined in `core/interfaces.py` (e.g., `IGenerator`, `IExplorer`). The concrete implementations of these services reside in the `services/` directory. This decouples the core orchestration logic from the specific implementation details, making the system extensible.
*   **Factory Pattern:** `core/factories.py` contains factories (e.g., `GeneratorFactory`) responsible for instantiating the correct concrete class from the `services/` directory based on the user's configuration. This centralises the object creation logic.
*   **Pydantic for Configuration and Data Models:** All configuration will be defined using strict Pydantic models in `domain/configuration.py`. This provides automatic type checking, validation, and clear documentation for all settings. Internal data structures, such as the results from a simulation, will also be represented by Pydantic models in `domain/data_models.py` to ensure type safety throughout the application.
*   **Separation of Concerns:**
    *   **`cli/`:** Handles only command-line argument parsing and user interaction, delegating all business logic to the `core`.
    *   **`core/`:** Contains the high-level orchestration logic and component interfaces, but no specific implementation details.
    *   **`services/`:** Contains all the business logic for the specific tasks (generation, exploration, etc.). These modules are self-contained and do not depend on each other directly.
    *   **`domain/`:** Purely defines the data structures and configuration schema for the entire application. It contains no business logic.
    *   **`utils/`:** Holds stateless, reusable utility functions, particularly those related to physical constraints and calculations.

This layered architecture ensures that the system is maintainable, testable, and extensible. Each component has a clearly defined responsibility, reducing complexity and making it easier to modify or replace parts of the system without affecting others.

## 5. Implementation Plan

The project will be developed over two distinct cycles. Cycle 1 focuses on establishing the core command-line pipeline and fundamental features. Cycle 2 builds upon this foundation to add advanced exploration capabilities, broader material support, and a user-friendly web interface.

**CYCLE01: Core Pipeline and CLI Foundation**

This initial cycle is focused on creating a fully functional, end-to-end command-line tool that can perform the entire four-stage pipeline for a single, well-defined class of materials: alloys. The primary goal is to establish a robust architectural foundation, including the CLI, the pipeline orchestrator, the database wrapper, and the configuration system. The exploration engine will use a basic Molecular Dynamics simulation, and sampling will be limited to a simple random selection.

**Key Features for CYCLE01:**
1.  **CLI Development:** Implement the main entrypoint using `click`, with commands to run the full pipeline.
2.  **Configuration:** Set up the Hydra configuration structure with Pydantic models for type-safe settings.
3.  **Pipeline Orchestrator:** Build the central `PipelineOrchestrator` that manages the four stages (Generation, Exploration, Sampling, Storage).
4.  **Database Wrapper:** Create a robust wrapper for the ASE database to handle all storage operations.
5.  **Alloy Generator:** Implement the first structure generator for creating multi-component alloy structures.
6.  **Basic MD Explorer:** Implement a basic MD exploration engine using an MLIP calculator (e.g., MACE). It will run NVT simulations at a fixed temperature.
7.  **Random Sampler:** Implement a simple sampler that randomly selects frames from the MD trajectories.
8.  **Unit and Integration Tests:** Develop a comprehensive suite of tests covering the core components to ensure reliability.

**CYCLE02: Advanced Exploration and User Interface**

This cycle enhances the capabilities of the core pipeline significantly, introducing more sophisticated physics and making the tool more versatile and accessible. The major focus is on the hybrid MD/MC exploration engine, which is critical for generating diverse and challenging structures. Support for more material types will be added, and a web-based GUI will be developed to provide an alternative, interactive way to use the tool.

**Key Features for CYCLE02:**
1.  **Hybrid MD/MC Engine:** Upgrade the explorer to a hybrid engine that combines MD with MC moves (e.g., atom swaps), which is crucial for efficient phase space exploration in alloys.
2.  **Advanced Generators:** Implement additional generators for other material classes, such as ionic crystals and surface adsorption systems, leveraging the factory pattern established in Cycle 1.
3.  **Automatic Ensemble Switching:** Implement the logic to automatically detect vacuum layers and switch between NPT (for bulk) and NVT (for slabs) ensembles during exploration.
4.  **Farthest Point Sampling (FPS):** Implement the FPS algorithm as a more intelligent sampling strategy to maximize the structural diversity of the final dataset.
5.  **Web UI:** Develop a graphical user interface using a web framework (e.g., Streamlit or Flask) that allows users to configure and run the pipeline from a browser and visualise the generated structures.
6.  **ZBL Potential Integration:** Add the capability to mix the MLIP with a ZBL potential to handle short-range atomic interactions correctly, preventing simulation failures at high temperatures.
7.  **Expanded Test Suite:** Add new tests to cover the advanced features, including the hybrid engine, new generators, and the FPS sampler.

## 6. Test Strategy

The project will employ a comprehensive, multi-layered testing strategy to ensure code quality, correctness, and robustness. This includes unit tests, integration tests, and user acceptance tests (UAT), with a strong emphasis on continuous integration.

**CYCLE01 Test Strategy:**
The focus in the first cycle is on building a solid foundation of tests for the core components.
*   **Unit Testing:** Each class and function in the `services`, `domain`, and `utils` modules will have dedicated unit tests. We will use `pytest` as the testing framework and `pytest-mock` to isolate components. For example, the `AlloyGenerator` will be tested to ensure it produces structures with the correct composition and cell size, without needing a real database or explorer. The `PipelineOrchestrator` will be tested by mocking the individual stage components to verify that it calls them in the correct order and handles file I/O properly. Pydantic models will be tested to ensure validation rules are correctly enforced.
*   **Integration Testing:** We will create a suite of integration tests that verify the interaction between components. A key integration test will be an end-to-end run of the CLI on a small, well-defined alloy system (e.g., CuAu). This test will use a fast, simple potential (like ASE's built-in EMT potential) to run the entire pipeline: generating a few structures, running a very short MD simulation, sampling a handful of frames, and storing them in a temporary database. The test will then assert that the final database contains the expected number of structures and that their properties are physically plausible. This will validate that the data flow between stages via the file system and database is working correctly.

**CYCLE02 Test Strategy:**
The strategy for Cycle 2 expands the test suite to cover the new, more complex features.
*   **Unit Testing:** New unit tests will be created for the advanced components. The `HybridMDMC_Engine` will be tested to ensure it correctly alternates between MD steps and MC moves. The logic for `FarthestPointSampling` will be tested with a known set of structures to verify that it selects the most diverse subset. The new generators (e.g., `IonicGenerator`) will have their own dedicated tests to check for specific constraints like charge neutrality. The `detect_vacuum` utility will be tested with various cell configurations to ensure it correctly identifies bulk vs. slab systems.
*   **Integration Testing:** The end-to-end integration tests will be expanded into a more comprehensive test suite. We will create specific tests for each major new feature. For instance, a dedicated test will run the pipeline with the hybrid MD/MC engine enabled and assert that atom swaps have actually occurred. Another test will run the pipeline with the FPS sampler and verify that the final sampled structures are indeed more diverse (e.g., by measuring the standard deviation of interatomic distances) than those from a random sample. We will also add tests for the new generators, running the full pipeline for a small ionic system, for example. The Web UI will be tested manually during development and could later be automated with tools like Playwright if needed.

**Overall Test Infrastructure:**
*   **Continuous Integration (CI):** All tests will be run automatically on every commit and pull request using a CI service (e.g., GitHub Actions).
*   **Code Coverage:** We will monitor code coverage using `pytest-cov` to ensure that all critical parts of the codebase are covered by tests.
*   **Static Analysis:** `ruff` will be used for linting and formatting to maintain a consistent and high-quality codebase.
