# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe system is an advanced, automated framework designed to address a critical bottleneck in the development of modern Machine Learning Interatomic Potentials (MLIPs): the generation of high-quality training data. The core purpose of this project is to create a robust and intelligent pipeline that can produce physically valid, diverse, and comprehensive datasets with minimal human intervention. This directly confronts the "human expert in the loop" problem, where the quality of an MLIP is often limited by the intuition, time, and manual effort of a domain expert in selecting and calculating training structures. Traditional methods often rely on a patchwork of scripts, manual selection of structures from known crystallographic databases, or simple random perturbations of equilibrium lattices. These approaches are not only labour-intensive but, more importantly, often fail to capture the full complexity of a material's potential energy surface. The resulting datasets can be sparse in critical regions, such as those corresponding to high-energy configurations, defect migration pathways, or phase transitions. This sparsity leads to MLIPs that, while accurate for configurations similar to their training data, are notoriously unreliable when used in simulations that explore new regions of the phase space, limiting their predictive power and scientific utility.

MLIP-AutoPipe systematically overcomes these challenges by employing a sophisticated, multi-stage pipeline grounded in the principles of statistical mechanics and materials science. The system's philosophy is to replace manual intuition with algorithmic rigour. Instead of merely generating random atomic configurations, the system intelligently explores the thermodynamically accessible phase space of a given material using proven simulation techniques, primarily Molecular Dynamics (MD) and Monte Carlo (MC) methods. This dynamic approach ensures that the generated dataset is naturally biased towards physically relevant configurations. It includes not only stable, low-energy structures but also the thermally-activated, high-energy structures and transition states that are frequently encountered in real-world simulations of material behaviour at finite temperatures. These "hard cases" are precisely the data points that are most crucial for training a robust and accurate MLIP that can generalize well to unseen configurations, leading to more reliable and transferable models.

The system is designed from the ground up for versatility, supporting a wide range of physical systems out-of-the-box. This includes multi-component alloys, where it can explore chemical ordering; ionic crystals with specific charge-balancing constraints; covalent materials; complex interfaces between different materials; and surface adsorption systems. Furthermore, it incorporates a knowledge-based approach that can leverage crystallographic databases and chemical formula information to generate plausible initial structures, providing a sensible starting point for the exploration. The architecture is modular, built around a clear, four-step pipeline: (1) **Generation** of a diverse set of initial seed structures; (2) **Exploration** of the potential energy surface through simulated dynamics to generate vast trajectories of configurations; (3) intelligent **Sampling** of the most informative and unique structures from these trajectories to create a compact and efficient dataset; and (4) final **Storage** of the curated dataset in a structured, queryable database format. Key engineering decisions, such as the use of a "late binding" calculator to manage computational resources efficiently in a parallel environment and the implementation of rigorous, physics-based validation checks at every stage, elevate MLIP-AutoPipe from a simple script collection to a production-ready, enterprise-grade scientific tool. Its successful implementation will significantly lower the barrier to entry for developing high-quality bespoke MLIPs, accelerating the pace of materials discovery and design.

## 2. System Design Objectives

The design and development of MLIP-AutoPipe are guided by a clear set of objectives, constraints, and measurable success criteria. These are essential for ensuring that the final product is not only functional but also robust, usable, and maintainable, meeting the high standards of scientific computing. The objectives are carefully defined to address the specific challenges of MLIP data generation and to deliver a tool that is both powerful and practical for its target audience of computational scientists and engineers.

**Primary Goals:**

*   **Full Automation:** The paramount goal is to create a "fire-and-forget" system. A user should be able to define a material system in a high-level configuration file, execute a single command, and receive a comprehensive training dataset without further manual input. This objective aims to democratise MLIP development by abstracting away the complex, low-level details of setting up and running dozens or hundreds of individual simulations. The system should manage the entire workflow, including error handling and recovery, making it suitable for running on unattended high-performance computing clusters.
*   **Data Diversity and Quality:** The system must generate datasets that are maximally informative. "Diversity" is not merely about randomness; it means producing a wide variety of structures that span a broad range of energies, forces, and stresses. This includes equilibrium states, strained configurations, structures with defects (like vacancies or interstitials), amorphous phases, and configurations near phase transitions. "Quality" is ensured by enforcing strict physical validity checks at all stages, guaranteeing that the final dataset is free from unphysical artefacts that could poison the MLIP training process.
*   **Physical Realism:** All generated structures must be physically plausible. The system will incorporate a multi-layered defence against unphysical configurations. This includes fundamental checks to prevent overlapping atoms, verification of correct coordination numbers where applicable, and enforcement of charge neutrality constraints in ionic systems. Crucially, the use of simulation techniques like MD and MC, which inherently follow the laws of physics, ensures that the generated structures represent states that are accessible under real thermodynamic conditions.
*   **Robustness and Fault Tolerance:** Long-running scientific simulations are prone to failure due to numerical instabilities, hardware issues, or incorrect input parameters. The system must be designed to be resilient. It will handle common errors (e.g., simulation crashes or "Coulomb explosions" in ionic systems) gracefully at the level of a single simulation without causing the entire pipeline to fail. It will feature robust checkpointing between each of the four main pipeline stages, allowing for the easy resumption of interrupted workflows and saving valuable computational time.
*   **Modularity and Extensibility:** The architecture must be modular to facilitate future expansion and community contributions. It should be straightforward for other developers to add new structure generation algorithms (e.g., for polymers or molecular crystals), alternative exploration techniques (e.g., simulated annealing, genetic algorithms, basin hopping), or more sophisticated sampling methods (e.g., clustering-based approaches). This will be achieved through a plugin-style architecture based on abstract base classes.
*   **Performance and Scalability:** The pipeline, particularly the exploration stage, must be computationally efficient. This includes native support for the parallel execution of simulations across multiple CPU cores to reduce the time-to-solution for large-scale data generation campaigns. The design must also be mindful of memory usage, especially when dealing with large MLIP models and long simulation trajectories.

**Key Constraints:**

*   **Technology Stack:** The system will be developed exclusively in Python 3.12+, leveraging the mature and powerful scientific Python ecosystem. Key libraries will include the Atomic Simulation Environment (ASE) for all atomistic manipulations, Pydantic for rigorous data modelling and configuration management, Hydra for flexible and composable configuration, and Typer for the command-line interface. This ensures that the tool is accessible to a wide audience of scientists who are already familiar with this stack.
*   **Development Methodology:** The project will be developed following a strict two-cycle plan. Cycle 1 will focus on building the core infrastructure and a minimal viable pipeline to validate the architecture. Cycle 2 will implement the advanced, computationally intensive exploration and sampling features. This incremental approach mitigates risk and allows for early feedback.
*   **Interface Definition:** The primary user interface will be a command-line tool. This is essential for automation, scripting, and integration into larger high-throughput computational workflows common in materials science research. A secondary, simple web-based GUI will be provided for interactive setup, visualization, and educational purposes.

**Success Criteria:**

*   **Functional Verification:** The system must be able to complete an end-to-end run, taking a configuration file for a standard binary alloy (e.g., Fe-Pt) and producing a final ASE database containing at least 1,000 unique, physically valid atomic structures, without any manual intervention after the initial command.
*   **Quality Verification:** An MLIP model (e.g., MACE) trained *exclusively* on a dataset generated by MLIP-AutoPipe must demonstrate a statistically significant improvement in prediction accuracy (defined as a minimum 15% reduction in the root-mean-square error on atomic forces for a held-out test set) compared to a model of the same architecture trained on a dataset of equivalent size containing only randomly perturbed and volumetrically strained conventional crystal structures.
*   **Performance Verification:** The pipeline, when running the exploration step on a standard multi-core workstation (e.g., 8 cores), should demonstrate effective parallel speed-up. Running with 4 parallel workers should be at least 3 times faster than running with a single worker, showcasing efficient resource utilization.

## 3. System Architecture

The MLIP-AutoPipe system is designed as a modular, pipeline-based application, adhering to the principle of separation of concerns. This architecture promotes maintainability, testability, and extensibility by ensuring that each component has a single, clearly defined responsibility. The overall data flow is orchestrated by a central `PipelineRunner`, which acts as a state machine, executing a sequence of four main stages: Generation, Exploration, Sampling, and Storage. This design provides a clear and logical progression for the data, from initial concept to final curated dataset.

**Core Components:**

1.  **User Interfaces (CLI & Web UI):** These are the primary entry points for the user. They are thin layers responsible for capturing user intent and delegating to the core pipeline logic.
    *   **CLI (Command-Line Interface):** Built with Typer and Rich, this is the primary interface for power users and automated, non-interactive workflows. It accepts a path to a configuration file, provides clear feedback on the pipeline's progress, and is designed to be easily scriptable.
    *   **Web UI:** A secondary interface, likely built with Streamlit or Gradio, for interactive use cases. Its purpose is to lower the barrier to entry, allowing users to visually configure a run, start the pipeline, and inspect the generated structures without leaving the browser.

2.  **Configuration Manager:** Utilising the combined power of Hydra and Pydantic, this component is responsible for managing all runtime settings. Hydra allows for a flexible and composable configuration system using YAML files, enabling users to easily mix and match different settings. Pydantic provides a rigorous schema validation layer on top of this, ensuring that all configuration parameters are of the correct type and within valid ranges before any computation begins. This catches user errors early and makes the system more robust.

3.  **Pipeline Orchestrator (`PipelineRunner`):** This is the heart of the application. It initialises all other components based on the validated configuration and manages the main data processing workflow. It ensures that each stage is executed in the correct order: Generation -> Exploration -> Sampling -> Storage. It is also responsible for state management, such as handling the flow of data (lists of atomic structures) between stages and implementing the checkpointing mechanism that saves intermediate results to disk, allowing the pipeline to be resumed after an interruption.

4.  **Structure Generators:** This component is a factory responsible for creating the initial "seed" structures that form the starting point for the exploration phase. It is designed for extensibility, consisting of an abstract base class (`BaseStructureGenerator`) that defines a common interface (a `generate` method). Multiple concrete implementations for different material types (e.g., `AlloyGenerator`, `IonicGenerator`, `InterfaceGenerator`) inherit from this base class. This polymorphic design allows the `PipelineRunner` to work with any type of generator, and new generators can be added without modifying the core orchestration logic.

5.  **Exploration Engine:** This is the primary computational engine of the pipeline, where the bulk of the computational effort is expended. It takes the seed structures and uses MD and/or hybrid MD/MC simulations to explore the potential energy surface, generating a vast number of diverse configurations in the form of simulation trajectories. It contains sophisticated logic for automatically selecting the appropriate simulation ensemble (NPT for bulk systems vs. NVT for surfaces) based on a physical analysis of the system's geometry. It also integrates auxiliary potentials (like the Ziegler-Biersack-Littmark potential, ZBL) to prevent unphysical simulation behaviour like atomic fusion at high temperatures.

6.  **Sampling Module:** After the exploration stage has produced a massive number of candidate structures (potentially millions), this module is responsible for intelligently selecting a small, information-rich subset. It implements various strategies, from simple random selection (as a baseline) to advanced techniques like Farthest Point Sampling (FPS). FPS uses structural descriptors (like the Smooth Overlap of Atomic Positions, SOAP) to select a subset of structures that are maximally different from one another, ensuring the final dataset is diverse and not redundant.

7.  **Storage Manager (`AseDBWrapper`):** This component handles all interactions with the output database. It provides a clean, high-level API for writing the final curated structures and their associated metadata (e.g., energy, forces, originating simulation parameters) into an ASE-compatible SQLite database. By abstracting the database logic into a dedicated component, the rest of the application is shielded from the low-level details of database connections and queries.

8.  **Physics Validator:** This is a cross-cutting utility that provides a library of functions for checking the physical plausibility of atomic structures. It is used as a quality gate by multiple components throughout the pipeline. For example, the Generators use it to ensure their initial structures do not have overlapping atoms, and the Exploration Engine uses it to detect and handle simulation crashes.

**Data Flow Diagram:**

```mermaid
graph TD
    subgraph User Input
        A[Config File .yaml] --> B{Configuration Manager};
    end

    subgraph Pipeline Core
        B -- Validated Pydantic Models --> C[Pipeline Orchestrator];
        C --> D[1. Structure Generators];
        D -- Seed Structures --> E[2. Exploration Engine];
        E -- Raw Trajectories --> F[3. Sampling Module];
        F -- Curated Structures --> G[4. Storage Manager];
    end

    subgraph Data & Utilities
        H[Physics Validator] -- Enforces Physical Constraints --> D;
        H -- Enforces Physical Constraints --> E;
        I[ASE Database .db] <--> G;
    end

    style User Input fill:#cde4ff
    style Pipeline Core fill:#d5e8d4
    style Data & Utilities fill:#fff2cc
```

## 4. Design Architecture

The software design of MLIP-AutoPipe is founded on modern object-oriented principles, a schema-first development approach, and a commitment to the SOLID principles of software design. This ensures the codebase is robust, easy to understand, maintainable, and extensible. Pydantic data models are used pervasively to define clear, unambiguous data contracts between all components, eliminating a common source of bugs in complex data processing pipelines.

**File Structure:**

The source code will be organised within a `src/mlip_autopipec` directory. This structure makes it an installable Python package, which simplifies dependency management and deployment. The modular layout separates concerns logically, making it easy for developers to locate relevant code.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py                 # Typer-based Command-Line Interface application
├── web_ui.py              # Web-based GUI application (e.g., Streamlit)
│
├── core/
│   ├── __init__.py
│   ├── orchestrator.py    # Contains the main PipelineRunner class
│   └── physics.py         # Advanced physics validation and utility functions (e.g., vacuum detection)
│
├── config/
│   ├── __init__.py
│   └── models.py          # All Pydantic models for configuration, defining the schema
│
├── generators/
│   ├── __init__.py
│   ├── base.py            # Abstract BaseStructureGenerator class, defining the interface
│   └── alloy.py           # Concrete implementation for creating alloy structures
│   └── ionic.py           # Concrete implementation for ionic materials
│
├── exploration/
│   ├── __init__.py
│   └── engine.py          # The MD/MC Exploration Engine, a computationally intensive component
│
├── sampling/
│   ├── __init__.py
│   └── samplers.py        # Implementation of sampling algorithms (e.g., Random, FPS)
│
└── storage/
    ├── __init__.py
    └── database.py        # AseDBWrapper for abstracting database interactions
```

**Key Class Definitions:**

*   **`config.models.FullConfig(pydantic.BaseModel)`:** The top-level Pydantic model that serves as the root of the entire configuration tree. It contains nested models for each major component, such as `SystemConfig`, `ExplorationConfig`, and `SamplingConfig`. Its primary responsibility is to provide a single, validated object representing the user's desired workflow.
*   **`core.orchestrator.PipelineRunner`:** The main orchestrator class. Its primary public method, `run()`, will execute the full Generation -> Exploration -> Sampling -> Storage sequence. It will be initialised with a `FullConfig` object, from which it will configure and instantiate all the necessary sub-components (generators, explorers, etc.). It holds the state of the pipeline as it executes.
*   **`generators.base.BaseStructureGenerator(abc.ABC)`:** An abstract base class that defines the interface for all structure generators. It will have a single abstract method, `generate() -> list[ase.Atoms]`. This use of an abstract base class is a key part of the extensible design, allowing new generator types to be added without any changes to the core `PipelineRunner`.
*   **`generators.alloy.AlloyGenerator(BaseStructureGenerator)`:** A concrete implementation for creating alloy structures based on parameters in the `SystemConfig`. It will contain the logic for setting up a lattice, populating it with atoms, and ensuring the final composition is correct.
*   **`exploration.engine.MDExplorer`:** This class will encapsulate the complex logic for running simulations. It will have a primary method like `explore(structures: list[ase.Atoms]) -> list[ase.Atoms]` that takes a list of seed structures and returns a much larger list of structures from the resulting trajectories. It will internally manage the ASE `Calculator` and dynamics objects, hiding these low-level details from the orchestrator.
*   **`storage.database.AseDBWrapper`:** This class will provide a simple, high-level interface for database operations. Crucially, it will be implemented as a context manager, so it can be used with a `with` statement to automatically handle the opening and closing of the database connection, preventing resource leaks. Its primary method will be `write_structures(structures: list[ase.Atoms])`.

**Data-Centric Design with Pydantic:**

The system's design is heavily reliant on Pydantic to act as the single source of truth for all configuration and data structures. This schema-first approach provides several profound advantages that contribute to the overall quality of the software:
*   **Proactive Data Validation:** Configuration files provided by the user are automatically and rigorously validated against the defined Pydantic models upon loading. This catches a wide range of common errors (typos in keys, incorrect data types, values outside of valid ranges) immediately, providing the user with clear, actionable feedback before any expensive computations are started. For instance, a validator on the `SystemConfig` model will ensure that the elemental compositions sum to 1.0.
*   **Clear and Enforced Contracts:** The Pydantic models serve as explicit, machine-readable documentation for the data that flows between different components of the system. A developer can look at a function's signature and know exactly the shape and type of the data it expects, because it will be defined by a Pydantic model. This eliminates ambiguity and prevents entire classes of bugs related to mismatched data assumptions.
*   **Improved Developer Experience:** Using typed models provides excellent autocompletion and static type-checking support in modern IDEs like VSCode or PyCharm. This improves developer productivity, reduces the cognitive load required to understand the codebase, and makes refactoring safer and easier. This is particularly important for ensuring the long-term maintainability of a complex scientific code.

## 5. Implementation Plan

The development of the MLIP-AutoPipe project will be divided into two distinct, sequential cycles. This phased approach allows for the incremental delivery of functionality, mitigating risk and allowing for a focused development effort in each phase. The first cycle establishes a solid foundation and a minimal viable product, while the second cycle builds the advanced, computationally intensive features that provide the core scientific value.

**Cycle 1: Core Pipeline, CLI, and Foundational Components**

*   **Objective:** The primary goal of this cycle is to build a minimal-viable, end-to-end working pipeline. This version will be fully functional but will deliberately defer the complex simulation and sampling logic by using simple placeholder components. The main purpose is to establish the project structure, define all the data models, and implement the core orchestration logic, thereby proving the viability of the overall architecture. By the end of this cycle, a user will be able to successfully generate a simple dataset for an alloy system using the command line, which serves as a crucial first milestone.
*   **Key Features & Detailed Implementation Steps:**
    1.  **Project Scaffolding:** Initialise the project with the defined file structure under `src/mlip_autopipec`. This includes creating all the necessary directories and `__init__.py` files. Configure the `pyproject.toml` file to recognise the new package, allowing it to be installed in editable mode for development.
    2.  **Pydantic Data Models:** Implement the complete set of Pydantic models for configuration in `config/models.py`. This is the first and most important step of the implementation. It includes `SystemConfig`, `ExplorationConfig`, `SamplingConfig`, and the top-level `FullConfig`. These models will include detailed field definitions, type hints, and validation logic from the very start to enforce data integrity.
    3.  **CLI Implementation:** Develop the main command-line interface in `cli.py` using Typer. The CLI will have a main `run` command that accepts a required `--config` option pointing to a YAML configuration file. This module will be responsible for loading the configuration, instantiating the `PipelineRunner`, triggering the run, and providing user-friendly progress updates using the Rich library.
    4.  **Database Wrapper:** Implement the `AseDBWrapper` in `storage/database.py`. This class will handle the creation of and connection to the output SQLite database. It will be implemented as a context manager and provide a simple, robust method to write a list of `ase.Atoms` objects to the database.
    5.  **Core Orchestrator:** Implement the `PipelineRunner` in `core/orchestrator.py`. Its `run` method will load the configuration and execute the four pipeline stages in the correct sequence. It will be responsible for passing the data (the list of `Atoms` objects) from one stage to the next.
    6.  **Basic Generator:** Implement the `BaseStructureGenerator` abstract class in `generators/base.py`. Then, implement a simple, functional `AlloyGenerator` in `generators/alloy.py`. This generator will be capable of creating a specified number of random solid-solution alloy structures in a given lattice, including a physical validity check to prevent atom overlaps.
    7.  **Placeholder Components:** For this cycle, the Exploration and Sampling stages will be implemented as simple pass-through components to complete the pipeline structure. The `MDExplorer` will have a method that simply returns the input list of structures without modification. The `RandomSampler` will be implemented with its final interface but will perform a simple random selection. This strategy is crucial for enabling early end-to-end testing without tackling the full complexity at once.
    8.  **Unit Testing:** Develop a comprehensive suite of unit tests for all the components created. This includes tests for the Pydantic models' validation logic, the CLI command parsing using a test runner, the database wrapper's ability to read and write, and the alloy generator's output.

**Cycle 2: Advanced Exploration Engine and Intelligent Sampling**

*   **Objective:** With the core pipeline and architecture now validated, this cycle will focus on implementing the sophisticated computational heart of the application. The goal is to replace the placeholder components from Cycle 1 with the full-featured, physics-based exploration and sampling engines. This will transform the tool from a simple structure generator into an intelligent, autonomous data generation framework capable of producing high-quality datasets.
*   **Key Features & Detailed Implementation Steps:**
    1.  **Advanced Exploration Engine:** Implement the full `MDExplorer` in `exploration/engine.py`. This is the most significant and complex task of the entire project.
        *   Integrate MLIP models (e.g., MACE, loaded via ASE's calculator interface) for calculating forces and energies during the simulation.
        *   Implement the logic for running MD simulations for a specified number of steps using ASE's dynamics modules (e.g., `Langevin` for NVT/NPT dynamics).
        *   Add support for hybrid MD/MC moves, specifically the atom swap (`SwapMove`), which is essential for efficiently exploring chemical order in alloys.
        *   Implement the ZBL potential mixing logic. This will likely require a custom ASE calculator that combines the forces and energies from the MLIP and the ZBL potential.
        *   Develop the vacuum detection algorithm in `core/physics.py` and use it within the `MDExplorer` to automatically and dynamically switch between NVT and NPT simulation ensembles, a key feature for robustly handling different types of material systems.
    2.  **Intelligent Sampling:** Implement the Farthest Point Sampling (FPS) algorithm in `sampling/samplers.py`. This will involve integrating a library (e.g., `dscribe`) to compute the SOAP (Smooth Overlap of Atomic Positions) structural descriptors, which will serve as the feature vectors for the FPS algorithm. The FPS algorithm itself will then be implemented to select a diverse subset of structures from the exploration trajectories.
    3.  **Additional Generators:** To broaden the applicability of the tool, implement the remaining structure generators, such as the `IonicGenerator`. This generator will include specific logic to handle oxidation states and ensure that the generated ionic crystal structures are charge-neutral.
    4.  **Web User Interface:** Develop a simple, proof-of-concept Web UI in `web_ui.py` using Streamlit. This UI will provide a graphical way for users to define run parameters, start the pipeline, and view visualizations of the output structures, making the tool more accessible to new users.
    5.  **Integration Testing:** Create a suite of integration tests that run the full pipeline on small, well-defined test cases. These tests are critical for ensuring that all the complex components work together correctly. They will use a fast, simple potential (like ASE's built-in EMT potential) to allow for execution within a reasonable time in a CI/CD environment. The tests will verify the entire workflow—from generation through exploration and sampling to storage—produces a valid and expected result.

## 6. Test Strategy

A comprehensive, multi-layered testing strategy is essential to ensure the correctness, robustness, and reliability of the MLIP-AutoPipe software. The strategy is designed to catch bugs at different levels, from isolated component logic to the interaction between all parts of the system. It is divided by cycle to align with the incremental development plan.

**Cycle 1 Test Strategy:**

The focus of testing in the first cycle will be on rigorous **unit tests**. Since the foundational components are being built from the ground up, ensuring each one works correctly in isolation is paramount before they are integrated. This approach helps to pinpoint the source of errors quickly and efficiently.

*   **Unit Testing Approach:**
    *   **Pydantic Models (`test_config.py`):** The validation logic in the configuration models is a critical line of defence against user error. Tests will be written to verify this logic thoroughly. This includes positive tests that load a valid YAML configuration and assert that the `FullConfig` object is created with the correct attributes. More importantly, it includes negative tests that use `pytest.raises(ValidationError)` to ensure that invalid configurations (e.g., compositions not summing to 1.0, negative temperatures, invalid lattice types) are correctly rejected with informative error messages.
    *   **Database Wrapper (`test_storage.py`):** The `AseDBWrapper` will be tested against a temporary, in-memory SQLite database to avoid disk I/O. Tests will cover the full lifecycle: creating a new database file, writing a list of known `ase.Atoms` objects, reading them back using the ASE API, and verifying that the data (atomic numbers, positions, cell) remains intact and has not been corrupted.
    *   **Alloy Generator (`test_generator.py`):** The generator will be tested to confirm that it produces structures that adhere to the user's specification. Tests will run the generator and assert that: (1) the number of returned `ase.Atoms` objects matches the `num_structures` parameter; (2) each structure has the correct number and ratio of atomic species as defined in the `composition`; and (3) all generated structures pass the physical validation checks, specifically that the minimum distance between any two atoms is greater than a reasonable physical cutoff.
    *   **Pipeline Orchestrator (`test_orchestrator.py`):** The `PipelineRunner`'s logic will be tested using mock objects from Python's `unittest.mock` library. The generator, explorer, sampler, and storage components will be replaced with mocks. The test will then call the `run` method and assert that each component's primary method was called exactly once, in the correct sequence, and that the data object (the list of atoms) was passed between them as expected. This verifies the orchestration logic without needing the actual components.
    *   **CLI (`test_cli.py`):** The command-line interface will be tested using `click.testing.CliRunner` (or its equivalent for Typer). These tests will simulate command-line invocations. A key test will invoke the `run` command with a valid configuration file while mocking the `PipelineRunner` to assert that the CLI successfully parses the arguments and calls the orchestrator. Another test will invoke the CLI with a path to a non-existent configuration file and assert that the application exits with a non-zero status code and prints a user-friendly "File not found" error message.

**Cycle 2 Test Strategy:**

With the core components unit-tested, the testing focus in Cycle 2 shifts towards **integration testing**, while still adding dedicated unit tests for the new, complex algorithms. The primary goal is to verify that the fully-featured pipeline works correctly from end to end.

*   **Unit Testing Approach (New Components):**
    *   **MD Explorer (`test_exploration.py`):** Testing the full `MDExplorer` is challenging. Unit tests will focus on its complex setup and decision-making logic, without running expensive simulations. For example, by providing known "bulk" and "slab" structures and mocking the `detect_vacuum` function, we can test that the explorer correctly selects the NPT and NVT ensembles, respectively. We will also mock the ASE calculator and dynamics objects to assert that they are initialised with the correct parameters (temperature, pressure, timestep) from the configuration.
    *   **Physics Utilities (`test_physics.py`):** The `detect_vacuum` algorithm will be unit-tested with a set of pre-defined `ase.Atoms` objects representing known cases (bulk FCC crystal, a slab with a vacuum layer, a single molecule) to ensure it correctly classifies each one.
    *   **FPS Sampler (`test_samplers.py`):** The FPS sampler's core algorithm will be tested with a small, deterministic set of input feature vectors to verify that it reproducibly selects the expected subset of points that are farthest apart in the feature space.

*   **Integration Testing Approach (`test_integration.py`):**
    *   **Mini-Pipeline Runs:** A suite of integration tests will be developed to run the entire pipeline on small, fast-to-simulate systems. A crucial choice here is to use ASE's built-in EMT (Effective Medium Theory) potential, which is computationally inexpensive and deterministic, allowing these tests to run efficiently within a CI/CD environment.
    *   **MD Workflow Verification:** A core test case will run a full workflow for a small binary alloy system (e.g., a 16-atom Cu-Au cell). The test will configure a short MD run and assert that: (1) the pipeline completes successfully without errors; (2) the final database is created and contains the expected number of structures; and (3) the average displacement of atoms between the initial and final structures is non-zero, proving the MD simulation ran.
    *   **Hybrid MC Workflow Verification:** A specific test will be designed to verify the `SwapMove`. It will run a short MD/MC simulation on a binary system and assert that the chemical ordering of the atoms in the final structure is different from the initial one, confirming that atom swaps occurred.
    *   **Failure Mode Testing:** An integration test will be designed to intentionally trigger a failure (e.g., by creating a starting structure with overlapping atoms that should cause the MD simulation to crash). The test will verify that the application handles this error gracefully within the parallel worker, logs the error correctly, and completes the pipeline for the other, valid structures, demonstrating the system's robustness.

*   **User Acceptance Testing (UAT):** For the final user-facing features, particularly the Web UI, UAT scenarios will be defined in the `UAT.md` documents. These will involve manually testing the Web UI and running the CLI with more complex, realistic example configuration files to confirm that the output is scientifically useful, correct, and meets user expectations.
