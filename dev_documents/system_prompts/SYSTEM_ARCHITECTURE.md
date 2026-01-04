# SYSTEM_ARCHITECTURE.md

## 1. Summary

The MLIP-AutoPipe project represents a paradigm shift in the generation of training data for Machine Learning Interatomic Potentials (MLIPs). Its core philosophy is to transition the process from a manually intensive, expert-driven art to a fully automated, reproducible scientific workflow. In the domain of computational materials science, the accuracy and generalisability of an MLIP are critically dependent on the quality and diversity of the training dataset. Historically, these datasets were compiled from existing databases, which often lack coverage of high-energy or out-of-equilibrium configurations, or were generated through simple perturbations of known ground-state structures. This approach frequently leads to "brittle" potentials that perform admirably within their training distribution but fail catastrophically when used in simulations that explore novel thermodynamic regimes, such as phase transitions, defect nucleation, or chemical reactions at surfaces. These failures are a significant bottleneck in the predictive power of computational simulations. MLIP-AutoPipe is architected to systematically address this fundamental limitation.

The framework's primary objective is to **remove the human expert from the loop**, not by simplifying the physics, but by encapsulating expert-level physical and algorithmic knowledge into a robust, automated pipeline. It intelligently explores the vast potential energy surface of a given chemical system to discover not just stable configurations, but also the physically plausible, high-energy structures that are essential for building robust models. By leveraging proven simulation techniques like Molecular Dynamics (MD) and Monte Carlo (MC), the system can simulate the effects of temperature and pressure, driving the material into challenging configurations that an MLIP must learn to handle correctly. These "hard" examples—such as structures near a phase transition boundary, strained crystal lattices, amorphous phases, and transition states for diffusion—are precisely the data points required to train MLIPs that can be trusted in predictive simulations. The system is explicitly designed to find the "edge cases" of materials science, which are often the most scientifically interesting.

The system's architecture is built upon a modular, four-stage pipeline that ensures a clear separation of concerns and promotes extensibility. The first stage, **Generation**, is responsible for creating a diverse set of initial seed structures. This goes far beyond simple crystal lattice generation, accommodating a wide array of material types, including multi-component alloys, ionic crystals with charge-neutrality constraints, complex interfaces between different materials, and molecular adsorption on surfaces. The second stage, **Exploration**, is the computational heart of the framework. It takes the seed structures and subjects them to simulated thermodynamic conditions, using a sophisticated hybrid MD/MC engine to evolve their geometries. This engine includes critical, physics-aware features such as automatic switching between thermodynamic ensembles (NPT for bulk, NVT for surfaces) to prevent simulation artifacts, and the integration of the classical ZBL potential to correctly model the strong repulsive forces at short interatomic distances, thereby preventing simulation crashes at high temperatures. The third stage, **Sampling**, intelligently curates the vast amount of data generated during exploration. Instead of simple random selection, it employs advanced techniques like Farthest Point Sampling (FPS) to select a final dataset that is maximally diverse and minimally redundant, ensuring efficient use of data in the training phase. The final stage, **Storage**, archives this curated dataset into a structured, queryable ASE (Atomic Simulation Environment) database, embedding rich metadata for provenance and reproducibility. The entire process is configurable via the Hydra framework and accessible through both a powerful Command-Line Interface (CLI) for batch processing and an intuitive Web-based User Interface (Web UI) for interactive use, making it a comprehensive solution for the materials science community.

## 2. System Design Objectives

The design and implementation of MLIP-AutoPipe are guided by a set of core objectives that ensure it is a powerful, reliable, and user-friendly scientific tool. These objectives address the technical, scientific, and usability aspects of the system.

*   **Physical Realism and Constraint Enforcement**: This is the most fundamental objective. The framework must generate atomic structures that are physically meaningful and adhere to the laws of chemistry and physics. This involves several layers of validation. At the most basic level, it must prevent unphysical configurations like atoms overlapping (being too close to each other). For specific material types, it must enforce additional constraints, such as ensuring charge neutrality in ionic crystals by maintaining the correct stoichiometry of cations and anions. Furthermore, when simulating systems under periodic boundary conditions, the simulation cell size must be sufficiently large relative to the MLIP's cutoff radius to avoid atoms interacting with their own periodic images. The system must be designed to rigorously check for and reject any structure that violates these fundamental physical constraints at every stage of the pipeline. The ultimate goal is to produce datasets that are a faithful representation of the material's true potential energy surface.

*   **Maximising Diversity and Robustness**: A key differentiator for MLIP-AutoPipe is its ability to generate diverse datasets that enhance the robustness of the trained MLIP. The design objective is to move beyond just sampling near the equilibrium state. The system must be capable of exploring a wide swath of the thermodynamic phase space. This means the exploration engine must be designed to operate at high temperatures and pressures, pushing the system into high-energy states, inducing phase transitions, and discovering various metastable configurations. The inclusion of hybrid MD/MC methods is central to this objective, as MC moves like atom swaps can overcome the high energy barriers that might trap a pure MD simulation in a single compositional ordering. The system should not just create stable structures; it should create a rich ensemble of configurations representing the material under a variety of conditions, including those close to points of failure or transition.

*   **Automation, Efficiency, and Scalability**: The system must be designed for full automation, enabling a "fire-and-forget" workflow. A user should be able to define a scientific problem through a configuration file, launch the pipeline, and have it run to completion without any further intervention. Efficiency is paramount, given the computationally expensive nature of the simulations. The design must incorporate parallelism at its core, using Python's modern multiprocessing capabilities (`ProcessPoolExecutor`) to run multiple independent simulations concurrently, thereby leveraging multi-core CPUs effectively. Furthermore, it must be mindful of memory usage and I/O overhead. Advanced techniques like the "late-binding" of the MLIP calculator (instantiating it inside worker processes) are critical design choices to avoid the performance pitfalls of serializing large ML models. The architecture should also be scalable, capable of managing hundreds or thousands of individual simulations in a single logical run.

*   **Modularity, Extensibility, and Maintainability**: The software architecture must be highly modular to facilitate long-term maintenance and future development. Each distinct function of the pipeline (Generation, Exploration, Sampling, Storage) is to be implemented as a separate, decoupled component. Communication between these components will be handled through well-defined interfaces (Abstract Base Classes). This design pattern makes the system extensible. For instance, adding a new type of material generator (e.g., a `PolymerGenerator`) would simply involve creating a new class that implements the `IStructureGenerator` interface, without requiring any changes to the core orchestration logic. A factory pattern will be used to dynamically select and instantiate the correct component at runtime based on the user's configuration. This approach ensures the codebase remains clean, organized, and easy for new developers to contribute to.

*   **User Accessibility and Dual Interfaces**: The system must be accessible to a broad range of users with different levels of technical expertise. To achieve this, the design specifies two primary user interfaces. First, a comprehensive Command-Line Interface (CLI) is required for power users, scripting, and integration into larger automated computational workflows. The CLI will be the primary interface for batch processing and reproducible research. Second, a browser-based Web User Interface (Web UI) is required for interactive use, visualization, and educational purposes. The Web UI will provide a guided, graphical way to construct simulation configurations, launch runs, monitor their progress in real-time, and, most importantly, visualize the atomic structures produced. This dual-interface approach ensures that the powerful backend of the system can be leveraged by the widest possible audience.

*   **Configuration Management and Scientific Reproducibility**: Absolute reproducibility is a cornerstone of good science. The design must ensure that any dataset generated by MLIP-AutoPipe can be perfectly reproduced at a later date. This is achieved by using the Hydra framework for configuration management. Every parameter that influences the data generation process, from the definition of the chemical system to the intricacies of the simulation and sampling algorithms, will be stored in structured, human-readable YAML files. A given run is entirely defined by its configuration. By archiving this configuration file alongside the generated dataset, the exact process can be repeated, ensuring full scientific provenance and satisfying the requirements of peer-reviewed publication and data sharing.

## 3. System Architecture

The MLIP-AutoPipe framework is architected as a modular, sequential pipeline driven by a central orchestrator. This design ensures a clean separation of concerns, where each component has a single, well-defined responsibility. The orchestrator, named `PipelineOrchestrator`, is responsible for managing the flow of data and control between the four primary stages of the process. This architecture is designed for clarity, robustness, and extensibility. The user interacts with the system through either the CLI or the Web UI, which translates user input into a single, comprehensive configuration object. This configuration object is the sole input to the `PipelineOrchestrator`, which then executes the pipeline.

```mermaid
graph TD
    subgraph User Interaction Layer
        A[Command-Line Interface (CLI)]
        B[Web User Interface (Web UI)]
        A -- Provides --> C(YAML Configuration File)
        B -- Generates --> C
    end

    subgraph Orchestration Layer
        D(PipelineOrchestrator)
        C -- Is Parsed and Validated by --> D
    end

    subgraph Core Pipeline Stages
        D -- Executes Stage 1 --> E{Generation}
        E -- Seed Structures --> F{Exploration}
        F -- Trajectory Data --> G{Sampling}
        G -- Curated Structures --> H{Storage}
    end

    subgraph Data Persistence Layer
        E --> I[Intermediate Files: initial_structures.xyz]
        F --> J[Intermediate Files: trajectories/*.traj]
        H --> K[Final Output: ASE Database (dataset.db)]
    end

    subgraph Pluggable Components
        E --> E1[AlloyGenerator]
        E --> E2[IonicGenerator]
        E --> E3[InterfaceGenerator]
        E --> E4[...]

        F --> F1[MD/MC Explorer Engine]

        G --> G1[RandomSampler]
        G --> G2[FPSSampler]
    end

    style D fill:#f9f,stroke:#333,stroke-width:2px
    style H fill:#ccf,stroke:#333,stroke-width:2px
```

**Component Breakdown:**

1.  **Orchestrator (`PipelineOrchestrator`)**: This class serves as the brain of the entire application. Its primary role is to execute the data generation workflow according to the user's configuration. It initializes all the necessary service components (generators, explorers, samplers) using a factory pattern, which reads the configuration and selects the appropriate implementations. The orchestrator manages the state of the pipeline, calling each stage in sequence and ensuring that the output of one stage is correctly passed as the input to the next. It is also responsible for handling intermediate file I/O, writing checkpoints like the initial structures and trajectory files to disk. This ensures that the state is preserved and allows for potential recovery or analysis of intermediate steps.

2.  **Generators (e.g., `AlloyGenerator`)**: These components are responsible for the very first step: creating the initial population of atomic structures (seeds). The architecture uses an interface-based design, with a `IStructureGenerator` abstract base class defining a common `generate()` method. This allows for multiple, specialized generator implementations for different types of materials. For example, the `AlloyGenerator` creates random solid solutions, while the `IonicGenerator` ensures charge balance. These generators are not just creating random atoms; they embed physical knowledge, applying initial volumetric strains and atomic "rattles" to create a diverse and physically plausible starting point for the exploration stage. They also perform critical validation to ensure no two atoms are too close.

3.  **Explorer (`MDMCExplorer`)**: This is the most computationally intensive component. It is a sophisticated service that takes the seed structures and runs MD or hybrid MD/MC simulations to explore the potential energy surface. The design emphasizes robustness and physical accuracy. A key feature is its ability to automatically detect the system's geometry (bulk vs. slab) and select the correct thermodynamic ensemble (NPT vs. NVT), a crucial step to avoid unphysical simulation artifacts. It also integrates a mixed potential, combining the MLIP with a classical ZBL potential, to handle high-energy collisions and prevent simulation failures. To achieve high throughput, the explorer is designed to run simulations in parallel using a pool of worker processes. A critical architectural choice is the "late binding" of the MLIP calculator, meaning the large model is loaded inside each worker process, avoiding costly serialization and inter-process communication bottlenecks.

4.  **Samplers (e.g., `FPSSampler`)**: After the exploration stage has generated potentially millions of structures within the trajectory files, the sampler's job is to select a small, high-quality subset. A simple `RandomSampler` provides a baseline capability. The more advanced `FPSSampler` (Farthest Point Sampling) implements a sophisticated algorithm to select a subset of structures that are maximally diverse. It does this by calculating a "fingerprint" (using SOAP descriptors) for each structure and then iteratively selecting the structure that is farthest (i.e., most different) from the ones already selected. This is a crucial step for creating efficient training datasets, as it reduces redundancy and ensures the final dataset covers the explored structural space as broadly as possible.

5.  **Storage (`AseDBWrapper`)**: The final component in the pipeline is responsible for persistent storage. It takes the curated list of structures from the sampler and saves them into an ASE (Atomic Simulation Environment) database. This is a standard, portable format in the computational materials science community, based on an SQLite file. The `AseDBWrapper` abstracts away the details of the database interaction, providing simple methods to write the `ase.Atoms` objects along with all their associated metadata, such as calculated energies, forces, stresses, and the configuration parameters used to generate them. This ensures the final dataset is self-contained and ready for downstream use in training MLIP models.

## 4. Design Architecture

The project's design architecture is based on modern, best-practice principles for building robust and maintainable Python applications. It employs a layered architecture that strictly separates concerns, a schema-first approach using Pydantic for data modeling, and interface-based design for service components. This ensures the system is type-safe, validated, and easily extensible.

**File Structure:**

The directory structure is organized to reflect the layered architecture, making the codebase intuitive to navigate.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py                  # Presentation Layer: Main CLI entry point (using Click)
├── web_ui.py               # Presentation Layer: Web UI application (using FastAPI)
│
├── core/                   # Application Core / Orchestration Layer
│   ├── __init__.py
│   ├── pipeline_orchestrator.py # Contains the main PipelineOrchestrator class
│   └── factories.py             # Responsible for creating service components based on config
│
├── domain/                   # Domain Layer: Core data structures and interfaces
│   ├── __init__.py
│   ├── models.py                # Pydantic models for configuration and data objects (e.g., DFTResult)
│   └── interfaces.py            # Abstract Base Classes defining service contracts (e.g., IStructureGenerator)
│
├── infrastructure/           # Infrastructure Layer: Wrappers for external systems
│   ├── __init__.py
│   ├── ase_db_wrapper.py        # Wrapper for all interactions with the ASE SQLite database
│   └── process_runner.py        # Generic wrapper for running external command-line processes
│
└── services/                 # Service Layer: Implementations of business logic
    ├── __init__.py
    ├── generation/              # Generation-related services
    │   ├── __init__.py
    │   └── alloy_generator.py   # Concrete implementation for alloys
    │   └── ...
    ├── exploration/             # Exploration-related services
    │   ├── __init__.py
    │   └── md_mc_explorer.py    # Concrete implementation of the exploration engine
    └── sampling/                # Sampling-related services
        ├── __init__.py
        ├── random_sampler.py    # Concrete implementation for random sampling
        └── fps_sampler.py       # Concrete implementation for Farthest Point Sampling
```

**Data Models (Pydantic):**

The "schema-first" design philosophy is central to the architecture. Before any logic is implemented, the data structures that flow through the system are rigorously defined using Pydantic. This provides several key advantages: automatic input validation, clear and self-documenting code, and easy serialization/deserialization from formats like YAML and JSON.

*   `SystemConfig`: This model encapsulates everything needed to define the chemical system itself. It includes fields for the list of elements, their relative compositions, the crystal lattice type (e.g., 'fcc', 'bcc'), and the number of initial structures to generate. Validators on this model will enforce invariants, such as ensuring the compositions sum to 1.0.
*   `ExplorationConfig`: This large model contains all the tunable parameters for the computationally intensive exploration stage. It includes fields for the path to the MLIP model, simulation parameters (temperature, pressure, number of steps), the choice of thermodynamic ensemble, and boolean flags to enable advanced features like hybrid MD/MC or ZBL potential mixing.
*   `SamplingConfig`: This model defines how the final dataset is curated from the raw trajectory data. It includes a field for the sampling method ('random' or 'fps') and the target number of samples. It also contains an optional sub-model for method-specific parameters, like the SOAP descriptor settings required for FPS.
*   `FullConfig`: This is the top-level Pydantic model that aggregates all the other configuration models into a single, cohesive object. The entire state of a data generation run is defined by this single model, which is what gets serialized to the YAML configuration file.
*   `DFTResult`: Beyond just configuration, Pydantic will also be used for data transfer objects. This model defines the structure for storing the results of a quantum mechanical calculation (energy, forces, stresses), which can then be consistently stored as a single entity in the ASE database's `data` field.

This schema-driven approach ensures that every component has a clear and validated contract for the data it consumes and produces, significantly reducing the likelihood of runtime errors and making the system as a whole more robust and reliable.

## 5. Implementation Plan

The project's implementation is strategically divided into two distinct, sequential cycles. This incremental approach allows for the delivery of a functional core product in the first cycle, followed by the addition of more advanced features and a graphical user interface in the second. This methodology manages complexity, reduces risk, and allows for feedback and testing on the core functionality before more complex features are added.

**Cycle 1: Core CLI Pipeline**

The exclusive focus of the first cycle is to build and deliver a fully functional, end-to-end data generation pipeline that is operated entirely from the command line. The goal is to create a robust and reliable tool for power users and automated workflows. This Minimal Viable Product (MVP) will validate the entire architecture and provide tangible value early in the development process. The user will define their desired dataset via a YAML configuration file and use a simple CLI command to execute the entire workflow.

*   **Features**:
    *   **CLI Development**: Implement the main command-line entry point using the `click` library. This will involve creating the main `run` command, which will accept arguments for the path to the configuration file and the desired output path for the database.
    *   **Configuration Modeling**: Develop the complete set of Pydantic models (`SystemConfig`, `ExplorationConfig`, `SamplingConfig`, `FullConfig`) that define the schema for the YAML configuration files. This includes adding appropriate validators to ensure data integrity.
    *   **Orchestration Logic**: Implement the `PipelineOrchestrator` class, which will serve as the central controller for the workflow. This class will be responsible for parsing the configuration, instantiating the necessary services using a factory, and executing the four stages (Generation, Exploration, Sampling, Storage) in the correct sequence.
    *   **Structure Generators**: Implement at least two fundamental structure generators: the `AlloyGenerator` for creating metallic alloys and solid solutions, and the `IonicGenerator` for creating charge-neutral ionic compounds. These will include essential physical validation checks.
    *   **Exploration Engine (Core)**: Develop the core functionality of the `MDMCExplorer`. This will involve setting up and running standard MD simulations using an MLIP calculator. A key part of this task is the implementation of the parallel execution logic using `ProcessPoolExecutor` to run multiple simulations concurrently.
    *   **Basic Sampler**: Implement a `RandomSampler` as the baseline method for selecting structures from the generated trajectories.
    *   **Database Storage**: Develop the `AseDBWrapper` class, which will encapsulate all interactions with the ASE database, providing a clean API for writing the final curated structures.
    *   **Testing Foundation**: Establish a comprehensive suite of unit and integration tests for all the components developed in this cycle, ensuring the reliability of the core pipeline.

**Cycle 2: Advanced Features and Web UI**

The second cycle builds upon the stable foundation laid in Cycle 1. Its focus is on enhancing the scientific sophistication of the exploration engine and dramatically improving user accessibility by adding a graphical Web UI.

*   **Features**:
    *   **Advanced Exploration Engine**: Significantly enhance the `MDMCExplorer` by adding more complex simulation capabilities.
        *   Implement **Hybrid MD/MC**: Add the ability to perform Monte Carlo moves, such as atom swaps and vacancy hops, during the simulation to more effectively explore configurational and compositional space.
        *   Implement **Automatic Ensemble Switching**: Develop the algorithm to detect vacuum slabs in a simulation cell and automatically switch from the NPT to the NVT ensemble to prevent simulation artifacts.
        *   Implement **ZBL Potential Integration**: Add the logic to mix the MLIP with the classical ZBL potential to provide a robust repulsive wall at short distances, stabilizing high-temperature simulations.
    *   **Intelligent Sampling**: Implement the `FPSSampler` (Farthest Point Sampling). This involves integrating a library for calculating SOAP descriptors and implementing the FPS algorithm to select a maximally diverse set of structures.
    *   **Web UI Development**: Develop a full-featured, browser-based graphical user interface using a modern Python web framework like FastAPI. The UI will provide an intuitive, form-based way for users to:
        *   Create and edit pipeline configurations without needing to write YAML manually.
        *   Launch new data generation runs and monitor their progress in real-time.
        *   Browse, search, and visualize the atomic structures stored in the output databases using an integrated 3D viewer.
    *   **GPU Resource Management**: Refine the parallel processing logic to be more aware of GPU resources, ensuring that parallel runs do not overload the available GPU memory when using GPU-accelerated MLIP models.
    *   **Specialized Generators**: Add new, more specialized structure generators, such as an `InterfaceGenerator` for creating heterostructures and an `AdsorptionGenerator` for placing molecules on surfaces.

## 6. Test Strategy

A multi-layered testing strategy is essential to ensure the correctness, reliability, and scientific validity of the MLIP-AutoPipe framework. The strategy encompasses unit tests for individual components, integration tests for the complete pipeline, and end-to-end tests for the user interfaces.

**Cycle 1 Test Strategy:**

The testing in Cycle 1 is focused on the backend components and the CLI-driven workflow.

*   **Unit Testing**: Each service and infrastructure component will be tested in complete isolation to verify its logic.
    *   **Generators**: The `AlloyGenerator` and `IonicGenerator` will be tested to confirm they produce `ase.Atoms` objects that strictly adhere to the requested physical constraints. Tests will involve mocking the random number generator for reproducibility and asserting that the output structures have the correct composition, lattice, and do not contain overlapping atoms. A key test will be to ensure that requesting an impossible structure (e.g., one that would force atoms to overlap) raises a specific `PhysicsViolationError`.
    *   **Explorer**: Unit testing the `MDMCExplorer` is challenging due to its complexity. The focus will be on the orchestration logic rather than the physics. We will use a very fast, simple classical potential (like EMT) instead of a real MLIP. Tests will verify that the `ProcessPoolExecutor` is correctly initialized, that simulations are launched for all seed structures, and that the "late-binding" of the calculator works as expected without causing pickling errors.
    *   **Samplers**: The `RandomSampler` will be tested with a dummy trajectory file. The test will simply assert that the sampler returns a list of `ase.Atoms` objects of the correct, specified length.
    *   **Storage**: The `AseDBWrapper` will be tested against a temporary SQLite database file. The tests will perform a full round-trip check: write a known list of `Atoms` objects to the database, then read them back and assert that the retrieved data (atomic numbers, positions, energy, forces) is identical to the original data.
*   **Integration Testing**: The primary goal of integration testing in this cycle is to verify that the `PipelineOrchestrator` correctly wires together all the individual services and manages the data flow from start to finish.
    *   A comprehensive end-to-end test will be created. It will use the `click.testing.CliRunner` to invoke the main CLI command with a minimal configuration file (e.g., a simple 2-atom Si system, 10 MD steps, random sampling). The test will assert that the command completes with an exit code of 0. It will then inspect the filesystem to verify that intermediate files were created and, most importantly, it will connect to the final output database and assert that it contains the expected number of structures. This single test provides high confidence in the correctness of the entire backend pipeline.

**Cycle 2 Test Strategy:**

Testing in Cycle 2 expands to cover the new advanced features and the Web UI.

*   **Unit Testing**:
    *   **Advanced Explorer Features**: The new, complex logic in the `MDMCExplorer` will require dedicated unit tests. The automatic ensemble switching feature will be tested by creating known bulk and slab `Atoms` objects and asserting that the correct ASE dynamics object (`NPT` vs. `NVT`) is chosen. The logic for mixing the ZBL and MLIP potentials will be tested by checking that the forces on a simple two-atom system are the correct sum of the forces from the individual potentials.
    *   **FPSSampler**: The `FPSSampler` will be tested with a carefully crafted dummy trajectory containing both redundant and unique structures. The test will mock the SOAP descriptor calculation to provide deterministic fingerprints and will assert that the sampler correctly selects the set of unique structures.
*   **Integration Testing**: The integration tests from Cycle 1 will be extended. New test cases will be added that enable the advanced exploration features (hybrid MD/MC, ZBL) and the `FPSSampler` in the configuration file, verifying that these more complex pipeline configurations also run to completion without errors.
*   **End-to-End (E2E) Testing for Web UI**: The introduction of the Web UI necessitates a new layer of testing.
    *   We will use a browser automation framework like `pytest-playwright`. The E2E tests will script a real user's interaction with the system through a browser.
    *   A "golden path" test will launch the web application, navigate to its URL, programmatically fill out the web form to configure a run (including enabling new Cycle 2 features), and click the "submit" button.
    *   The test will then assert that the UI transitions to a "running" state. It will poll a backend API status endpoint until the run is complete.
    *   Finally, it will navigate to the results page and assert that the output data is correctly displayed in the results table and that the 3D structure viewer can be loaded. This ensures that the entire stack, from the frontend JavaScript to the backend web server and the core scientific pipeline, is correctly integrated.
