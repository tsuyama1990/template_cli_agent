# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe project is a state-of-the-art computational framework designed to automate the generation of high-quality training datasets for Machine Learning Interatomic Potentials (MLIPs), such as MACE and SevenNet. The core philosophy of the system, as derived from the requirements, is to "remove the human expert from the loop" by creating a fully automated, physically robust, and scientifically rigorous pipeline. The primary challenge in developing accurate MLIPs is the creation of a dataset that is not only diverse but also encompasses the vast and complex potential energy surface (PES) of a material. This includes stable low-energy configurations, high-energy states, transition states, and even defect structures, all of which are critical for the model's predictive power and robustness. The MLIP-AutoPipe system directly addresses this challenge by moving beyond simple, random structure generation. It employs a sophisticated, multi-stage pipeline that simulates realistic thermodynamic conditions to explore the accessible phase space of a given material system.

The pipeline begins with a 'Generation' stage, where initial seed structures are created based on user-defined physical parameters (e.g., chemical composition, crystal structure). The system is designed to handle a wide variety of material classes, including alloys, ionic crystals, covalent materials, interfaces, and surfaces with adsorbates. This stage incorporates fundamental physical constraints, such as ensuring minimum atomic separation distances and appropriate supercell sizes to avoid periodic boundary condition artifacts. The second stage, 'Exploration', is the scientific core of the framework. It uses molecular dynamics (MD) and hybrid Monte Carlo (MC) simulations to evolve the initial structures. This is not a simple simulation; the engine is engineered with advanced features like automatic ensemble switching (NPT for bulk, NVT for surfaces) to apply the correct thermodynamic conditions, and it integrates classical potentials like ZBL to handle high-energy atomic collisions gracefully. This prevents common simulation failures and ensures the generation of physically meaningful, albeit sometimes "unlikely," configurations that are crucial for training a robust MLIP.

The third stage, 'Sampling', intelligently selects a subset of structures from the vast trajectory data generated during exploration. Instead of random selection, it employs techniques like Farthest Point Sampling (FPS) using SOAP descriptors to maximize structural diversity and ensure the final dataset is information-rich. Finally, the 'Storage' stage archives these carefully selected structures and their associated metadata into a structured, queryable ASE database. The entire process is orchestrated by a robust runner that manages parallelism, checkpointing, and error handling, ensuring reliability even during long, computationally intensive runs. The system is configured via the Hydra framework, offering exceptional flexibility, and is designed with two primary user interfaces: a powerful Command-Line Interface (CLI) for automated, batch-driven workflows and a user-friendly Web UI for interactive exploration, configuration, and visualization. This dual-interface approach caters to both power users scripting large-scale campaigns and researchers who prefer a more visual, hands-on approach. The project's engineering quality is evident in its thoughtful design, which includes features like late-binding of MLIP calculators to optimize parallel performance and physics-based validation at every step, making it a truly expert-level tool for materials science research.

## 2. System Design Objectives

The principal objective of the MLIP-AutoPipe system is to democratise and accelerate the development of high-fidelity MLIPs by automating the most laborious and expertise-intensive part of the workflow: dataset generation. The design is guided by a set of clear goals, constraints, and success criteria.

**Goals:**
1.  **Automation:** To create a "turnkey" solution where a user can define a material system in a high-level configuration file and receive a high-quality, ready-to-use training dataset with minimal manual intervention. The system must handle all intermediate steps, from structure creation to simulation and data curation.
2.  **Physical Realism:** The generated structures must be physically meaningful. This means they must adhere to fundamental chemical and physical constraints (e.g., charge neutrality in ionic systems, correct bond distances) and be representative of configurations accessible under realistic temperature and pressure conditions.
3.  **Diversity and Robustness:** The dataset must be diverse, spanning a wide range of the potential energy surface. The system is explicitly designed to find and include "hard" cases—high-energy configurations, transition states, and defect structures—that are essential for training a robust MLIP capable of handling unseen configurations.
4.  **Extensibility and Modularity:** The architecture must be modular, allowing for the easy addition of new structure generators, exploration algorithms, or sampling methods. The use of a `BaseGenerator` class and a factory pattern is a direct implementation of this goal.
5.  **User-Friendliness and Accessibility:** The system must be accessible to a broad range of users. The CLI provides a powerful interface for automation and scripting, while the Web UI will offer an intuitive, visual way to interact with the pipeline, lowering the barrier to entry for non-experts.

**Constraints:**
1.  **Computational Resources:** The exploration phase is computationally expensive. The system must be designed to run efficiently on modern multi-core CPUs and, where applicable, GPUs. It must manage parallel processes intelligently, avoiding common pitfalls like memory bottlenecks or GPU resource contention.
2.  **Dependency Management:** The project relies on a complex stack of scientific libraries (ASE, PyMatGen, MACE, etc.). The system must have a clear and reproducible dependency management strategy to ensure that it can be deployed reliably across different environments.
3.  **Data Integrity:** The pipeline generates large amounts of data. It must be robust against failures, using checkpointing and progressive saving to prevent data loss during long-running simulations. The final database must be structured and contain all necessary metadata for reproducibility.

**Success Criteria:**
1.  **Successful Pipeline Execution:** A user can successfully run the entire pipeline from the command line for a standard alloy system (e.g., FePt), generating an ASE database containing at least 100 unique structures.
2.  **MLIP Training Viability:** A dataset generated by the system for a well-studied material can be used to train an MLIP (e.g., MACE) that achieves a level of accuracy comparable to literature benchmarks for key properties like forces and energies.
3.  **Web UI Functionality:** A user can use the web interface to define a simple workflow, launch the calculation, monitor its progress, and visualize the final generated structures.
4.  **Robustness Demonstration:** The system can successfully complete a high-temperature (e.g., 2000K) simulation run without crashing due to "Coulomb explosion" or other physical violations, demonstrating the effectiveness of its safety mechanisms.

## 3. System Architecture

The MLIP-AutoPipe system is designed as a modular, pipeline-based architecture. The core logic is decoupled from the user interface, allowing the same backend engine to be driven by either the CLI or the Web UI. The architecture is composed of four main stages, orchestrated by a central `PipelineRunner`.

**Components:**
1.  **Configuration (Hydra):** This is the entry point for the user. Hydra provides a flexible and powerful way to define all aspects of the workflow, from the material system to the simulation parameters and sampling strategy. This component is responsible for parsing the user's input and creating a structured configuration object.
2.  **Pipeline Runner:** This is the central orchestrator. It takes the configuration object and executes the four main stages of the pipeline in sequence. It manages the overall workflow, handles checkpointing and restart logic, and controls the parallel execution of the exploration stage.
3.  **Stage 1: Generation:** This component is a factory that produces initial seed structures. It contains a collection of `Generator` modules, each specialized for a different type of material (e.g., `AlloyGenerator`, `IonicGenerator`). The appropriate generator is selected based on the configuration. This stage applies initial physical validation and constraints.
4.  **Stage 2: Exploration:** This is the computational core. The `Explorer` component takes the seed structures and runs MD or hybrid MD/MC simulations. It is responsible for managing the physics of the simulation, including the choice of thermodynamic ensemble (NPT/NVT), the integration of MLIP and classical potentials (e.g., ZBL), and handling physical violations. It leverages Python's `ProcessPoolExecutor` for parallelism.
5.  **Stage 3: Sampling:** This component processes the raw trajectory data from the exploration stage. It contains different `Sampler` implementations (e.g., `RandomSampler`, `FPSSampler`) that select a representative subset of structures based on the chosen strategy.
6.  **Stage 4: Storage:** The final component, `Storage`, takes the sampled structures and saves them into a persistent database. It uses the ASE DB module to create a structured SQLite database, storing not only the atomic coordinates but also relevant metadata like potential energy, forces, and simulation parameters.

**Data Flow:**
The data flow is sequential and unidirectional. The `Configuration` object is passed to the `PipelineRunner`. The `Generation` stage outputs a set of `ASE.Atoms` objects, which are saved to disk and passed to the `Exploration` stage. The `Exploration` stage produces large trajectory files (e.g., `.xyz`). The `Sampling` stage reads these files and outputs a smaller, curated list of `ASE.Atoms` objects. Finally, the `Storage` stage takes this list and commits it to the final ASE database. This clear, staged flow with state isolation (saving intermediate results to disk) enhances robustness and modularity.

```mermaid
graph TD
    subgraph User Interface
        A[CLI / Web UI]
    end

    subgraph Core Pipeline
        B[Configuration Parser <br> (Hydra)]
        C[Pipeline Runner]
        D[Stage 1: Generation Engine <br> (Factory Pattern)]
        E[Stage 2: Exploration Engine <br> (MD/MC)]
        F[Stage 3: Sampling Engine <br> (FPS/Random)]
        G[Stage 4: Storage Engine <br> (ASE DB)]
    end

    subgraph Data
        H[YAML Config Files]
        I[Intermediate Files <br> (initial_structures.xyz)]
        J[Trajectory Files <br> (trajectory.xyz)]
        K[Final Database <br> (structures.db)]
    end

    A --> H
    H --> B
    B --> C
    C --> D
    D -- ASE.Atoms objects --> I
    I --> E
    E -- Trajectory Data --> J
    J --> F
    F -- Curated Atoms --> G
    G -- Structured Data --> K

    style A fill:#cde4ff
    style K fill:#d5f0d5
```

## 4. Design Architecture

The design architecture of MLIP-AutoPipe emphasizes modularity, abstraction, and adherence to good software engineering principles. The codebase will be organized into a clear, hierarchical package structure within the `src/mlip_autopipec/` directory.

**File Structure:**
```
src/mlip_autopipec/
├── __init__.py
├── cli.py                  # Main entry point for the Command-Line Interface (Typer/Click)
├── web_ui.py               # Main entry point for the Web UI (e.g., Streamlit/FastAPI)
│
├── config/                 # Hydra configuration files (.yaml)
│   ├── config.yaml
│   └── system/
│       ├── alloy.yaml
│       └── ionic.yaml
│
├── pipeline/
│   ├── __init__.py
│   └── runner.py             # Contains the PipelineRunner orchestrator class
│
├── generators/
│   ├── __init__.py
│   ├── base.py               # Abstract Base Class for all generators
│   ├── alloy.py              # Generator for alloy structures
│   ├── ionic.py              # Generator for ionic crystals
│   └── factory.py            # Factory to select the appropriate generator
│
├── explorers/
│   ├── __init__.py
│   └── md_engine.py          # Core MD/MC simulation engine
│
├── sampling/
│   ├── __init__.py
│   ├── base.py               # Abstract Base Class for samplers
│   ├── random_sampler.py
│   └── fps_sampler.py
│
├── storage/
│   ├── __init__.py
│   └── ase_db_writer.py      # Module for writing to ASE DB
│
└── common/
    ├── __init__.py
    ├── atoms_validator.py    # Physics-based validation logic
    └── exceptions.py         # Custom exception classes
```

**Class/Function Definitions Overview:**

*   **`cli.py`**: Will contain a main function decorated with `@app.command()` (from Typer) that initializes the Hydra configuration and invokes the `PipelineRunner`.
*   **`pipeline/runner.py`**:
    *   `PipelineRunner`: The main class with a `run()` method. It will instantiate and call the necessary classes from the `generators`, `explorers`, `sampling`, and `storage` modules in sequence.
*   **`generators/`**:
    *   `BaseStructureGenerator` (ABC): Defines the interface for all generators, including a `generate()` abstract method and common validation logic.
    *   `AlloyGenerator(BaseStructureGenerator)`: Implements the logic for creating alloy structures based on composition and crystal structure.
    *   `GeneratorFactory`: A class with a static method `get_generator(config)` that returns an instance of the correct generator class.
*   **`explorers/md_engine.py`**:
    *   `MDEngine`: A class that encapsulates the logic for running a single MD/MC simulation for one `Atoms` object. It will handle calculator setup (MACE + ZBL), thermostat/barostat integration, and MC moves.
*   **`sampling/`**:
    *   `BaseSampler` (ABC): Defines the `sample(trajectory_file)` interface.
    *   `RandomSampler` and `FPSSampler`: Concrete implementations of the sampling logic.
*   **`storage/ase_db_writer.py`**:
    *   `save_to_db(atoms_list, db_path)`: A function that takes a list of `Atoms` objects and writes them to the specified SQLite database path.

**Data Models (Pydantic):**
While Hydra manages configuration, Pydantic models will be used internally for data validation and to define clear data structures, especially for complex configuration objects that are passed between components. For instance, the configuration for a material system could be encapsulated in a Pydantic model to ensure all required fields (e.g., elements, concentrations) are present and have valid types before being passed to a generator. This adds a layer of robustness on top of Hydra's flexibility.

## 5. Implementation Plan

The project will be developed over two distinct cycles, ensuring a focused, iterative delivery of functionality.

**Cycle 1: Core CLI Pipeline and Foundational Components**
This cycle focuses on building the complete, end-to-end automated pipeline accessible from the command line. The goal is to deliver a fully functional tool for power users that can generate a high-quality dataset for a simple alloy system. This cycle lays the foundational codebase for all future development.
*   **Features:**
    1.  **Project Scaffolding:** Set up the repository structure, `pyproject.toml`, and dependency management with `uv`.
    2.  **Configuration:** Implement the Hydra configuration structure for defining a simple alloy system (`elements`, `composition`, `crystal_structure`) and basic exploration parameters (`temperature`, `steps`).
    3.  **CLI Entry Point:** Create the main CLI entry point using Typer (`cli.py`).
    4.  **Generation Module:** Implement the `BaseStructureGenerator` and a concrete `AlloyGenerator`. Implement the `GeneratorFactory`.
    5.  **Exploration Module:** Implement the core `MDEngine` to run a basic NVT molecular dynamics simulation using an existing MLIP like MACE. Focus on the core simulation loop and calculator integration.
    6.  **Sampling Module:** Implement a simple `RandomSampler` as the initial sampling strategy.
    7.  **Storage Module:** Implement the `ase_db_writer` to save the final structures to an ASE database.
    8.  **Pipeline Orchestration:** Implement the `PipelineRunner` to connect all the modules and execute the full workflow.
*   **Outcome:** A functional CLI tool (`mlip-autopipec run --config-name=alloy`) that takes a YAML file and produces an ASE database of atomic structures.

**Cycle 2: Advanced Features and Web UI**
This cycle builds upon the stable foundation of Cycle 1. It introduces more sophisticated scientific features and develops the user-friendly Web UI. The goal is to make the tool more powerful for expert users and more accessible for non-experts.
*   **Features:**
    1.  **Advanced Exploration:** Enhance the `MDEngine` with hybrid MD/MC capabilities (e.g., atom swaps), automatic ensemble switching (NPT/NVT), and the integration of the ZBL potential for improved physical realism at high energies.
    2.  **Advanced Sampling:** Implement the `FPSSampler` using SOAP descriptors to enable intelligent, diversity-driven sampling.
    3.  **Expanded Generators:** Add more generator classes to handle other material types, such as `IonicGenerator` and a knowledge-based generator for specific chemical formulas.
    4.  **Web UI Scaffolding:** Set up a basic web application using a framework like Streamlit or FastAPI.
    5.  **Interactive Configuration:** Create UI components (sliders, text inputs, dropdowns) that allow a user to build a configuration for a workflow interactively.
    6.  **Job Execution and Monitoring:** Implement the logic for the Web UI to trigger the `PipelineRunner` as a background process and provide real-time feedback on its progress (e.g., current step, simulation progress).
    7.  **Results Visualization:** Add a component to the Web UI that can read the final ASE database and display the generated structures, allowing for 3D visualization and inspection.
*   **Outcome:** A comprehensive suite of tools, including an enhanced CLI with advanced physics capabilities and a fully functional Web UI for interactive workflow design, execution, and analysis.

## 6. Test Strategy

The test strategy will be multi-layered, encompassing unit, integration, and user acceptance testing to ensure the reliability, correctness, and usability of the framework across both cycles.

**Cycle 1: CLI and Core Logic Testing**
*   **Unit Testing:** Each component will be tested in isolation.
    *   **Generators:** Tests will verify that generators produce the correct number of structures, adhere to physical constraints (e.g., minimum atomic distances), and have the correct chemical composition. Mocking will be used to isolate the generator from filesystem dependencies.
    *   **Samplers:** The `RandomSampler` will be tested to ensure it selects the correct number of frames from a sample trajectory file.
    *   **Validators:** The `atoms_validator` logic will have dedicated tests to check its ability to detect invalid structures (e.g., overlapping atoms).
    *   **Runner:** The `PipelineRunner`'s logic will be tested by mocking the individual pipeline stages (Generate, Explore, etc.) to verify that they are called in the correct order and with the correct parameters.
*   **Integration Testing:** Tests will cover the interaction between components.
    *   A primary integration test will run the entire pipeline via the CLI (`CliRunner` from Typer/Click's testing tools) on a minimal, well-defined test case (e.g., generating 2 structures of a 4-atom Si cell).
    *   This test will use a fast, simple calculator (like ASE's EMT) instead of a full MLIP to ensure the test runs quickly.
    *   The test will assert that the final database file is created and contains the expected number of structures. It will check the integrity of the data flow from configuration parsing to final database writing.

**Cycle 2: Advanced Features and UI Testing**
*   **Unit Testing:**
    *   **Advanced Exploration:** The logic for automatic ensemble switching (`detect_vacuum`) will be tested with sample "bulk" and "slab" structures. The MC move logic will be tested to ensure it correctly modifies `Atoms` objects (e.g., swaps atoms of different species).
    *   **Advanced Sampling:** The `FPSSampler` will be tested with a known set of structures to ensure its output is deterministic and selects the most diverse structures as expected.
    *   **Web UI Backend:** If using FastAPI, the API endpoints will be tested using an HTTP client library to verify correct request handling and response formatting.
*   **Integration Testing:**
    *   The CLI-based integration test from Cycle 1 will be extended to cover the new, advanced exploration and sampling features.
    *   A new set of integration tests will be developed for the Web UI. These might use tools like Selenium or Playwright to simulate user interactions (filling forms, clicking buttons) and verify that these actions correctly trigger the backend pipeline and that the UI updates as expected.
*   **User Acceptance Testing (UAT):** For the Web UI, UAT will be critical. Test scenarios will be defined (as per `UAT.md`) where a user performs a task, such as "Configure and launch a simulation for amorphous Silica," and the final output and user experience are evaluated against the requirements. This ensures the UI is not only functional but also intuitive and effective.
