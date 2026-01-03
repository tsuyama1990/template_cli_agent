# System Architecture: MLIP-AutoPipe

## 1. Summary

The "MLIP-AutoPipe" project is a specialised software framework designed to automate the generation of high-quality, physically realistic datasets for training modern Machine Learning Interatomic Potentials (MLIPs), such as MACE and SevenNet. The core philosophy of this project is to "remove the human expert from the loop," creating a fully automated pipeline that can be operated by materials scientists without requiring deep expertise in simulation-craft or manual data curation. Traditional methods for generating training data often involve a laborious and error-prone process of manual structure creation, simulation setup, and data filtering. This manual intervention introduces bias, limits the diversity of the dataset, and creates a significant bottleneck in the development of new, more accurate MLIPs. MLIP-AutoPipe directly addresses these challenges by providing a robust, reproducible, and highly automated workflow.

The system is engineered to explore the thermodynamic phase space of a given material system efficiently. Instead of generating simple, random atomic configurations, it employs sophisticated simulation techniques, including Molecular Dynamics (MD) and Monte Carlo (MC) methods. This allows the system to discover and collect energetically challenging configurations, such as high-energy states, transition states, and defect structures, which are critical for training a robust potential that can accurately predict material behaviour across a wide range of conditions. The framework is designed to handle a diverse array of physical systems, from simple alloys and ionic crystals to complex interfaces and adsorption scenarios. This versatility is achieved through a modular, factory-based architecture where specific "Generator" components can be selected based on user-defined configuration. The ultimate output is a curated, well-structured database of atomic configurations, complete with relevant metadata, ready for direct use in MLIP training codes. By automating this critical step, MLIP-AutoPipe aims to significantly accelerate the pace of materials discovery and design. The project's success will be measured by its ability to produce datasets that lead to MLIPs with demonstrably superior accuracy and generalisability compared to those trained on manually curated data.

## 2. System Design Objectives

The primary objective of MLIP-AutoPipe is to deliver a fully automated, reliable, and efficient pipeline for MLIP training data generation. The design is guided by several key principles and goals.

**Goals:**
- **Automation:** The entire workflow, from initial structure generation to final database storage, must be executable via a single command-line instruction. This minimises the need for user intervention and ensures reproducibility.
- **Physical Realism:** All generated structures must be physically plausible. The system will incorporate rigorous validation checks, such as minimum interatomic distances and charge neutrality for ionic systems, to discard unrealistic configurations automatically.
- **Data Diversity:** The pipeline must generate a rich and diverse set of atomic configurations. This will be achieved by combining deterministic structure generation methods with stochastic exploration techniques (MD/MC) that sample a wide area of the potential energy surface.
- **Modularity and Extensibility:** The architecture must be modular, allowing for the easy addition of new structure generators, exploration algorithms, or sampling methods without requiring significant changes to the core framework. This ensures the tool can adapt to new scientific challenges and techniques.
- **User-Friendliness:** While the core is a powerful CLI tool, the project also aims to provide a simple web-based user interface (Web UI) for interactive configuration, monitoring, and visualisation of results. This lowers the barrier to entry for researchers who may not be comfortable with command-line environments.

**Constraints:**
- **Dependency Management:** The project will rely on established open-source scientific libraries, primarily the Atomic Simulation Environment (ASE), Pydantic for data validation, and Hydra for configuration management. The choice of dependencies will be carefully managed to avoid conflicts and ensure long-term maintainability.
- **Computational Resources:** The exploration phase (MD/MC simulations) is computationally intensive. The system must be designed to run efficiently on multi-core CPUs and, where applicable, leverage GPUs. It should also handle resource management gracefully, especially in parallel execution environments.
- **Platform Independence:** The core application will be delivered as a containerised CLI tool, ensuring it runs consistently across different operating systems (Linux, macOS, Windows) without complex local environment setup.

**Success Criteria:**
- The final system can successfully generate a database of at least 1,000 unique and physically valid structures for a given alloy system within a reasonable timeframe (e.g., a few hours on a standard research workstation).
- The datasets generated by the pipeline, when used to train a standard MLIP model (e.g., MACE), result in a model with prediction errors (energy, forces) that are at least 10% lower than a model trained on a dataset of equivalent size generated by purely random methods.
- The project is accompanied by a comprehensive test suite with at least 80% code coverage, ensuring the reliability and correctness of the implementation.

## 3. System Architecture

The MLIP-AutoPipe follows a sequential, four-stage pipeline architecture. Each stage is a distinct, self-contained module that performs a specific task and passes its output to the next stage. This design ensures a clear separation of concerns and facilitates independent development and testing of each component.

**Pipeline Stages:**
1.  **Generation:** This initial stage creates a set of "seed" structures based on user configuration. It uses a factory pattern to select the appropriate generator (e.g., `AlloyGenerator`, `IonicGenerator`) for the specified physical system. These initial structures are validated for physical realism and saved to an intermediate store.
2.  **Exploration:** The seed structures are fed into the exploration engine. This module runs MD or hybrid MD/MC simulations to perturb the initial structures, exploring the local energy landscape and generating a large trajectory of configurations. This stage is designed for parallel execution to accelerate the process.
3.  **Sampling:** The raw trajectory data from the exploration stage can be vast and redundant. The sampling module intelligently selects a smaller, more diverse subset of structures that are most valuable for training. It will support multiple strategies, from simple random sampling to more advanced techniques like Farthest Point Sampling (FPS).
4.  **Storage:** The final, curated set of structures is stored in a structured format. An ASE database (backed by SQLite) will be used to save the atomic configurations along with essential metadata, such as potential energy, forces, and the parameters of the simulation from which they were derived.

**Data Flow:**
The data flows linearly through these four stages. The user provides a single configuration file (managed by Hydra). The `PipelineRunner` orchestrates the execution of each stage, ensuring that the output of one stage is correctly formatted and passed as the input to the next.

**Mermaid Diagram:**
```mermaid
graph TD
    A[User Configuration (.yaml)] --> B{PipelineRunner};
    B --> C[Stage 1: Generation];
    C --> D[Stage 2: Exploration];
    D --> E[Stage 3: Sampling];
    E --> F[Stage 4: Storage];
    F --> G[ASE Database (.db)];

    subgraph "Generators"
        C1[AlloyGenerator]
        C2[IonicGenerator]
        C3[...]
    end

    subgraph "Exploration Engines"
        D1[MD Engine]
        D2[Hybrid MD/MC Engine]
    end

    subgraph "Sampling Methods"
        E1[Random Sampling]
        E2[FPS Sampling]
    end

    C -- uses --> C1;
    D -- uses --> D1;
    E -- uses --> E2;
```

This modular architecture ensures that each part of the process can be refined or replaced independently. For instance, a new exploration technique could be added by simply creating a new engine class that conforms to the required interface, without altering the surrounding pipeline.

## 4. Design Architecture

The project will be structured within the `src/mlip_autopipec` directory, following standard Python packaging conventions. The design emphasizes a clean separation between configuration, business logic, and data persistence.

**File Structure:**
```
src/mlip_autopipec/
├── __init__.py
├── cli.py                # Main CLI entry point (using Click)
├── config/
│   ├── __init__.py
│   └── models.py         # Pydantic models for configuration
├── database/
│   └── ase_db.py         # Wrapper for ASE database interactions
├── generators/
│   ├── __init__.py
│   ├── base.py           # Abstract base class for generators
│   └── alloy.py          # Concrete implementation for alloys
├── explorers/
│   ├── __init__.py
│   └── md_engine.py      # MD and MD/MC simulation engine
├── samplers/
│   ├── __init__.py
│   ├── base.py           # Abstract base class for samplers
│   └── fps.py            # Farthest Point Sampling implementation
├── pipeline/
│   └── orchestrator.py   # The main orchestrator (PipelineRunner)
└── shared/
    └── physics.py        # Physics validation functions
```

**Class/Function Definitions Overview:**
- **`cli.py`**: Will contain the main `click` command group and subcommands for running the pipeline and managing the UI.
- **`config/models.py`**: This file is central to the project's robustness. It will define a set of Pydantic models (`FullConfig`, `SystemConfig`, `ExplorationConfig`, etc.) that parse and validate the user's YAML configuration file. This ensures that all parameters are of the correct type and within valid ranges before any computation begins.
- **`database/ase_db.py`**: The `AseDBWrapper` class will encapsulate all interactions with the ASE database, providing methods to connect, write, and read atomic structures.
- **`generators/base.py`**: An abstract base class `BaseStructureGenerator` will define the common interface for all generators (e.g., a `generate()` method).
- **`pipeline/orchestrator.py`**: The `WorkflowOrchestrator` class will be the heart of the application. It will take the validated configuration object, instantiate the necessary components (generator, explorer, etc.) using factories, and execute the four pipeline stages in sequence.
- **`shared/physics.py`**: This module will contain pure functions for performing physics-based validation, such as `check_interatomic_distance`.

This structure ensures that components are loosely coupled and highly cohesive. For example, a `Generator` does not need to know about the database; its only job is to produce a list of `ase.Atoms` objects. The `WorkflowOrchestrator` is responsible for handling the data flow between these components.

## 5. Implementation Plan

The project will be developed over two distinct cycles, allowing for an iterative and manageable implementation process.

**Cycle 1: Core CLI Pipeline & Foundational Components**
This cycle focuses on building the essential, non-interactive backbone of the application. The primary goal is to have a functional command-line tool that can execute the entire data generation pipeline from start to finish for a simple case (e.g., an alloy system). Key deliverables for this cycle include:
-   **Configuration Management:** Implement the full Pydantic model structure in `config/models.py` to parse and validate a user's YAML configuration file using Hydra.
-   **CLI Entry Point:** Create the initial CLI using `click` in `cli.py`, with a single command `run-pipeline` that accepts the path to the configuration file.
-   **Structure Generation:** Implement the `BaseStructureGenerator` and a concrete `AlloyGenerator`. This generator will be capable of creating randomized alloy structures in a supercell.
-   **Basic Exploration Engine:** Implement a simplified MD engine in `explorers/md_engine.py` that can run simulations using a basic potential like EMT (Effective Medium Theory) provided by ASE. The advanced features (hybrid MD/MC, potential mixing) will be deferred to Cycle 2.
-   **Simple Sampling:** Implement a basic `RandomSampler` that selects a specified number of frames randomly from the MD trajectory.
-   **Database Storage:** Implement the `AseDBWrapper` to connect to a SQLite database and save the final sampled structures.
-   **Orchestration:** Develop the `WorkflowOrchestrator` to tie all the above components together and execute the pipeline sequentially.

**Cycle 2: Advanced Features, Web UI & Optimisation**
This cycle builds upon the foundation of Cycle 1, introducing the more complex scientific algorithms and a user-friendly interface that distinguish the project. The goal is to enhance the quality and diversity of the generated data and improve the user experience.
-   **Advanced Exploration Engine:** Enhance the `md_engine.py` to support the hybrid MD/MC scheme. This includes implementing Monte Carlo moves like atom swaps and vacancy hops. It will also include the logic for mixing MLIPs with the ZBL potential for accurate short-range repulsion.
-   **Advanced Sampling:** Implement the Farthest Point Sampling (FPS) algorithm in `samplers/fps.py`. This requires integrating a SOAP descriptor calculator to characterize the local atomic environment of each atom.
-   **Parallel Processing:** Refactor the `WorkflowOrchestrator` and `md_engine.py` to use `ProcessPoolExecutor` for running multiple MD simulations in parallel, significantly speeding up the exploration stage. This will include implementing the "late binding" calculator pattern to avoid pickling issues with large ML models.
-   **Web User Interface:** Develop a simple web application (using a framework like Streamlit or Flask) that allows a user to:
    -   Build the configuration file through an interactive form.
    -   Launch a pipeline run.
    -   Monitor the progress of the run.
    -   Visualize the structures in the final database.
-   **Expanded Generators:** Implement additional generators for other physical systems, such as `IonicGenerator`, which will include charge-balancing logic.

## 6. Test Strategy

Testing will be a critical component of each cycle to ensure the correctness and robustness of the framework.

**Cycle 1 Test Strategy:**
-   **Unit Testing:** Each component will be tested in isolation.
    -   **Configuration:** Tests will ensure the Pydantic models correctly validate good configurations and raise `ValidationError` for bad ones (e.g., negative temperatures, invalid element symbols).
    -   **Generators:** The `AlloyGenerator` will be tested to confirm that it produces the correct number of structures, with the correct atomic composition and within plausible physical constraints. Mocks will be used for any external dependencies.
    -   **Database:** The `AseDBWrapper` will be tested to verify that it can correctly connect to, write to, and read from a temporary test database.
-   **Integration Testing:** An end-to-end test for the `run-pipeline` CLI command will be created. This test will use a minimal, non-trivial configuration file. It will run the entire pipeline, mocking the computationally expensive MD simulation (`md_engine.py`) to return a pre-defined, small trajectory. The test will then assert that the final ASE database contains the expected number of structures, correctly sampled from the mock trajectory. This verifies that all components are correctly wired together and that data flows between them as expected.

**Cycle 2 Test Strategy:**
-   **Unit Testing:**
    -   **Advanced Explorer:** The logic for the hybrid MD/MC moves will be tested to ensure atom swaps are performed correctly and that charge safety mechanisms work as intended.
    -   **FPS Sampler:** The FPS algorithm will be tested with a small, known set of vectors to ensure it correctly selects the most diverse subset. The SOAP descriptor integration will be mocked.
    -   **Web UI:** The backend logic of the web application (form handling, configuration generation) will be unit tested.
-   **Integration Testing:**
    -   The integration test from Cycle 1 will be expanded to cover the new features. A separate test case will be added to invoke the pipeline using the FPS sampler and assert the correctness of the final selection.
    -   Tests for the parallel execution logic will be added to ensure that the `ProcessPoolExecutor` is called correctly and that results from multiple processes are aggregated properly.
-   **End-to-End (E2E) Testing:**
    -   A simple E2E test for the Web UI will be implemented using a tool like Playwright. The test will simulate a user interacting with the web form, submitting a job, and verifying that the resulting configuration is correct. This ensures the frontend and backend are correctly integrated.
