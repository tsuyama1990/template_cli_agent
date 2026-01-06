# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) is a sophisticated, high-throughput software framework designed to automate the generation of high-quality, physically valid, and diverse training data for modern Machine Learning Interatomic Potentials (MLIPs), such as MACE and SevenNet. The core philosophy of the project is to fundamentally address a critical bottleneck in the materials science workflow: the creation of robust and comprehensive training datasets. Historically, this process has been a manual, time-consuming, and often ad-hoc endeavor, heavily reliant on the intuition and experience of domain experts. MLIP-AutoPipe aims to replace this paradigm with a robust, reproducible, and intelligent pipeline that removes the human expert from the loop. The system is engineered to explore the thermodynamic phase space of a given chemical system efficiently and systematically. Instead of merely generating random atomic configurations, which often results in physically impossible or irrelevant structures, it leverages proven simulation techniques like molecular dynamics (MD) and Monte Carlo (MC) to discover thermodynamically accessible structures. This includes not only stable, low-energy configurations but also critical high-energy configurations and transition states. These "hard" examples are often difficult to sample but are absolutely essential for training MLIPs that are accurate, robust, and can generalize well outside of their initial training distribution.

The framework is designed from the ground up to be a versatile tool for a wide array of materials science research. Its capabilities include the generation of structures for complex multi-component alloys, ionic crystals with charge-balance constraints, covalent materials, technologically important material interfaces, and surface adsorption systems for catalysis research. A particularly powerful feature is its knowledge-based structure generation, which can infer and create plausible initial crystal structures from a simple chemical formula (e.g., Fe3Pt) by intelligently utilizing crystallographic databases and symmetry principles. This allows researchers to rapidly bootstrap their simulations without needing to manually build initial structures. The entire workflow is orchestrated by a central `PipelineRunner`, which executes a clear, four-stage process: Generation, Exploration, Sampling, and Storage. This modular architecture is a deliberate design choice that ensures a strong separation of concerns, which in turn enhances maintainability, testability, and allows for future extensibility. Each stage is designed with state isolation in mind; intermediate results are aggressively checkpointed to a database and the filesystem. This design is critical for the long-running simulations common in this field, preventing catastrophic data loss in the event of a system crash and allowing the pipeline to be cleanly resumed from the last successful stage. Performance is also a key consideration, and the system is designed to leverage modern multi-core architectures by parallelizing the computationally intensive exploration phase.

## 2. System Design Objectives

The primary objective of MLIP-AutoPipe is to provide a fully automated, reliable, and efficient pipeline for MLIP training data generation that is accessible to both computational experts and experimentalists. The key goals, constraints, and success criteria are detailed below:

**Goals:**
*   **Full Automation:** The system aims to minimize, and ideally eliminate, the need for manual intervention in the data generation process. The user's interaction should be limited to providing a high-level configuration file that specifies the chemical system, desired physical conditions, and computational parameters. The pipeline should handle everything else, from initial structure creation to final database archival, without requiring intermediate user input. This frees up valuable researcher time for data analysis and scientific discovery, rather than tedious data preparation.
*   **High-Quality, Physically Valid Data:** The generation of physically nonsensical structures is a common pitfall of naive data generation methods. MLIP-AutoPipe will enforce strict physical validity at all stages. This includes enforcing minimum interatomic distances to prevent atomic overlap, ensuring correct charge balance in ionic systems, and maintaining simulation cell dimensions that are compatible with the potential's cutoff radius under periodic boundary conditions to avoid self-interaction artifacts. The goal is to produce data that is as close to experimental reality as possible.
*   **Diverse and Informative Data:** A key goal is to produce a rich dataset that covers a wide and relevant area of the potential energy surface. Data diversity is far more important than sheer data volume for training robust MLIPs. This objective is achieved by combining the complementary strengths of Molecular Dynamics (for exploring local energy basins) and Monte Carlo methods (for making large-scale, topology-altering moves). This hybrid approach allows the system to escape kinetic traps and discover a much broader range of configurations. Furthermore, the use of advanced sampling techniques like Farthest Point Sampling (FPS) ensures that the final dataset is not redundant and is rich in structural information.
*   **Computational Robustness and Resilience:** Computational materials science simulations can run for hours or even days. The system must therefore be resilient to failures. This includes robust error handling within the simulation workers (e.g., catching and logging "Coulomb explosion" events without crashing the entire pipeline), automatic recovery mechanisms where possible, and the progressive saving of simulation trajectories to disk to prevent data loss. The pipeline's state will be managed transactionally, so that a failed run can be cleanly audited and resumed.
*   **Modularity and Extensibility:** The architecture must be modular to foster innovation and adaptation. It should be straightforward for other developers to add new structure generators for different material classes (e.g., polymers, 2D materials), new exploration algorithms (e.g., basin hopping, metadynamics), new sampling methods (e.g., clustering-based approaches), and new MLIP models. The use of factory patterns and abstract base classes will be central to achieving this goal by enforcing a common, well-defined interface for all pluggable components.
*   **Accessibility and Usability:** The tool will offer two primary user interfaces to cater to different user profiles. A powerful command-line interface (CLI) will be the primary entry point for batch processing, scripting, and integration into larger high-throughput computational workflows. Complementing the CLI, a web-based graphical user interface (GUI) will be provided for interactive use, allowing for easy parameter tuning, real-time monitoring of running jobs, and immediate visualization of the generated structures. This dual-interface approach lowers the barrier to entry and makes the tool useful for a wider scientific audience.

**Constraints:**
*   The system will be built upon a foundation of well-established, high-quality open-source libraries, including ASE (Atomic Simulation Environment) for core atomic data structures and I/O, Pydantic for rigorous data validation and configuration management, and Hydra for flexible configuration overrides. This avoids reinventing the wheel and leverages the strength of the existing scientific Python ecosystem.
*   The initial implementation will provide out-of-the-box support for widely used MLIPs such as MACE and SevenNet, alongside ASE's built-in EMT potential for fast testing and demonstration purposes. The design will ensure a clear, documented path for adding other MLIP calculators.
*   The computationally intensive exploration stage must be designed for parallel execution on multi-core CPU architectures, as this is the most common hardware available in academic and research environments.

**Success Criteria:**
*   The final system can successfully generate a database of at least 1,000 unique, physically valid structures for a moderately complex ternary alloy system within a reasonable timeframe (e.g., under 12 hours) on a standard 16-core workstation.
*   The generated dataset, when used to train a standard MACE model, results in a potential with demonstrably lower prediction errors for forces and energies (at least a 15% reduction in RMSE) compared to a model trained on a dataset of purely random structures of the same size. This will be the ultimate validation of the "intelligent" data generation approach.
*   The CLI and Web UI are intuitive and well-documented. A new graduate student with basic knowledge of materials science should be able to successfully install the software and run a complete data generation pipeline within 30 minutes, relying only on the provided documentation.

## 3. System Architecture

The MLIP-AutoPipe framework is designed as a modular, four-stage pipeline. A central orchestrator, the `PipelineRunner`, manages the flow of data between these stages. This orchestrator is responsible for initializing each stage, passing the required data, and handling stage transitions. To ensure data integrity, facilitate checkpointing, and allow for easy post-hoc analysis, all persistent data—including initial structures, final datasets, and metadata—are stored in an ASE database, which uses an SQLite backend for portability and ease of use. This centralized data store acts as the single source of truth for the pipeline's state.

```mermaid
graph TD
    subgraph User Input
        A[Config File .yaml]
    end

    subgraph Pipeline Orchestration
        B[PipelineRunner]
    end

    subgraph Stage 1: Generation
        C[GeneratorFactory]
        D{Select Generator}
        E[AlloyGenerator]
        F[IonicGenerator]
        G[KnowledgeGenerator]
    end

    subgraph Stage 2: Exploration
        H[ExplorationEngine]
        I[Parallel MD/MC Runner]
        J[Hybrid MD/MC Logic]
        K[MLIP/ZBL Potential]
    end

    subgraph Stage 3: Sampling
        L[SamplingEngine]
        M{Select Sampler}
        N[RandomSampler]
        O[FPSSampler]
    end

    subgraph Stage 4: Storage
        P[StorageEngine]
    end

    subgraph Data Store
        Q[ASE Database .db]
        R[File System .xyz]
    end

    A --> B
    B --> C
    C --> D
    D -- Alloy --> E
    D -- Ionic --> F
    D -- Knowledge --> G
    E -- writes initial structures --> Q
    F -- writes initial structures --> Q
    G -- writes initial structures --> Q

    B -- triggers --> H
    H -- reads initial structures <-- Q
    H --> I
    I --> J
    J --> K
    K -- writes raw trajectory --> R
    R -- is read by --> L

    B -- triggers --> L
    L --> M
    M -- Random --> N
    M -- FPS --> O
    N -- provides sampled structures --> P
    O -- provides sampled structures --> P

    P -- writes final dataset --> Q

    classDef stage fill:#f9f,stroke:#333,stroke-width:2px;
    class C,D,E,F,G stage;
    class H,I,J,K stage;
    class L,M,N,O stage;
    class P stage;
```

**Data Flow Explained:**
1.  **Configuration:** The process begins with the user providing a YAML configuration file. This file is the complete blueprint for the desired pipeline run, specifying everything from the chemical elements to the simulation temperature and the sampling algorithm.
2.  **Orchestration:** The `PipelineRunner`, the system's central nervous system, parses and validates this configuration file using Pydantic models. This ensures that all parameters are valid before any computation begins.
3.  **Generation:** The `PipelineRunner` instructs the `GeneratorFactory` to create the appropriate generator (e.g., `AlloyGenerator`) based on the configuration. This generator produces a set of initial seed structures. These are not just random atoms; they are physically plausible configurations that serve as the starting points for exploration. These seed structures are immediately saved to the ASE database, marking the first successful checkpoint.
4.  **Exploration:** The `ExplorationEngine` reads the seed structures from the database. It then uses a `ProcessPoolExecutor` to farm out these structures to multiple worker processes, where the MD/MC simulations are run in parallel. Each simulation evolves the initial structure, exploring the nearby phase space and generating a trajectory of thousands of atomic configurations. A key feature here is the automatic ensemble switching (NVT for surfaces to prevent vacuum collapse, NPT for bulk to allow cell relaxation), a piece of embedded expertise that prevents common simulation artifacts. The raw trajectories are continuously streamed to `.xyz` files on disk to avoid holding large datasets in memory.
5.  **Sampling:** Once the exploration phase is complete, the `PipelineRunner` triggers the `SamplingEngine`. This engine reads the raw trajectory files from disk. Based on the user's configuration, it selects a sampling strategy via the `SamplerFactory`. It then processes the trajectories to extract a smaller, more diverse, and more informative subset of structures.
6.  **Storage:** The `StorageEngine` takes the final, curated list of sampled structures and permanently archives them in the ASE database. These structures are stored along with relevant metadata, such as their calculated potential energy, atomic forces, and the simulation parameters from which they originated, ensuring the final dataset is self-contained and reproducible.

## 4. Design Architecture

The software will be designed with a clear separation of concerns, organized into a Python package named `mlip_autopipec`. The file structure is designed to be modular and extensible, reflecting the pipeline's stages. This clean architecture makes the codebase easier to understand, maintain, and extend.

**File Structure:**
```
.
├── pyproject.toml
├── src
│   └── mlip_autopipec
│       ├── __init__.py
│       ├── cli.py                # Main CLI entry point (Typer)
│       ├── config.py             # Pydantic models for configuration
│       ├── database.py           # ASE DB wrapper for safe transactions
│       ├── main_gui.py           # Web UI entry point (Streamlit/Flask)
│       │
│       ├── pipeline
│       │   ├── __init__.py
│       │   └── runner.py         # PipelineRunner orchestrator
│       │
│       ├── generators
│       │   ├── __init__.py
│       │   ├── base.py           # BaseGenerator abstract class
│       │   ├── alloy.py          # AlloyGenerator implementation
│       │   ├── ionic.py          # IonicGenerator implementation
│       │   └── factory.py        # GeneratorFactory to select generator
│       │
│       ├── exploration
│       │   ├── __init__.py
│       │   └── engine.py         # ExplorationEngine, parallel MD/MC logic
│       │
│       ├── sampling
│       │   ├── __init__.py
│       │   ├── base.py           # BaseSampler abstract class
│       │   ├── random.py         # RandomSampler implementation
│       │   ├── fps.py            # FPSSampler implementation
│       │   └── factory.py        # SamplerFactory to select sampler
│       │
│       └── storage
│           ├── __init__.py
│           └── engine.py         # StorageEngine for database archival
│
└── tests
    ├── __init__.py
    ├── test_cli.py
    └── ... (a dedicated test file for each source module)
```

**Data Models (Pydantic):**
The `config.py` file is the most critical piece of the design architecture. It defines the "schema" for the entire application through a hierarchy of Pydantic models. This schema-first approach is a core principle of the project. It allows us to leverage Pydantic's powerful validation capabilities to catch configuration errors before any expensive computation is performed. This makes the system more robust and user-friendly.

```python
# Example Pydantic Models in config.py
from pydantic import BaseModel, Field, PositiveInt, model_validator

class SystemConfig(BaseModel):
    generator_type: str
    elements: list[str]
    composition: dict[str, float]
    num_structures: PositiveInt

    @model_validator(mode='after')
    def check_composition_sum(self) -> 'SystemConfig':
        if abs(sum(self.composition.values()) - 1.0) > 1e-6:
            raise ValueError("Composition values must sum to 1.0")
        return self

class ExplorationConfig(BaseModel):
    temperature_k: float = Field(gt=0)
    pressure_gpa: float | None = None
    timestep_fs: float = Field(gt=0, default=1.0)
    mc_swap_probability: float = Field(ge=0, le=1, default=0.0)

class FullConfig(BaseModel):
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: dict # To be defined with its own model
    database_path: str
```
This schema serves as a contract for all data that flows into the system, ensuring its integrity. It is also self-documenting; the Pydantic models provide a clear and unambiguous definition of all possible configuration parameters.

## 5. Implementation Plan

The project will be developed over two distinct, sequential cycles. This iterative approach allows us to build and stabilize the core functionality first before adding more complex features.

**CYCLE 01: Core Pipeline and Foundation**
This initial cycle is focused on building the essential scaffolding of the application. The goal is to deliver a functional, end-to-end command-line pipeline that can perform the basic workflow: generate initial structures based on a config file and store them in a database. The more complex exploration and sampling logic will be deferred to the next cycle. This approach mitigates risk by allowing us to validate the core architecture early on.
*   **Features & Sub-tasks:**
    *   **Project Scaffolding:** Set up the project structure with `pyproject.toml`, define dependencies, and create the initial directory layout as specified in the architecture.
    *   **Pydantic Configuration:** Implement the complete set of Pydantic models in `config.py`. This includes defining all fields, types, and custom validators for the `system` and global configurations. This is the first implementation task.
    *   **Database Wrapper:** Implement the `AseDBWrapper` class in `database.py`. This class will abstract all `ase.db` interactions, providing clean, high-level methods like `add_structures` and `get_structures_by_group`. It will handle connection management and transaction safety.
    *   **CLI Entry Point:** Implement the main CLI in `cli.py` using Typer. This will include defining the `run-pipeline` command, handling the `--config-path` argument, and adding user-friendly output using the `rich` library.
    *   **Core Pipeline Runner:** Implement the `PipelineRunner` in `runner.py` with a simplified, linear workflow for this cycle (Generation -> Storage). It will be responsible for orchestrating the interactions between the other components.
    *   **Generator Implementation:** Implement the `BaseGenerator` abstract class and a simple but functional `AlloyGenerator` that can create random alloy structures based on a given composition. The associated `GeneratorFactory` will also be implemented.
    *   **Storage Engine:** Implement the `StorageEngine` in `storage/engine.py`. For this cycle, its responsibility will be to take the list of generated structures and pass them to the `AseDBWrapper` for saving.
    *   **Unit & Integration Tests:** Write comprehensive unit tests for every component (config models, database wrapper, generator) and at least one integration test that runs the simplified CLI pipeline from end to end.

**CYCLE 02: Advanced Exploration, Sampling, and UI**
This cycle builds directly upon the stable foundation of Cycle 01. It introduces the sophisticated, computationally intensive parts of the pipeline and the user-facing graphical interface, transforming the basic tool into a powerful research platform.
*   **Features & Sub-tasks:**
    *   **Exploration Engine:** Implement the `ExplorationEngine` in `exploration/engine.py`. This is the most complex backend task. It involves implementing the parallel MD execution logic using `ProcessPoolExecutor`. A key sub-task is to develop the `run_single_md_process` worker function, which will encapsulate the logic for a single simulation, including the "late-binding" of the ASE calculator to avoid pickling errors. The hybrid MD/MC logic, including atom swaps and vacancy hops, will also be implemented here.
    *   **Sampling Engine:** Implement the `FPSSampler` in `sampling/fps.py`. This requires integrating a library like `dscribe` to compute the SOAP descriptors that serve as the input for the FPS algorithm. The `BaseSampler` and `SamplerFactory` will also be created to ensure a modular design.
    *   **Pipeline Integration:** The main `PipelineRunner` will be modified to incorporate the `ExplorationEngine` and `SamplingEngine`, expanding the workflow to the full four stages. This involves managing the flow of data (e.g., lists of trajectory file paths) between the new stages.
    *   **Web UI Development:** Develop the Web UI in `main_gui.py` using Streamlit. The UI will be designed to allow users to intuitively configure all pipeline parameters using interactive widgets (sliders, dropdowns, text inputs). A key feature will be a "Run Pipeline" button that executes the `PipelineRunner` in a background process. The UI will provide real-time feedback and display the final, visualized structures.
    *   **Advanced Testing:** Write new unit tests for the complex logic in the exploration and sampling engines, using mocking extensively. Develop a comprehensive integration test that runs the full four-stage pipeline on a small test case to verify that all components work together correctly.

## 6. Test Strategy

Testing will be a continuous and integral part of each development cycle to ensure the reliability, correctness, and scientific validity of the code. A multi-layered testing approach will be used.

**CYCLE 01 Test Strategy:**
*   **Unit Testing:** Unit tests will be written using `pytest` and will focus on testing each module in complete isolation.
    *   **Configuration (`test_config.py`):** The Pydantic models will be tested exhaustively. This will include tests for the "happy path" with valid configurations, as well as numerous failure cases to ensure that edge cases (e.g., zero values, empty lists) and invalid configurations (e.g., compositions not summing to 1.0, misspelled keys) raise the expected `ValidationError` with informative messages.
    *   **Database (`test_database.py`):** The `ase.db.connect` function will be mocked using `pytest-mock`. This allows us to test the logic of our `AseDBWrapper` (e.g., that it calls the underlying `connection.write()` method with the correct arguments) without any actual disk I/O, making the tests fast and reliable.
    *   **Generator (`test_generators.py`):** The `AlloyGenerator` will be tested to verify that its output is correct. We will assert that it produces the correct number of structures, that the chemical composition of the generated structures matches the request, and, critically, that the internal physical validation checks (e.g., `overlap_check`) correctly identify and reject invalid structures.
*   **Integration Testing:** A small number of integration tests will be written to verify that the core components work together as expected.
    *   The primary integration test will use the `Typer.testing.CliRunner` to invoke the main CLI command. This test will create a temporary configuration file and a temporary database path, run the minimal pipeline (Generation -> Storage), and then assert that the database file is created and contains the correct number of structures. This validates the entire data flow for Cycle 1.

**CYCLE 02 Test Strategy:**
*   **Unit Testing:** The focus will be on the new, complex components.
    *   **Exploration Engine (`test_exploration.py`):** This component is computationally expensive and has external dependencies, making it a prime candidate for mocking. We will not run actual MD simulations in the unit tests. Instead, the ASE dynamics runner (`dyn.run()`) will be mocked. Tests will focus on verifying the correctness of the setup logic: for a given configuration, was the ASE dynamics object initialized with the correct temperature, pressure, and ensemble? We will also test the parallelization logic by asserting that the mocked worker function is called the expected number of times.
    *   **Sampling (`test_sampling.py`):** The `FPSSampler` will be tested with a small, deterministic set of input structures for which the diversity is known. The SOAP descriptor calculation will be mocked to return pre-computed feature vectors. The test will then assert that the FPS algorithm correctly selects the most diverse subset of structures based on these known feature vectors, proving the correctness of the implementation.
*   **Integration Testing:** The integration tests will be expanded to cover the full, four-stage pipeline.
    *   A full end-to-end integration test will be created. It will be carefully configured to run on a very small system (e.g., a 4-atom cell) for a minimal number of MD steps using a fast potential like EMT. The test will run the entire pipeline (Generation -> Exploration -> Sampling -> Storage) and assert the final state of the database. This test is essential for finding bugs in the data contracts and interactions between the different pipeline stages. It will be marked as a "slow" test and run less frequently than the unit tests (e.g., only on pull requests to the main branch).
*   **UI Testing:**
    *   The Web UI will be manually tested across major browsers to ensure all widgets are responsive and correctly update the configuration. Given the complexity of setting up automated UI tests, the initial focus will be on manual testing, with a plan to potentially add automated end-to-end tests using a tool like Playwright in the future if the project grows in complexity.
