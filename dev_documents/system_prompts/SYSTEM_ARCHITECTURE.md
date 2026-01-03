# SYSTEM ARCHITECTURE: MLIP-AutoPipe

## 1. Summary

The Machine Learning Interatomic Potential Automatic Pipeline (MLIP-AutoPipe) represents a paradigm shift in computational materials science, moving from a manually intensive, expert-driven process to a fully automated, intelligent, and self-correcting system for generating high-fidelity interatomic potentials. The foundational philosophy of this platform is to **"remove the human expert from the loop"**, directly targeting the primary bottleneck that currently impedes the speed and reproducibility of materials research. The development of Machine Learning Interatomic Potentials (MLIPs) has traditionally been an artisanal craft, where the quality of the final model is inextricably linked to the scientist's intuition and years of experience. Critical decisions—such as selecting a diverse and informative set of training structures, defining appropriate quantum mechanical calculation parameters, and choosing the optimal simulation conditions—rely heavily on heuristics and implicit knowledge. This dependence on human expertise not only creates a steep learning curve for new researchers but also introduces significant variability, making it difficult to systematically compare results across different studies. MLIP-AutoPipe is engineered to dismantle these barriers by codifying this expert knowledge into a robust, autonomous software pipeline.

The system is designed to be a "zero-touch" platform. A user provides the most minimal possible input: the chemical elements of their material (e.g., 'Fe', 'Pt'). From this, the pipeline autonomously orchestrates the entire complex workflow. It intelligently generates a diverse set of initial atomic configurations, submits them for high-precision Density Functional Theory (DFT) calculations, uses the resulting data to train an initial MLIP, and then deploys that MLIP in a closed active learning loop. This loop is the core of the system's intelligence. The MLIP is used to run large-scale simulations, and during these simulations, the system continuously monitors the model's own uncertainty. When the model encounters an atomic environment it is unfamiliar with, it flags that configuration, sends it back for DFT calculation, and incorporates the new information to retrain and improve itself. This iterative refinement process allows the system to build a comprehensive model of the material's potential energy surface, focusing expensive computational resources only on the data points that provide the most significant new information. The system is engineered for maximum modularity and performance, integrating state-of-the-art open-source tools like Quantum Espresso, LAMMPS, and modern ML frameworks such as MACE and ACE. It wraps these powerful tools in a sophisticated orchestration layer that handles automated parameterization, robust error recovery, and intelligent decision-making. Performance-critical algorithms, such as the calculation of atomic descriptors for data sampling, are accelerated using Just-In-Time compilation with Numba, ensuring the pipeline is both intelligent and computationally efficient. The ultimate vision is to create a tool that empowers scientists and engineers to perform high-throughput materials discovery and to simulate complex, long-timescale phenomena like phase transformations and chemical reactions with unprecedented speed, accuracy, and ease of use.

## 2. System Design Objectives

The design of MLIP-AutoPipe is guided by a set of core objectives that collectively aim to create a robust, efficient, and fully autonomous platform for materials potential generation. These objectives address the key pain points of the traditional MLIP development workflow.

**1. Complete Autonomy and the Codification of Expert Knowledge:** The highest-priority objective is the elimination of manual intervention and heuristic decision-making from the user's workflow. The system must encapsulate the knowledge of a domain expert. This means it must be able to automatically analyze a chemical system, determine its fundamental properties (e.g., metallic, covalent, ionic), and select the most appropriate algorithms for each stage of the pipeline. For example, it must know to use Special Quasirandom Structures (SQS) for a disordered alloy but Normal Mode Sampling (NMS) for a molecule. It must automatically select appropriate, high-precision DFT parameters from vetted protocols like SSSP without requiring the user to have any knowledge of pseudopotentials or energy cutoffs. The success of this objective is measured by the system's ability to transform the simplest possible user input (chemical species) into a high-quality, validated potential without requiring any intermediate adjustments or decisions. This objective directly serves the goal of making advanced simulation techniques accessible to a much broader range of scientists and engineers.

**2. Maximisation of Computational Efficiency and Strategic Resource Allocation:** The generation of training data via DFT is, by far, the most computationally expensive part of the process. A central design objective is to achieve the desired model accuracy with the minimum possible number of DFT calculations. MLIP-AutoPipe is designed as a "data-frugal" system. It achieves this through a hierarchical strategy. At the lowest level, it uses a pre-trained, universal surrogate potential (the MACE model) to perform massive, ultra-fast explorations of the material's configuration space. This allows the system to scan millions of potential structures at a tiny fraction of the cost of DFT. From this vast pool, the DIRECT sampling algorithm intelligently selects a small, diverse, and information-rich subset for further analysis. This ensures that the initial training set is not redundant but covers a wide swath of the potential energy surface. Finally, the active learning loop ensures that subsequent DFT calculations are not chosen randomly but are precisely targeted at the specific configurations where the current MLIP is most uncertain. This intelligent, multi-stage data acquisition strategy is designed to concentrate expensive computational resources where they will have the most impact on improving the model.

**3. Unwavering Robustness and Automated Fault Tolerance:** For an autonomous system designed to run for hours or days without supervision, robustness is not a feature but a necessity. The pipeline must be resilient to the myriad of transient and recoverable errors that are common in complex scientific computing workflows. This objective dictates that every interaction with an external program, especially the DFT engine, must be wrapped in a sophisticated error-handling layer. For example, the system will not simply fail if a DFT calculation does not converge. Instead, it will trigger a pre-defined, multi-step recovery protocol, automatically trying a sequence of interventions such as reducing the electronic mixing parameter, increasing the smearing temperature, or even attempting a different initial magnetic configuration. Only after this exhaustive, automated recovery process fails will the system flag the structure as problematic and move on. This built-in fault tolerance is critical to ensuring the pipeline's reliability and preventing entire runs from being wasted due to minor, fixable issues.

**4. Commitment to Physical Fidelity and Predictive Accuracy:** The ultimate measure of the system's success is the quality of the potentials it produces. The final MLIP must be as accurate as the underlying DFT data and must obey fundamental physical laws. The design incorporates several features to ensure this. The "delta learning" approach, where the model learns the difference between DFT and a simpler, classical potential, helps to enforce physically correct repulsive behavior at very short interatomic distances. Furthermore, the active learning loop's fragment-based retraining process includes advanced boundary condition handling. When a high-uncertainty configuration is extracted from a large simulation, the system uses techniques like force masking and hydrogen passivation to ensure that the small, periodic fragment used for the DFT calculation does not contain unphysical artifacts from the artificial boundaries. This meticulous attention to physical detail ensures that the resulting MLIP is not just a statistical fit but a genuinely predictive model of the material's atomic interactions.

**5. Architectural Modularity and Future-Proof Extensibility:** The fields of machine learning and materials science are evolving at a breathtaking pace. A monolithic, inflexible architecture would quickly become obsolete. Therefore, a core design objective is to build the system as a collection of loosely coupled, swappable modules. Each major component—the structure generator, the labeling engine, the training engine—is defined by a clean, abstract interface. The central `WorkflowOrchestrator` interacts with these interfaces, not their concrete implementations. This powerful design pattern means that new technologies can be easily integrated. For instance, if a new, more efficient DFT code becomes available, a new `VaspEngine` class can be written that conforms to the `ILabelingEngine` interface, and it can be dropped into the pipeline without modifying any other part of the system. This modularity ensures the long-term viability of the MLIP-AutoPipe platform, allowing it to adapt and incorporate the best new tools and algorithms as they emerge.

## 3. System Architecture

The MLIP-AutoPipe system is architected as a modular, pipeline-driven application, orchestrated by a central control unit and underpinned by a robust data persistence layer. The design philosophy emphasizes a clear separation of concerns, ensuring that each component has a single, well-defined responsibility. This approach promotes maintainability, testability, and the extensibility outlined in the design objectives. The system's brain is the **Workflow Orchestrator**, a stateful component responsible for managing the sequence of operations, moving data between modules, and controlling the flow of the active learning loop. All data generated and consumed by the workflow, from atomic structures to final DFT results and trained models, is systematically tracked in a database layer abstracted by the **ASE DB Wrapper**. This guarantees provenance and reproducibility for every step of the process.

The logical flow of the application begins when the user provides a minimal `input.yaml` file. This file is immediately consumed by the **Configuration Engine**. This is the first and most critical component in the "expert in a box" design. It contains a sophisticated set of physics-based heuristics that translate the user's high-level request into a comprehensive, low-level execution plan. It automatically determines hundreds of technical parameters—such as DFT k-point densities, plane-wave energy cutoffs, smearing types, and initial magnetic moments—that would typically require an expert's judgment. This expanded and validated configuration is then passed to the Workflow Orchestrator, which commences the main pipeline execution.

The pipeline itself is composed of five core modules, each representing a distinct stage of the potential generation process:

1.  **Module A: Structure Generator:** This module serves as the initial source of atomic configurations. It is designed as a factory that selects the most physically appropriate generation method based on the material type identified by the Configuration Engine. For metallic alloys, it employs the Special Quasirandom Structures (SQS) method to model disorder. For molecular systems, it uses Normal Mode Sampling (NMS) to explore vibrational degrees of freedom. For ionic and covalent crystals, it uses techniques like Ab Initio Random Structure Searching (AIRSS) and deep rattling, respectively. The guiding principle of this module is to create an initial dataset that is both diverse and physically plausible, providing a solid foundation for the subsequent learning stages without the immense cost of ab initio molecular dynamics.

2.  **Module B: Explorer & Sampler:** This module is responsible for the cost-effective exploration of the potential energy surface. It leverages a pre-trained, general-purpose surrogate potential (MACE) that is orders of magnitude faster than DFT. The module runs large-scale molecular dynamics simulations with this surrogate to generate massive trajectories, potentially containing millions of atomic configurations. This vast dataset is then processed by the DIRECT (Diversity Ranking by Interactive Clustering of Trajectories) sampler. The sampler computes a physical descriptor for each configuration and then uses a clustering algorithm to select a small subset of structures that are maximally diverse and representative of the entire trajectory. This ensures that the computationally expensive DFT calculations are reserved for the most informative and unique structural motifs.

3.  **Module C: Labeling Engine:** This component acts as the robust, fault-tolerant interface to the external quantum mechanics code (Quantum Espresso). It receives a set of structures and is responsible for generating the correct input files, executing the DFT calculation, parsing the output to extract the energy, forces, and stresses, and handling any errors that occur. Its design is centered around reliability, incorporating a multi-stage automated error recovery protocol to handle common issues like SCF convergence failure.

4.  **Module D: Training Engine:** This module consumes the labeled data from the database and trains the MLIP. It implements the delta learning strategy, fitting the model to the residual between the DFT data and a baseline physical potential to improve robustness. It also encapsulates the logic for automated hyperparameter optimization, allowing it to fine-tune key model parameters for the specific material under study.

5.  **Module E: Simulation Engine:** This is the heart of the active learning loop. It takes a newly trained MLIP and deploys it in a full-scale simulation. The key function of this engine is to perform "on-the-fly" uncertainty quantification. It continuously monitors the confidence of the MLIP's predictions as the simulation progresses. When the uncertainty exceeds a dynamic threshold, it pauses the simulation, extracts the challenging configuration, and sends it back to the Workflow Orchestrator to be labeled by the Labeling Engine. This feedback loop is the mechanism by which the system autonomously identifies and corrects its own deficiencies.

```mermaid
graph TD
    A[User: input.yaml] --> B(Configuration Engine);
    B --> C{exec_config_dump.yaml};
    C --> D(Workflow Orchestrator);

    subgraph "Autonomous Learning Loop"
        D --> E[Module A: Structure Generator];
        E --> F[Initial Structures];
        F --> G[Module B: Explorer & Sampler];
        G --> H[Candidate Structures];
        H --> I[Module C: Labeling Engine];
        I --> J[Labelled Data (DFT)];
        J --> K[Module D: Training Engine];
        K --> L[Trained MLIP];
        L --> M[Module E: Simulation Engine];
        M --> N{High Uncertainty?};
        N -- Yes --> O[High-Uncertainty Structure];
        O --> I;
        N -- No --> P[Simulation Results];
    end

    subgraph "Data & Model Persistence"
        Q[ASE Database]
        I --> Q;
        K --> Q;
    end

    style B fill:#f9f,stroke:#333,stroke-width:2px
    style D fill:#ccf,stroke:#333,stroke-width:2px
```

## 4. Design Architecture

The software design of MLIP-AutoPipe is rigorously based on modern, best-practice software engineering principles to ensure the system is maintainable, scalable, and robust. The architecture is built upon a foundation of Dependency Inversion and a schema-first design philosophy, utilizing Pydantic for data validation and contract enforcement between components. This approach minimizes tight coupling and allows different parts of the system to be developed and tested independently. The project is structured as a standard Python package, managed by the `uv` toolchain, with a clear separation between application logic, configuration, and data definitions.

**File Structure:**

```
mlip_autopipec/
├── pyproject.toml
├── uv.lock
└── src/
    └── mlip_autopipec/
        ├── __init__.py
        ├── cli.py
        ├── config/
        │   ├── __init__.py
        │   ├── models.py      # Pydantic models for configuration
        │   └── loader.py      # Logic for loading and expanding config
        ├── data/
        │   ├── __init__.py
        │   └── models.py      # Pydantic models for scientific data (Atoms, DFTResults)
        ├── database/
        │   ├── __init__.py
        │   └── ase_db_wrapper.py # Wrapper for ASE DB interactions
        ├── workflows/
        │   ├── __init__.py
        │   └── orchestrator.py  # Main workflow control logic
        ├── interfaces/
        │   ├── __init__.py
        │   └── engines.py     # Abstract Base Classes for all engines
        └── engines/
            ├── __init__.py
            ├── structure_generator.py # Module A
            ├── explorer_sampler.py    # Module B
            ├── labeling_engine.py     # Module C
            ├── training_engine.py     # Module D
            └── simulation_engine.py   # Module E
```

**Core Components Breakdown:**

1.  **Configuration (`config/`):** This package is the single source of truth for all runtime parameters. The `models.py` file defines a series of Pydantic models, most importantly the `MinimalConfig` (reflecting the user's simple input) and the `FullConfig` (the comprehensive, validated configuration object). The `loader.py` file contains the `ConfigExpander` class. This class is responsible for the critical logic of transforming the `MinimalConfig` into the `FullConfig` by applying physical heuristics and populating hundreds of default values. This ensures that the rest of the application deals with a complete, type-safe, and immutable configuration object.

2.  **Data Models (`data/`):** This package defines the Pydantic models for all scientific data that is passed between modules or stored in the database. Key models include `AtomsModel` (representing an atomic structure with symbols, positions, and cell vectors) and `DFTResult` (representing the outputs of a DFT calculation, including energy, forces, and stress). By defining these as Pydantic models, we decouple our internal data representation from any specific external library (like ASE). This allows us to enforce strict validation rules (e.g., ensuring the shape of the forces array matches the number of atoms) and provides a clear, self-documenting contract for the data that flows through the pipeline.

3.  **Interfaces (`interfaces/`):** This is the cornerstone of the system's modularity. Following the Dependency Inversion Principle, this package defines a set of Abstract Base Classes (ABCs) for each of the core engine types (e.g., `ILabelingEngine`, `ITrainingEngine`). The `WorkflowOrchestrator` is programmed against these abstract interfaces, not their concrete implementations. This means the orchestrator's code does not know or care whether it is talking to a `QuantumEspressoEngine` or a hypothetical `VaspEngine`; it only knows it is talking to an object that fulfills the `ILabelingEngine` contract. This decoupling is what makes the system so extensible.

4.  **Engines (`engines/`):** This package contains the concrete implementations of the interfaces defined in the `interfaces` package. Each file (`labeling_engine.py`, `training_engine.py`, etc.) contains a class that inherits from the corresponding ABC and implements its methods. This is where the actual scientific "work" happens: generating input files for external codes, running subprocesses, parsing outputs, and calling machine learning libraries. Computationally intensive algorithms, like the descriptor calculations in `explorer_sampler.py`, will be implemented here and optimized with Numba.

5.  **Workflow Orchestrator (`workflows/orchestrator.py`):** This is the conductor of the entire process. It is initialized with concrete implementations of the engines (an example of Dependency Injection). Its methods contain the high-level logic of the pipeline, such as "for each new structure, call the labeling engine, then save the result." In Cycle 2, its logic will become more complex as it will manage the state of the active learning loop, iteratively calling the simulation and training engines.

6.  **Database Wrapper (`database/ase_db_wrapper.py`):** This class provides a clean, high-level API for all database operations. It completely abstracts the underlying database technology (SQLite via ASE). The rest of the application does not write SQL queries; it calls methods like `db.save_dft_result(structure_id, result_model)`. This isolates database-specific code and makes it much easier to test the components that use the database, as the wrapper can be easily mocked.

7.  **CLI (`cli.py`):** This is the thin user-facing layer of the application. Built with a library like Typer, it is responsible for parsing command-line arguments, loading the `input.yaml` file, instantiating the necessary objects (like the `WorkflowOrchestrator`), and kicking off the pipeline. It also provides user feedback on the progress of the run.

## 5. Implementation Plan

The development of the MLIP-AutoPipe system is strategically divided into two distinct cycles. This phased approach ensures that a solid, functional core is established before the more complex, intelligent features are added. This mitigates risk and allows for early validation of the core architectural decisions.

**CYCLE 01: Core Framework and Automated DFT Pipeline**

This inaugural cycle is focused on building the bedrock of the application. The primary goal is to create a robust, linear, one-shot pipeline that proves the viability of the core automation concept. This cycle corresponds to the first two phases of the project roadmap outlined in `ALL_SPEC.md`. The key deliverable will be a system that can take a minimal configuration file, generate an initial set of structures, automatically calculate their properties using Quantum Espresso, and train a basic MLIP from this data. The emphasis is on architectural soundness and the reliability of the most complex external interface—the `LabelingEngine`.

The implementation will begin by setting up the project structure, defining all the Pydantic data models (`config` and `data`), and establishing the abstract engine interfaces (`interfaces`). This schema-first approach ensures that all subsequent components are built on a solid, contract-driven foundation. The next major task is the implementation of the `QuantumEspressoEngine`, including its sophisticated error-handling and recovery logic. This is the most critical and highest-risk component, so it will be addressed early. Concurrently, the `StructureGenerator` (Module A) and the `TrainingEngine` (Module D) will be developed. The `AseDBWrapper` will be implemented to handle all data persistence needs from the outset. Finally, the `WorkflowOrchestrator` will be created to tie these components together into a single, executable pipeline, which is then exposed to the user via a simple command-line interface in `cli.py`. By the end of this cycle, we will have a functional, end-to-end system, albeit one that lacks the advanced learning capabilities planned for the next phase.

**CYCLE 02: Advanced Sampling, Active Learning, and Performance Optimisation**

The second cycle takes the solid foundation from Cycle 01 and builds the intelligent, autonomous features upon it. This cycle is about closing the loop and transforming the system from a simple automation script into a self-improving learning agent. It corresponds to Phases 3, 4, and 5 of the original roadmap. The two main components to be implemented are the `Explorer & Sampler` (Module B) and the `SimulationEngine` (Module E).

Work in this cycle will commence with the implementation of the `Explorer & Sampler`. This involves integrating the MACE universal surrogate model to run fast, exploratory molecular dynamics simulations. A key part of this task is the development of the performance-critical descriptor calculation code for the DIRECT sampler, which will be heavily optimized with Numba to prevent it from becoming a bottleneck. The second, and most significant, part of this cycle is the development of the `SimulationEngine` and the refactoring of the `WorkflowOrchestrator` to manage the active learning loop. The `SimulationEngine` will be responsible for running simulations with the newly trained MLIP, monitoring its uncertainty, and yielding high-uncertainty structures. The `WorkflowOrchestrator` will be upgraded to handle this new, cyclical workflow: it will receive a structure from the simulation engine, send it to the labeling engine, and then trigger a retraining event with the training engine, thus completing the loop. This cycle will also involve polishing the CLI and providing more detailed feedback to the user about the progress of the active learning generations. At the conclusion of Cycle 02, the MLIP-AutoPipe system will be feature-complete, delivering on the promise of a fully autonomous, "zero-touch" potential generation platform.

## 6. Test Strategy

The test strategy for MLIP-AutoPipe is comprehensive and multi-layered, designed to ensure correctness, reliability, and physical accuracy at every level of the system. The strategy combines rigorous unit testing to validate individual components in isolation, integration testing to ensure they work together correctly, and end-to-end testing to validate the entire workflow on realistic physical systems.

**CYCLE 01: Core Framework and Engine Testing**

During the first cycle, the testing focus is on the foundational components and the linear data pipeline. **Unit testing** will be the primary method of verification. Every class and function will be tested independently of its dependencies using mocks and stubs. The `ConfigExpander` will be a key focus, with tests designed to verify that it correctly infers a wide range of physical parameters from various minimal inputs. The Pydantic data models will be tested to ensure their validation logic is strict and correct. The most intensive unit testing will be applied to the `QuantumEspressoEngine`. We will build a suite of tests that mock the external `subprocess.run` call, allowing us to test the engine's parser against a library of pre-saved, canonical Quantum Espresso output files—including successful runs, various types of convergence failures, and other error conditions. This will allow us to thoroughly validate the parser and the automated error recovery logic without the enormous overhead of actually running DFT calculations.

**Integration testing** in Cycle 01 will focus on verifying the data flow and interactions between the major components. The primary integration test will simulate the entire linear pipeline using a temporary SQLite database file. In this test, the `QuantumEspressoEngine` will be mocked to return predictable results, preventing the test from being slow and flaky. The test will verify that the `WorkflowOrchestrator` correctly calls the `StructureGenerator`, writes the initial structures to the database, passes them to the mocked `LabelingEngine`, correctly receives the mock results, updates the database, and finally passes the complete, collated dataset to the `TrainingEngine`. This test confirms that the "plumbing" of the application is sound, that data is being correctly persisted and retrieved, and that the components are communicating according to their contracts. These tests will provide a high degree of confidence in the core architecture before the more complex, stateful logic of the active learning loop is introduced.

**CYCLE 02: End-to-End and Performance Testing**

In the second cycle, the test strategy will be expanded to cover the new, intelligent components and the full, closed-loop active learning workflow. **Unit testing** will be extended to the new `Explorer & Sampler` and `SimulationEngine` modules. The `SimulationEngine`, in particular, will be tested by mocking the MD simulation itself, allowing us to feed it a pre-determined sequence of atomic structures and mock uncertainty values. This will enable us to test the core logic of the uncertainty thresholding and structure yielding mechanism in a controlled and reproducible way. We will also introduce **performance testing** for the first time. The Numba-optimized descriptor calculation function within the `Explorer & Sampler` will be benchmarked against a large, realistic dataset to verify that its performance meets the required level and does not act as a bottleneck for the entire pipeline.

**Integration testing** for Cycle 02 will be significantly more complex, designed to validate the entire active learning feedback loop. A key test will involve a real `WorkflowOrchestrator` and mocked versions of the `SimulationEngine` and `LabelingEngine`. The mocked `SimulationEngine` will be programmed to yield a specific structure after a few steps. The test will then assert that this exact structure is passed to the `LabelingEngine` and that the result is subsequently used to trigger a retraining event in the `TrainingEngine`. This will validate the data flow through the feedback loop. The ultimate validation will come from **end-to-end (E2E) testing**. These tests will run the *entire*, un-mocked pipeline on a small but physically meaningful system (e.g., a SiGe alloy). It will use the real Quantum Espresso and LAMMPS executables. The final assertion of an E2E test will not be on the internal state of the code, but on the physical properties of the generated potential. For example, the test will conclude by using the final, trained MLIP to calculate a physical property, like the lattice constant of the alloy, and assert that it falls within an acceptable tolerance of the known DFT value. This final layer of testing ensures that the system, as a whole, is not just running correctly but is also producing physically accurate and scientifically valid results.
