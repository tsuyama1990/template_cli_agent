# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) is a system designed to fully automate the creation and validation of high-fidelity interatomic potentials. These potentials are mathematical functions that describe the energy of a system of atoms, and they are crucial for performing large-scale atomistic simulations in materials science. The primary goal of this project is to eliminate the human expert from the loop, a concept we refer to as "removing the human expert from the loop." Traditionally, creating a reliable MLIP requires deep domain knowledge and intuition. A researcher must decide which atomic configurations to generate for the training dataset, what temperature and pressure conditions to simulate, and how to assess the quality of the resulting potential. This process is often slow, expensive, and difficult to reproduce because it relies heavily on the implicit knowledge of an experienced scientist.

MLIP-AutoPipe replaces this manual, intuition-driven process with a systematic and autonomous workflow. The system's core philosophy is to combine physics-based heuristics with modern machine learning techniques, such as uncertainty quantification. This dual approach allows the system to make intelligent decisions about what data to generate and learn from. The pipeline starts with minimal user input: typically just the chemical elements and their composition. From this, it autonomously generates a diverse set of initial atomic structures, calculates their properties using high-precision quantum mechanics simulations (specifically, Density Functional Theory or DFT), trains an initial MLIP, and then iteratively improves it through a process called active learning. The active learning loop involves using the current MLIP to run molecular dynamics (MD) or kinetic Monte Carlo (kMC) simulations, monitoring for areas where the model is uncertain, and then generating new DFT data for those specific configurations to retrain and improve the model. This ensures that computational effort is focused only on the most informative data points, dramatically increasing efficiency. The final output is a robust MLIP that can predict energies and forces with near-DFT accuracy but at a fraction of the computational cost, enabling simulations of complex phenomena over long timescales, such as phase transformations, chemical reactions, and material degradation.

## 2. System Design Objectives

The design of the MLIP-AutoPipe system is guided by several key objectives aimed at making it powerful, user-friendly, and future-proof.

**Primary Objectives:**
- **Autonomy:** The highest priority is to create a fully autonomous pipeline. The system must be able to proceed from initial structure generation to final potential validation without requiring any manual intervention. This involves automating DFT parameter selection, convergence checks, error recovery, and the active learning feedback loop.
- **Accuracy:** The generated MLIP must achieve an accuracy level comparable to the underlying DFT calculations. This means that the forces, energies, and stresses predicted by the MLIP should have a very low error margin relative to the ground truth data.
- **Efficiency:** The pipeline must be computationally efficient. This is achieved by minimising the number of expensive DFT calculations. The system employs a "foundation model" surrogate for initial exploration and uses uncertainty quantification to intelligently select new data points, ensuring that computational resources are spent wisely.
- **Reproducibility:** Every step of the workflow must be logged and version-controlled. The system will use a two-tier configuration strategy, where a minimal user input file is expanded into a complete, detailed execution configuration. This full configuration file is saved, ensuring that any result can be perfectly reproduced at a later date.

**Secondary Objectives:**
- **Modularity and Extensibility:** The system is designed as a series of loosely coupled modules. This modular architecture makes it easy to update, replace, or add new algorithms in the future. For example, a new structure generation method or a different DFT code could be integrated with minimal effort.
- **Performance:** Key computational bottlenecks, particularly in descriptor calculation and data processing, will be optimised using high-performance computing techniques. We will use tools like Numba for Just-In-Time (JIT) compilation of Python code to achieve near-native performance for critical loops.
- **Usability:** Although the system is autonomous, it must be easy for a user to set up and run. The initial configuration will be minimal and intuitive, requiring only the essential chemical information. The system will handle the complex details internally.

**Success Criteria:**
The success of the MLIP-AutoPipe system will be measured by the following criteria:
- The system can successfully generate a validated MLIP for a variety of material types (alloys, molecules, ionic crystals) starting from only a chemical formula.
- The root-mean-square error (RMSE) of the final MLIP's force predictions is below a predefined threshold (e.g., < 0.05 eV/Å) when compared against a held-out DFT test set.
- The total number of DFT calculations required is significantly less (e.g., by an order of magnitude) than a traditional, non-autonomous approach for a comparable system.
- The entire workflow, from input to output, is completed without fatal errors or the need for user intervention.

## 3. System Architecture

The MLIP-AutoPipe system is designed as a modular pipeline orchestrated by a central workflow manager. The architecture is built around five core modules, each responsible for a specific stage of the MLIP generation process. Data is passed between modules, and its state is tracked in a central database to ensure traceability.

```mermaid
graph TD
    A[User Input: input.yaml] --> B{Config Expander};
    B --> C[Full Config: exec_config_dump.yaml];
    C --> D[Module A: Structure Generator];
    D --> E[Initial Structures];
    E --> F{Module B: Explorer & Sampler};
    F --> G[Candidate Structures];
    G --> H{Module C: Labeling Engine (DFT)};
    H --> I[Labeled Data (Energy, Forces)];
    I --> J{Module D: Training Engine};
    J --> K[MLIP];
    K --> L{Module E: Simulation Engine};
    L --> M{Uncertainty Detection};
    M -- Uncertain --> G;
    M -- Certain --> N[Validated MLIP];

    subgraph Workflow Orchestrator
        B; C; D; F; H; J; L; M;
    end

    subgraph Data Store
        E; G; I; K; N;
    end
```

**Component Descriptions:**

1.  **Config Expander (Heuristic Engine):** This is the entry point of the system. It takes the user's minimal `input.yaml` file and uses a physics-based heuristic engine to generate a complete `exec_config_dump.yaml`. This engine automatically determines optimal parameters for DFT calculations (like cutoff energies), simulation settings (temperature steps), and model choices based on the input chemical system.

2.  **Module A: Structure Generator:** This module creates the initial seed dataset of atomic structures. It does not perform any expensive DFT calculations. Instead, it uses fast, physics-aware algorithms based on the type of material identified by the Config Expander (e.g., Special Quasirandom Structures (SQS) for alloys, Normal Mode Sampling (NMS) for molecules).

3.  **Module B: Explorer & Sampler:** This module takes the initial structures and explores a wider configuration space using a pre-trained, universal "foundation model" (like MACE-MP or M3GNet). This allows for rapid, low-cost exploration. It then uses a technique called DIRECT sampling to select a diverse and informative subset of structures to be passed to the DFT engine, avoiding redundancy and maximising learning efficiency.

4.  **Module C: Labeling Engine:** This is the automated DFT calculation module. It takes candidate structures and computes their energies, forces, and stresses using Quantum Espresso. It handles all aspects of the DFT calculation automatically, including input file generation, job submission, convergence checks, and error recovery.

5.  **Module D: Training Engine:** This module trains the MLIP using the labeled data generated by Module C. It employs a "Delta Learning" approach, where the model learns the difference between a simple, physical baseline potential (like a Lennard-Jones potential) and the true DFT values. This technique improves data efficiency and ensures the model behaves physically even in regions with sparse data.

6.  **Module E: Simulation Engine:** This module uses the trained MLIP to run large-scale simulations (MD or kMC). Its key feature is the "on-the-fly" (OTF) uncertainty quantification. It continuously monitors the model's predictions during the simulation. If it encounters a configuration where it is uncertain, it pauses the simulation, extracts that structure, and sends it back to Module C for DFT labeling. This active learning loop progressively refines the MLIP.

7.  **Data Store:** All artifacts—structures, DFT results, trained models, configurations—are stored and tracked, likely using a database interface like ASE's DB module, ensuring data provenance and traceability.

## 4. Design Architecture

The software will be a Python-based command-line application, built with modern development practices. The project will be managed using `pyproject.toml` and the `uv` package manager for fast and reproducible environments.

**File Structure:**
```
.
├── dev_documents/
│   ├── ALL_SPEC.md
│   ├── SYSTEM_ARCHITECTURE.md
│   └── CYCLE01/
│       ├── SPEC.md
│       └── UAT.md
├── src/
│   └── mlip_autopyke/
│       ├── __init__.py
│       ├── cli.py          # Main entry point (using Typer/Click)
│       ├── config/         # Handles input.yaml parsing and expansion
│       │   ├── expander.py
│       │   └── models.py   # Pydantic models for config files
│       ├── workflow/       # Orchestrator and state management
│       │   └── orchestrator.py
│       ├── modules/        # Core modules (A-E)
│       │   ├── a_structure_generator.py
│       │   ├── b_explorer_sampler.py
│       │   ├── c_labeling_engine.py
│       │   ├── d_training_engine.py
│       │   └── e_simulation_engine.py
│       └── common/         # Shared utilities, data models, etc.
│           ├── atoms_handler.py
│           └── db_interface.py
└── tests/
    ├── conftest.py
    └── test_*.py
```

**Class/Function Overview:**

-   `cli.py`: Will contain the main CLI application logic. A function `run_pipeline(config_path: Path)` will initiate the process.
-   `config/models.py`: Pydantic models (`MinimalConfig`, `FullConfig`) will define the schema for the configuration files, providing validation and type safety.
-   `config/expander.py`: The `ConfigExpander` class will contain the heuristic logic to transform a `MinimalConfig` object into a `FullConfig` object.
-   `workflow/orchestrator.py`: The `WorkflowOrchestrator` class will manage the overall pipeline. It will have methods like `execute_cycle()` that call the different modules in sequence and manage the flow of data.
-   `modules/*.py`: Each module will be implemented as a class (e.g., `StructureGenerator`, `LabelingEngine`). Each class will have a primary public method, e.g., `generate()`, `label()`, `train()`, which takes data from the previous step and returns the processed artifact. This promotes a clean, consistent interface between modules.
-   `common/db_interface.py`: A `DatabaseInterface` class will provide a simple API (e.g., `save_atoms()`, `get_labeled_data()`) to abstract the underlying database operations.

**Data Models:**
Data transfer between modules will use well-defined objects, primarily leveraging the `ase.Atoms` object from the Atomic Simulation Environment (ASE) library, which is a standard for representing atomic structures in Python. DFT results and model predictions will be stored in the `atoms.info` and `atoms.arrays` dictionaries within these objects, ensuring that data and metadata are tightly coupled.

## 5. Implementation Plan

The project is decomposed into five logical cycles, each building upon the last. This iterative approach allows for the incremental delivery of functionality and reduces risk.

**Cycle 1: The Core Engine.** The first cycle focuses on establishing the absolute core of the pipeline: the ability to perform a DFT calculation and train a model from it. This cycle will implement the `Labeling Engine (Module C)` to automate Quantum Espresso calculations and the `Training Engine (Module D)` to train a basic MLIP using the generated data. The main goal is to create a functional, albeit manual, workflow for `Data In -> MLIP Out`. This cycle will also establish the modern Python project structure with `pyproject.toml` and `uv`.

**Cycle 2: Structure Generation and Configuration.** This cycle builds the user-facing entry point of the system. It will implement the `Structure Generator (Module A)` to create initial atomic configurations based on physical heuristics (SQS, NMS, etc.). Crucially, it will also develop the `Config Expander`, which translates a simple `input.yaml` into the detailed `exec_config_dump.yaml`. By the end of this cycle, a user can provide a minimal input file and the system will automatically generate a diverse set of initial structures and the full configuration needed for the subsequent steps.

**Cycle 3: Advanced Sampling and Optimisation.** This cycle focuses on computational efficiency. It will implement the `Explorer & Sampler (Module B)`. This involves integrating a universal foundation model to perform rapid, large-scale explorations of the configuration space without DFT. It will also develop the DIRECT sampling algorithm to intelligently select structures, ensuring maximum information gain for each DFT calculation. Performance bottlenecks identified in previous cycles, especially in data processing and descriptor calculation, will be addressed using Numba for JIT compilation.

**Cycle 4: Active Learning and Advanced Simulation.** This cycle closes the main automation loop. It implements the `Simulation Engine (Module E)`. This module will integrate the MLIP with a simulation code like LAMMPS. The core task is to develop the on-the-fly (OTF) uncertainty monitoring system. This will enable the pipeline to perform active learning: running a simulation, detecting new and uncertain atomic configurations, and automatically sending them back to Module C for labeling and retraining. This makes the system self-improving.

**Cycle 5: User Interface and Finalisation.** The final cycle polishes the system for release. This involves refining the command-line interface (CLI) to make it more user-friendly, improving logging and error reporting, and writing comprehensive user documentation. The focus will be on ensuring the end-to-end workflow is robust, reliable, and easy for a non-expert to use. This cycle also includes final validation and benchmarking of the complete pipeline on several test cases to demonstrate its effectiveness.

## 6. Test Strategy

Testing will be a continuous process throughout all cycles, with a multi-layered approach to ensure the reliability and correctness of the system.

**Unit Testing:** Each function and class will be accompanied by unit tests. We will use the `pytest` framework. For example, in **Cycle 1**, we will write unit tests for the `LabelingEngine` that mock the external Quantum Espresso executable and verify that the correct input files are generated. For the `TrainingEngine`, we will test with mock training data to ensure the model training function completes and returns a model object of the correct format.

**Integration Testing:** We will write tests that verify the interaction between different modules. For **Cycle 2**, an integration test will check that the `StructureGenerator` output can be correctly processed by the `ConfigExpander` and that the resulting configuration is valid. For **Cycle 3**, we will test the full chain from the `Explorer` to the `Sampler` to ensure the data flows correctly and the final selection of structures is in the expected format.

**End-to-End (E2E) Testing:** For each cycle, we will create small, self-contained E2E tests that run the entire pipeline available up to that point on a very simple and fast-to-calculate system (e.g., a silicon dimer). For **Cycle 4**, a key E2E test will be to demonstrate that the active learning loop can successfully execute at least one iteration: run a short MD simulation, trigger an uncertainty event, generate a new DFT calculation, and retrain the model.

**Regression Testing:** All tests will be automated and run via a Continuous Integration (CI) server. This will ensure that new changes do not break existing functionality. We will aim for high test coverage of the codebase.

**Specific Cycle Test Plans:**
-   **Cycle 1:** Test that the `LabelingEngine` correctly parses QE output and that the `TrainingEngine` can fit a model to mock data.
-   **Cycle 2:** Test the `ConfigExpander` with various minimal `input.yaml` files (for alloys, molecules, etc.) and assert that the generated `exec_config_dump.yaml` is complete and correct. Test each method of the `StructureGenerator`.
-   **Cycle 3:** Test the `Explorer` with a mock foundation model. Test the DIRECT sampling implementation to ensure it correctly identifies diverse structures from a large trajectory file.
-   **Cycle 4:** Create a test with a deliberately "bad" initial potential that is guaranteed to trigger an uncertainty event quickly, and verify that the OTF loop executes correctly.
-   **Cycle 5:** E2E tests will be the primary focus, running the full pipeline for several representative materials and checking that the final output is a valid, usable MLIP. We will also test the CLI's error handling for invalid user inputs.
