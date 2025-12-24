# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) system is a comprehensive, end-to-end platform designed to automate the generation and validation of high-fidelity machine learning interatomic potentials. The core philosophy of the project is to "remove the human expert from the loop," addressing the significant bottlenecks in traditional computational materials science. Currently, the development of MLIPs relies heavily on the intuition and experience of researchers to generate training data, particularly through computationally expensive methods like Ab Initio Molecular Dynamics (AIMD). This dependency on expert knowledge makes the process slow, costly, and difficult to reproduce. MLIP-AutoPipe aims to replace this manual, intuition-driven process with a fully automated, physically-grounded, and robust pipeline.

The system takes minimal user input—typically just a material's composition (e.g., "Fe2O3")—and autonomously executes a complex workflow that spans initial structure generation, active learning, DFT-based data labelling, model training, and on-the-fly simulation. By integrating a suite of modern, high-performance tools and physics-based heuristics, the pipeline intelligently explores the vast chemical and structural space to generate diverse and informative training data. It is designed to produce MLIPs that rival the accuracy of Density Functional Theory (DFT) but at a fraction of the computational cost, enabling simulations on time and length scales previously inaccessible. This automation democratises the creation of bespoke potentials, allowing a wider range of scientists and engineers to investigate complex material phenomena such as phase transitions, defect dynamics, and chemical reactions without requiring deep expertise in quantum chemistry or simulation methodologies. The entire process is designed for maximum efficiency, leveraging modern software practices like `uv` for environment management, `Numba` and `JAX` for computational acceleration, and a modular architecture that ensures future extensibility.

## 2. System Design Objectives

The primary objective of MLIP-AutoPipe is to fully automate the MLIP creation workflow, from initial user input to a validated, simulation-ready potential. This overarching goal is supported by several key design objectives:

*   **Autonomy and Usability**: The system must operate with minimal human intervention. The "Two-Tier Configuration Strategy" is central to this. A user provides a simple `input.yaml` with only the essential information (e.g., elements, composition). The system's "Config Expander" heuristic engine then generates a complete, deterministic `exec_config_dump.yaml`, automatically inferring optimal parameters for DFT calculations, MD simulations, and model training. This lowers the barrier to entry and ensures reproducibility.
*   **Computational Efficiency**: The pipeline must be significantly more efficient than traditional AIMD-based data generation. This is achieved by minimising expensive DFT calculations. A "DIRECT Sampling" approach uses pre-trained, universal foundation models (e.g., MACE-MP, M3GNet) to rapidly explore the potential energy surface. DFT is only used to label a small, information-rich subset of structures identified by this exploration, dramatically reducing the overall computational budget.
*   **Physical Robustness**: The generated potentials must be physically realistic and reliable, even in regions of the configuration space not explicitly covered by DFT data. The "Delta Learning" strategy is key here. The model learns the *difference* between a baseline physical potential (like Lennard-Jones or ZBL) and the true DFT energies. This ensures correct asymptotic behaviour (e.g., core repulsion at short distances) and prevents unphysical predictions.
*   **Self-Improving (Active Learning)**: The system must be capable of iteratively refining its own potential. The "On-the-fly (OTF)" simulation engine (Module E) continuously monitors the model's uncertainty during MD or kMC simulations. When the simulation enters an unknown configuration, the system automatically pauses, extracts the relevant structure, sends it to the DFT Labeling Engine (Module C), and retrains the model (Module D) before resuming. This creates a closed feedback loop for autonomous potential improvement.
*   **Modularity and Extensibility**: The architecture is composed of five loosely-coupled modules (Structure Generator, Explorer & Sampler, Labeling Engine, Training Engine, Simulation Engine). This design allows individual components to be updated or replaced as new algorithms and models become available without requiring a complete system overhaul. For instance, the specific MLIP framework (e.g., ACE, MACE) or DFT code can be swapped out with minimal disruption.
*   **Performance and Scalability**: The entire software stack is built on modern, high-performance tools. `uv` provides rapid and reproducible environments. `Numba` is used to JIT-compile performance-critical Python code (e.g., descriptor calculations, kMC event loops) to achieve C/C++ level speed. Where applicable, `JAX` will be considered for large-scale numerical tasks. The system is designed to leverage multi-core CPUs and GPUs effectively.

## 3. System Architecture

The MLIP-AutoPipe system is orchestrated as a modular pipeline, where data flows sequentially through distinct processing stages. A central Workflow Orchestrator, guided by the `exec_config_dump.yaml`, manages the execution of the five core modules. Data persistence and provenance are handled by a central database (e.g., ASE DB), which tracks every structure, calculation, and model, ensuring full reproducibility.

```mermaid
graph TD
    A[User Input: input.yaml] --> B{Config Expander};
    B --> C[Full Config: exec_config_dump.yaml];
    C --> D[Workflow Orchestrator];

    subgraph Core Modules
        E[A: Structure Generator]
        F[B: Explorer & Sampler]
        G[C: Labeling Engine]
        H[D: Training Engine]
        I[E: Simulation Engine]
    end

    D --> E;
    E --> J[Initial Structures];
    J --> F;
    F --> K[Candidate Structures];
    K --> G;
    G --> L[Labelled Data (DFT)];
    L --> H;
    H --> M[MLIP Model v1];
    M --> I;
    I --> N{Uncertainty Check};
    N -- High Uncertainty --> O[Extract Structure];
    O --> G;

    subgraph Data & Models
        J
        K
        L
        M
    end

    L -- update --> P[Database (ASE DB)];
    J -- write --> P;
    M -- write --> P;

    N -- Low Uncertainty --> Q[Simulation Complete];

    style A fill:#f9f,stroke:#333,stroke-width:2px
    style C fill:#ccf,stroke:#333,stroke-width:2px
    style Q fill:#9f9,stroke:#333,stroke-width:2px
```

**Data Flow:**

1.  **Initialization**: The user provides a minimal `input.yaml`. The `Config Expander` heuristic engine parses this, determines the system type (alloy, molecular, etc.), and generates a complete `exec_config_dump.yaml` with optimized default parameters for the entire workflow.
2.  **Initial Seeding (Module A)**: The `Structure Generator` creates a diverse set of initial atomic configurations without DFT. It uses physically-motivated techniques based on the system type: SQS for alloys, Normal Mode Sampling for molecules, AIRSS for ionic materials, and Deep Rattling/Melt-Quench for covalent systems.
3.  **Exploration (Module B)**: The `Explorer & Sampler` takes the initial structures and uses a fast, pre-trained universal potential (the "surrogate model") to run large-scale MD simulations. It generates millions of frames, calculates structural descriptors (optimised with `Numba`), performs dimensionality reduction, and uses stratified sampling to select a small, diverse, and information-rich set of candidate structures for DFT labelling.
4.  **Labeling (Module C)**: The `Labeling Engine` receives the candidate structures and performs automated, robust DFT calculations using Quantum Espresso. It handles parameter automation (k-points, cutoffs via SSSP), error recovery, and ensures high-quality Force/Stress/Energy labels are generated. All results and metadata are stored in the central database.
5.  **Training (Module D)**: The `Training Engine` uses the DFT-labelled data to train the first version of the MLIP. It employs a "Delta Learning" approach, training the model to predict the residual between a physical baseline potential (e.g., ZBL) and the DFT values.
6.  **Simulation & Active Learning (Module E)**: The trained MLIP is passed to the `Simulation Engine` for dynamics simulations (MD or kMC). The engine performs "On-the-fly (OTF)" inference, continuously monitoring the model's uncertainty. If the uncertainty exceeds a predefined threshold, the simulation is paused, the problematic structure is extracted (using a periodic embedding approach to avoid surface effects), and it is sent back to the Labeling Engine (Step 4) to be added to the training set. The model is retrained, and the simulation resumes. This loop continues until the potential is robust over the desired simulation landscape.

## 4. Design Architecture

The project will be structured as a modern Python package, managed by `pyproject.toml` and the `uv` package manager. The source code will reside in a `src/mlip_autoprope` directory.

**File Structure:**

```
.
├── dev_documents/
├── src/
│   └── mlip_autoprope/
│       ├── __init__.py
│       ├── cli.py              # Main entry point (using Typer/Click)
│       ├── config/             # Pydantic models for config files
│       │   ├── minimal.py
│       │   └── full.py
│       ├── core/               # Core data structures (e.g., Atoms, Calculation)
│       ├── database/           # ASE DB interface and schema
│       ├── heuristics/         # The "Config Expander" logic
│       └── modules/
│           ├── __init__.py
│           ├── a_structure_generator.py
│           ├── b_explorer_sampler.py
│           ├── c_labeling_engine.py
│           ├── d_training_engine.py
│           └── e_simulation_engine.py
├── tests/
└── pyproject.toml
```

**Key Classes/Functions Overview:**

*   `cli.py`: Implements the main command-line interface, e.g., `uv run mlip-pipe run input.yaml`.
*   `config/`: Defines Pydantic models for strict validation and parsing of `input.yaml` and `exec_config_dump.yaml`.
*   `heuristics/ConfigExpander`: The core class responsible for the two-tier configuration strategy. It contains methods to analyze the minimal input and generate the full, executable configuration.
*   `modules/`: Each module will be implemented as a separate Python file containing a primary class (e.g., `StructureGenerator`) with a public `run()` method that accepts the system configuration and returns the generated artifacts.
*   `database/DatabaseInterface`: A wrapper around the ASE database to provide methods for querying, storing, and retrieving calculation data, ensuring transactional integrity.

**Data Models:**

The primary data objects will be based on the Atomic Simulation Environment (ASE) `Atoms` object, which provides a standard representation for atomic structures, including cell, positions, and momenta. All DFT calculations and MLIP training results will be linked to these objects within the database, providing a clear chain of provenance. Pydantic models will be used extensively for configuration and API validation to ensure data integrity throughout the pipeline.

## 5. Implementation Plan

The project is decomposed into five logical cycles, each building upon the last to deliver a functional, progressively more sophisticated system.

**CYCLE01: Core Engine and Workflow Foundation**
This cycle focuses on establishing the project's backbone. The primary goal is to create a workflow that can take a pre-existing set of atomic structures, run DFT calculations on them using Quantum Espresso, and train a basic MLIP model. This involves implementing the core components of the Labeling Engine (Module C) and the Training Engine (Module D). We will also establish the modern Python project structure with `pyproject.toml` and `uv`, and define the initial database schema for storing calculation results. The Two-Tier configuration strategy will be designed, but the heuristic engine will be minimal, mostly passing through user-defined values.

**CYCLE02: Automated Structure Generation and Configuration Expansion**
This cycle focuses on removing the need for pre-existing structures by implementing the Structure Generator (Module A). This module will be able to take a chemical formula and generate a diverse set of initial configurations based on the material type (alloy, molecular, etc.). The core of the `Config Expander` heuristic engine will also be developed in this cycle. It will automatically determine bond types, estimate necessary DFT parameters (cutoffs, k-points), and populate the `exec_config_dump.yaml` from a minimal `input.yaml`, realizing the core vision of user-friendliness.

**CYCLE03: High-Throughput Exploration and Optimisation**
This cycle addresses the computational bottleneck of data generation by implementing the Explorer & Sampler (Module B). The focus will be on integrating a pre-trained universal potential (e.g., MACE-MP) to perform fast surrogate-based MD simulations. A key part of this cycle is the development of performance-critical code for descriptor calculation and sampling, which will be heavily optimised using `Numba`. This will enable the system to efficiently scan vast configuration spaces and intelligently select only the most valuable data points for expensive DFT calculations.

**CYCLE04: Active Learning and Simulation Integration**
This cycle closes the loop by implementing the Simulation Engine (Module E) and the On-the-fly (OTF) active learning capability. This involves integrating LAMMPS as the primary MD/kMC engine and developing the real-time uncertainty monitoring system. When high uncertainty is detected, the system will be able to pause the simulation, extract the structure, send it for DFT labelling, retrain the model, and resume. This cycle will also focus on JIT-compiling the kMC event loop with `Numba` for high-speed exploration of rare events.

**CYCLE05: Finalisation, UI, and Documentation**
The final cycle will focus on hardening the system, improving error handling, and refining the user experience. This includes developing a more user-friendly command-line interface, potentially with progress bars and informative logging. Comprehensive documentation will be written, covering user tutorials, API references, and best practices. The entire system will undergo rigorous testing to ensure stability and reproducibility across different material systems.

## 6. Test Strategy

Testing will be a critical component of each cycle, ensuring the reliability and correctness of the complex scientific workflow. The strategy is multi-layered, combining unit, integration, and user acceptance testing.

**CYCLE01: Core Engine and Workflow Foundation**
*   **Unit Tests**: Mock the Quantum Espresso executable. Test the `LabelingEngine`'s ability to correctly generate QE input files based on various configurations. Test the `TrainingEngine`'s data handling and its interface with the underlying MLIP framework (mocked). Test the database interface for correct data insertion and retrieval.
*   **Integration Tests**: Run tests with a real (but small) Quantum Espresso calculation on a simple system (e.g., a silicon dimer). Verify that the entire flow from a structure input to a saved MLIP model in the database completes successfully and that the numerical outputs are within an expected tolerance.

**CYCLE02: Automated Structure Generation and Configuration Expansion**
*   **Unit Tests**: For each generator (SQS, NMS, etc.), provide input compositions and verify that the output structures are physically plausible (e.g., no overlapping atoms, correct stoichiometry). Test the `ConfigExpander` heuristic engine by providing various minimal `input.yaml` files and asserting that the generated `exec_config_dump.yaml` contains correct and complete parameters.
*   **Integration Tests**: Test the full workflow from a minimal `input.yaml` (e.g., "Si2") to the generation of initial structures and the successful execution of the (previously tested) DFT and training workflow on those structures.

**CYCLE03: High-Throughput Exploration and Optimisation**
*   **Unit Tests**: Test the descriptor calculation functions (optimised with `Numba`) against reference implementations for correctness. Test the sampling algorithms to ensure they select a diverse subset of structures from a known dataset. Mock the universal potential's inference call to test the data flow within the sampler.
*   **Integration Tests**: Run the sampler on a short, real surrogate MD trajectory. Verify that it correctly identifies and extracts distinct structures and passes them to the DFT engine. Performance benchmarks will be established to ensure the `Numba`-optimised code provides the expected speed-up.

**CYCLE04: Active Learning and Simulation Integration**
*   **Unit Tests**: Mock the LAMMPS executable. Test the uncertainty quantification logic by providing mock model outputs and verifying that the system correctly triggers the retraining loop. Test the periodic embedding extraction to ensure it correctly carves out cells without surface artifacts.
*   **Integration Tests**: Run a short, end-to-end active learning loop. Start with a poorly-trained model, run a small MD simulation, and verify that the system detects high uncertainty, pauses, retrains with a new DFT point, and resumes with a measurably improved model.

**CYCLE05: Finalisation, UI, and Documentation**
*   **Unit Tests**: Test the CLI argument parsing and user feedback mechanisms (e.g., error messages, progress indicators).
*   **User Acceptance Testing (UAT)**: Define a set of end-to-end test cases for different material types (e.g., an alloy, a molecule). Run the full pipeline for each case and verify that it produces a valid, reasonably accurate potential without crashing or requiring user intervention. The final potential's accuracy will be validated against known physical properties (e.g., lattice constant, bulk modulus).
