# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) system represents a paradigm shift in computational materials science, aiming to create a fully autonomous, end-to-end solution for the generation of high-fidelity, first-principles-accurate machine learning interatomic potentials (MLIPs). The foundational philosophy of this project is to "remove the human expert from the loop." This principle directly confronts the most persistent and critical bottleneck in the field: the laborious, time-consuming, and often intuition-driven process of curating training data for MLIPs. By seamlessly integrating and automating every stage of the workflow—from the initial generation of diverse atomic structures and their subsequent labelling with quantum mechanical data, through to iterative model refinement via active learning and deployment in long-timescale simulations—the system empowers researchers. It allows them to transcend the tedious mechanics of potential development and focus their expertise on high-level scientific inquiry. The ultimate vision is a system where a user can simply specify a chemical composition and target thermodynamic conditions, and the pipeline will autonomously deliver a production-quality MLIP, ready to unlock the secrets of complex material phenomena such as phase transformations, chemical reactions, diffusion, and defect dynamics.

The traditional approach to developing MLIPs is fraught with challenges that limit both its speed and accessibility. It relies profoundly on the deep domain knowledge and chemical intuition of an experienced researcher. This expert must decide which atomic configurations are most salient for training, a process that often involves running computationally prohibitive Ab Initio Molecular Dynamics (AIMD) simulations and then painstakingly hand-picking representative snapshots. This methodology is not only slow and incredibly expensive in terms of supercomputing resources, but it is also fundamentally subjective. The quality and transferability of the resulting potential are intrinsically tied to the biases and experience of the researcher who trained it, often leading to models that perform well within their narrow training domain but fail catastrophically when used to explore new phenomena. MLIP-AutoPipe is designed to systematically dismantle this reliance on human intervention. It achieves this by operationalising expert knowledge through a sophisticated combination of physics-based heuristics, cutting-edge machine learning techniques, and a robust active learning framework. This allows the system to intelligently, efficiently, and, most importantly, objectively explore the vast, high-dimensional potential energy surface of a given material.

The system's architecture is conceived as a modular, five-stage, orchestrated pipeline that mirrors a logical, scientific workflow. The process begins with the **Structure Generator (Module A)**, which employs physically-sound, non-DFT methods (such as Special Quasirandom Structures for alloys and Normal Mode Sampling for molecules) to algorithmically produce an initial dataset that is both diverse and physically relevant. This seed data then informs the **Explorer & Sampler (Module B)**, which leverages a pre-trained, general-purpose surrogate model (the MACE architecture) to rapidly and cheaply scan the material's phase space, identifying unique and informative configurations that warrant further investigation without resorting to costly DFT calculations. These carefully selected structures are then passed to the **Labelling Engine (Module C)**, a robust wrapper around the Quantum Espresso code. This engine performs automated, high-precision DFT calculations to generate the ground-truth energy, forces, and stresses that form the training data. Subsequently, these labels are fed into the **Training Engine (Module D)**, which employs a physically-motivated delta-learning approach to train a highly accurate and stable MLIP. The final component is the **Simulation Engine (Module E)**, which deploys the newly trained potential in large-scale simulations. This engine is the heart of the active learning loop; it actively monitors the model's uncertainty and, upon encountering a novel atomic environment, automatically triggers a retraining cycle. This entire, complex workflow is managed seamlessly by a central orchestrator, with all data, metadata, and models captured in a structured database to guarantee full scientific provenance and reproducibility.

## 2. System Design Objectives

The design of the MLIP-AutoPipe system is governed by a set of core objectives that, together, aim to produce a tool that is not only powerful but also autonomous, efficient, reliable, and accessible. The primary goal is to democratise the creation of high-quality MLIPs, making this capability available to a broader scientific community beyond elite computational groups.

*   **Autonomy and Usability**: The foremost objective is to create a system that operates with minimal human supervision. The ideal user interaction should be limited to specifying the chemical system and the scientific objective (e.g., "create a potential for FePt stable up to 1000K"). To achieve this, the system must automate a host of complex decisions traditionally made by experts. This includes selecting appropriate DFT convergence parameters (k-point meshes, plane-wave cutoffs), choosing pseudopotentials, defining a robust active learning strategy, and even tuning the hyperparameters of the MLIP itself. This automation is realised through the "Two-Tier Configuration" strategy, where a sophisticated heuristic engine translates a user's high-level intent into a detailed, explicit execution plan. This design principle is the cornerstone of making the technology accessible to experimentalists and students, not just computational experts.

*   **Computational Efficiency**: A relentless focus on minimising computational cost, particularly the use of the DFT engine, is a central design tenet. The pipeline is being engineered from the ground up to be "DFT-frugal." While DFT provides ground-truth accuracy, it is thousands of times slower than an MLIP. Therefore, its use must be judicious and targeted. The primary mechanism for achieving this efficiency is the intelligent, multi-stage data generation process. The use of a pre-trained universal surrogate model (MACE) in the exploration phase allows the system to sift through millions of potential configurations at a trivial cost, ensuring that the valuable DFT resources are exclusively allocated to the most informative atomic configurations—those that patch weaknesses in the model or expand its domain of applicability. This is a stark contrast to the brute-force AIMD approach, and it is expected to reduce the computational cost of generating a high-quality potential by at least an order of magnitude.

*   **Robustness and Reliability**: The system must be engineered for long, unsupervised runs on high-performance computing clusters. This demands a high degree of robustness and the ability to recover from common failures. The workflow is not a simple script; it is a complex, stateful process. The DFT Labelling Engine, for example, will not simply fail if a calculation does not converge. It will incorporate a multi-stage, automated error recovery protocol that attempts less aggressive mixing parameters or other common fixes, mimicking the actions of an experienced user. The active learning loop will feature a dynamic uncertainty threshold that adapts as the MLIP becomes more mature, preventing the system from getting stuck in a loop of retraining on minor variations of known structures. Furthermore, advanced treatments for periodic boundary conditions are mandated to ensure that data extracted from simulations is free from unphysical artefacts that could poison the training process.

*   **Physical Fidelity**: The pipeline must not be a mere black-box function approximator; it must be grounded in the principles of physics and chemistry. This objective influences the design of every module. The Structure Generator uses physically-motivated algorithms tailored to the bonding type of the material. More importantly, the Training Engine will implement a delta-learning approach. Instead of learning the absolute DFT energy, the MLIP will be trained to predict the *difference* between the DFT energy and the energy from a simple, classical, physics-based potential (like Lennard-Jones or ZBL). This technique forces the model to learn only the complex, quantum-mechanical component of the interactions, while the simple potential guarantees correct physical behaviour (e.g., strong repulsion at very short interatomic distances), leading to significantly more stable and accurate molecular dynamics simulations.

*   **Modularity and Extensibility**: The system is designed as a collection of five distinct, loosely-coupled modules, each with a well-defined responsibility. This modular architecture is a strategic choice to ensure the long-term viability and adaptability of the project. It allows individual components to be upgraded or replaced as new scientific methods and software tools emerge. For instance, the Labelling Engine could be extended with a VASP module to complement the initial Quantum Espresso one. The Training Engine could be adapted to incorporate novel MLIP architectures beyond ACE. This plug-and-play design facilitates parallel development, simplifies testing, and makes the system future-proof.

*   **Modern Software Engineering Practices**: The project will adhere to the highest standards of modern software development to ensure the codebase is maintainable, performant, and easy for new developers to contribute to. This includes using `uv` for fast, deterministic dependency management, defining the project structure and metadata declaratively in `pyproject.toml`, and leveraging high-performance libraries like Numba and JAX to accelerate computationally-intensive bottlenecks. All code will be statically typed, linted, and formatted, and the system will have a comprehensive, multi-layered test suite. This professional approach is critical for building a reliable scientific tool, not just a proof-of-concept script.

## 3. System Architecture

The MLIP-AutoPipe system is architected as a modular, state-driven workflow, orchestrated by a central control unit that manages the flow of data and execution between five specialized modules. This design ensures a clean separation of concerns, where the orchestrator handles the "what" and "when," while the individual modules handle the "how." All persistent data—atomic structures, DFT results, trained models, and simulation trajectories—is stored in a centralized ASE Database, which acts as the system's single source of truth, ensuring data provenance and enabling asynchronous operation.

```mermaid
graph TD
    A[User Input: input.yaml] --> B{Workflow Orchestrator};
    B --> C{Config Expander};
    C --> D[Full Config: exec_config_dump.yaml];

    subgraph "Cycle 1: Core Engine"
        E[Module C: Labelling Engine<br>(Quantum Espresso)]
        F[Module D: Training Engine<br>(ACE Delta Learning)]
    end

    subgraph "Cycle 2: Structure Generation"
        G[Module A: Structure Generator<br>(SQS, NMS, etc.)]
    end

    subgraph "Cycle 3: Exploration"
        H[Module B: Explorer & Sampler<br>(MACE Surrogate + DIRECT)]
    end

    subgraph "Cycle 4: Active Learning Simulation"
        I[Module E: Simulation Engine<br>(OTF MD/kMC)]
    end

    J[ASE Database<br>(SQLite)]

    B --> G;
    D --> E;
    D --> F;
    D --> H;
    D --> I;

    G -- Initial Structures --> J;
    H -- Candidate Structures --> J;
    J -- Structures to Label --> E;
    E -- Labelled Data (E, F, S) --> J;
    J -- Training Data --> F;
    F -- Trained MLIP --> I;
    I -- Simulation --> K{Uncertainty Quantification};
    K -- High Uncertainty Structure --> J;
    I -- Trajectory Data --> L[Analysis & Output];

    style J fill:#f9f,stroke:#333,stroke-width:2px
```

**Detailed Component Descriptions:**

1.  **Workflow Orchestrator**: This is the brain of the pipeline, containing the main business logic that dictates the sequence of operations. It is responsible for initiating the workflow, managing the state of the system, and directing the flow of data between the database and the various processing modules. For example, it decides when to transition from the exploration phase to the active learning phase and manages the iterative retraining loop. It is the active component that drives the entire process forward.

2.  **Config Expander**: This is the primary interface between the user and the system's complex internal machinery. It embodies the "Two-Tier Configuration" strategy. Its sole responsibility is to ingest the user's high-level, minimal `input.yaml` file and, using a rich, built-in library of physics-based heuristics, expand it into a comprehensive, explicit `exec_config_dump.yaml` file. This generated file contains every parameter needed for the run, from the precise Quantum Espresso convergence settings to the learning rate for the MLIP training. This module effectively encapsulates a significant amount of expert knowledge, making the system accessible to non-experts.

3.  **Module A: Structure Generator**: This module is responsible for creating the initial "seed" data for the learning process. It is the first step in replacing the human's intuition. The module intelligently analyses the user's input composition to determine the material's fundamental bonding character (e.g., metallic alloy, covalent network, molecule). It then dispatches to the most physically appropriate algorithmic structure generator—for instance, using the Special Quasirandom Structures (SQS) method to model disordered alloys or Normal Mode Sampling (NMS) to explore the vibrational space of a molecule. Its output is a diverse yet physically plausible set of initial structures that provide a solid foundation for the subsequent exploration phase.

4.  **Module B: Explorer & Sampler**: This module is the heart of the system's "DFT-frugal" design. Its purpose is to perform a vast and rapid exploration of the material's conformational space at a negligible computational cost. It achieves this by using a pre-trained, general-purpose surrogate model (MACE) to run large-scale molecular dynamics simulations. It then processes the resulting massive trajectories using the DIRECT sampling algorithm. This algorithm uses atomic environment descriptors and clustering to identify a small, diverse, and highly informative subset of structures that are most likely to improve the MLIP. This module acts as a powerful information filter, ensuring that only the most valuable structures are passed on for expensive DFT analysis.

5.  **Module C: Labelling Engine**: This module serves as a robust, automated interface to the "ground truth" engine: the Quantum Espresso DFT code. It takes an atomic structure as input, automatically generates the precise input file required by QE based on the parameters in the full configuration, executes the `pw.x` binary, and then parses the complex text output to extract the essential energy, forces, and stress tensor. Crucially, it also includes a sophisticated error-handling and recovery subsystem designed to manage common DFT convergence failures without crashing or requiring human intervention.

6.  **Module D: Training Engine**: This module is responsible for training the machine learning potential. It retrieves the labelled training data (structures and their corresponding DFT results) from the database. It then implements the delta-learning scheme, calculating the contribution from a simple baseline potential and subtracting it from the DFT ground truth. It then trains a high-performance ACE (Atomic Cluster Expansion) model on these residual values. This module also includes logic for automated hyperparameter optimisation, ensuring the final model is as accurate as possible.

7.  **Module E: Simulation Engine**: This is the module that closes the active learning loop. It takes the latest trained MLIP from the Training Engine and uses it to run production-style simulations (either Molecular Dynamics or Kinetic Monte Carlo). Its most critical function is the on-the-fly (OTF) uncertainty quantification. During the simulation, it constantly assesses the model's confidence in its own predictions. If it encounters an atomic configuration that it deems novel or uncertain, it pauses the simulation, carefully extracts the structure, and sends it back to the Labelling Engine for DFT calculation and subsequent retraining. This iterative refinement process is what allows the MLIP to grow its domain of validity and become progressively more robust.

8.  **ASE Database**: This is the central, passive repository for all data generated and consumed by the pipeline. Implemented using the flexible and powerful ASE (Atomic Simulation Environment) database interface (backed by SQLite), it stores every atomic structure, its current state in the workflow, its calculated DFT properties, and the metadata for every trained MLIP. This centralized, structured approach is essential for ensuring data provenance, allowing for full reproducibility of results, and enabling the system to be stopped and restarted without losing its state.

## 4. Design Architecture

The software will be developed as a modern, high-performance Python command-line application. The design prioritises maintainability, testability, and extensibility through a clear, modular structure and adherence to best practices in software engineering. The project will be packaged and distributed using standard Python tools, making it easy to install and use.

**Detailed File Structure:**

The code will be organised within the `src/` directory to follow standard Python packaging conventions.

```
src/
└── mlip_autopipec/
    ├── __init__.py
    ├── main.py             # CLI entry point using the Click library. Defines commands and parses arguments.
    ├── orchestrator.py     # Contains the main Orchestrator class, holding the high-level business logic.
    ├── config/
    │   ├── __init__.py
    │   ├── expander.py     # Contains the ConfigExpander class and its heuristic logic.
    │   └── models.py       # Defines the Pydantic models for MinimalConfig and FullConfig.
    ├── data/
    │   ├── __init__.py
    │   ├── database.py     # Contains the AseDB wrapper class for all database interactions.
    │   └── models.py       # Defines Pydantic models for internal data structures like DFTResult.
    ├── modules/
    │   ├── __init__.py
    │   ├── structure_generator.py # Implements Module A.
    │   ├── explorer.py            # Implements Module B.
    │   ├── labelling_engine.py    # Implements Module C.
    │   ├── training_engine.py     # Implements Module D.
    │   └── simulation_engine.py   # Implements Module E.
    └── utils/
        ├── __init__.py
        ├── dft_utils.py       # Low-level helper functions for parsing Quantum Espresso files.
        ├── parallel.py      # Utilities for managing parallel execution of tasks.
        └── descriptors.py   # Numba-optimised functions for descriptor calculations (e.g., SOAP).
pyproject.toml              # Project definition, dependencies, and entry points.
input.yaml                  # Example user-facing minimal configuration file.
```

**Key Class and Data Model Definitions:**

*   `main.py`: This file will be kept lean and will only contain the user interface logic. It will use the `click` library to define a clean command-line interface, for instance, `mlip-pipe run --config input.yaml`. Its role is to parse user input and delegate the actual work to the `Orchestrator`.

*   `Orchestrator`: This is the core stateful class of the application. It will be initialised with the configuration and will have high-level methods that correspond to the major phases of the workflow, such as `run_pipeline()`, `run_initial_seeding()`, `run_exploration_phase()`, and `run_active_learning_loop()`. It holds the logic for transitioning between these states.

*   `ConfigExpander`: A largely stateless class whose primary method, `expand(minimal_config)`, will be a pure function that deterministically transforms the minimal user input into the full, explicit execution configuration. This makes its logic highly testable.

*   `AseDB`: This class will provide a high-level API for database operations, abstracting away the underlying `ase.db` syntax. It will have methods like `add_structure_for_labelling(atoms)`, `get_completed_training_pairs()`, and `save_model_metadata(model)`, which will make the code in the modules cleaner and more readable. The database schema will include tables for `structures` (with columns for `id`, `state`, `ase_atoms_object`), `dft_results` (linked to `structures`), and `mlip_models`.

*   `LabellingEngine`: An object-oriented wrapper around the `pw.x` command-line tool. An instance of this class will manage the entire lifecycle of a DFT calculation: generating the input from an `ase.Atoms` object, executing the subprocess, monitoring for errors, and parsing the output.

*   `TrainingEngine`: This class will abstract the details of the specific MLIP framework being used (e.g., `mace-torch` or `pacemaker`). Its `train(dataset)` method will contain the logic for the delta-learning calculation and the boilerplate code required to run the underlying training library.

**Data Models (Pydantic):**

The extensive use of Pydantic for both configuration and internal data models is a key design choice. It provides several benefits:
1.  **Validation**: Pydantic models automatically validate incoming data, ensuring that, for example, the user's `input.yaml` is correctly formatted and that all necessary fields are present. This prevents runtime errors deep within the application logic.
2.  **Type Safety**: They provide static type hints, which improves code readability and allows for static analysis, catching potential bugs before the code is even run.
3.  **Self-Documentation**: The Pydantic model definitions themselves serve as a clear and explicit schema for the data structures used throughout the application.

*   `MinimalConfig`: Will define the simple, user-facing `input.yaml` structure with just a few required keys.
*   `FullConfig`: Will be a large, nested Pydantic model defining the entire execution plan, with sub-models for `DFTConfig`, `MDConfig`, `MLIPConfig`, etc.
*   `DFTResult`: A structured model for passing the results of a DFT calculation between the `LabellingEngine` and the `AseDB`, ensuring no data is accidentally lost or misinterpreted.

This structured, modular, and type-safe design is essential for building a complex scientific workflow application that is both reliable and maintainable.

## 5. Implementation Plan

The project's implementation is strategically decomposed into five logical, sequential development cycles. Each cycle delivers a functional, self-contained increment of the total system, allowing for iterative development, testing, and refinement. This approach ensures that a working, albeit partial, system exists at the end of each cycle, providing tangible progress and a solid foundation for subsequent work.

**CYCLE 01: The Core Engine** (Words: 520)
This inaugural cycle is foundational, focused on establishing the absolute core functionality of the pipeline: the ability to perform a DFT calculation for a given atomic structure and subsequently train a basic MLIP from the resulting data. This cycle will deliberately ignore user-facing features, automation, and advanced learning loops. The singular goal is to build and validate the essential data processing backbone. The primary deliverables will be `Module C (Labelling Engine)` and `Module D (Training Engine)`. The Labelling Engine will be implemented as a robust Python wrapper for Quantum Espresso. Its key feature will be the ability to take an ASE Atoms object, automatically generate a valid QE input file, execute the `pw.x` subprocess, and parse the text output to extract energy, forces, and stress. A crucial part of this implementation will be the development of a multi-stage error-recovery protocol to handle common SCF convergence failures, making the module resilient. Concurrently, the Training Engine will be developed to interface with the ACE training library (e.g., `pacemaker` or `mace-torch`). It will be responsible for taking a dataset of labelled structures and producing a serialised, trained MLIP file. The critical delta-learning mechanism, which is key to the physical accuracy of the MLIP, will be implemented within this module. To manage the data flow between these engines, a basic `AseDB` wrapper class will be created to handle the storage and retrieval of structures and their corresponding DFT results in an SQLite database. The initial project structure, including the `pyproject.toml` file and the `uv` virtual environment setup, will also be established. The cycle will culminate in a simple orchestrator script that executes a hardcoded, linear workflow: take one predefined structure, label it, and train a model from this single data point. This will provide a crucial end-to-end validation of the core concept.

**CYCLE 02: Initial Structure Generation & Configuration** (Words: 515)
The second cycle builds directly upon the core engine, focusing on making the system usable and automated by removing the need for hardcoded inputs. This is achieved by introducing the user-friendly configuration system and the automated initial structure generation. The centerpiece of this cycle is the development of the "Two-Tier Configuration" strategy, embodied in the `Config Expander` component. This heuristic engine will be designed to read a minimalistic, human-friendly `input.yaml` file and expand it into the comprehensive `exec_config_dump.yaml` required for execution. This involves creating the Pydantic models for both configuration files and implementing the heuristic logic to make intelligent, physics-based decisions (e.g., selecting SSSP pseudopotentials and cutoff energies based on the elemental composition). The second major deliverable is `Module A (Structure Generator)`. This module will implement the logic for analysing the input composition from the configuration file to automatically determine the material's bonding type (alloy, molecular, ionic, or covalent). Based on this classification, it will then dispatch to the appropriate structure generation algorithm, such as SQS for alloys via the `icet` library, or Normal Mode Sampling for molecules. The output of this module—a diverse list of physically plausible seed structures—will be saved to the database. The `Orchestrator` will be significantly updated to integrate this new initialisation phase. The main workflow will now start by reading the user's config file, calling the Config Expander, and then using the full configuration to drive the Structure Generator. By the end of this cycle, the pipeline will be capable of initiating a full data generation workflow from a single user command and a simple input file, marking a major step towards true automation.

**CYCLE 03: Efficient Exploration & Optimisation** (Words: 530)
This cycle tackles the critical challenge of computational efficiency, aiming to make the pipeline "DFT-frugal" through the implementation of an intelligent exploration phase. This is a performance-sensitive cycle that introduces `Module B (Explorer & Sampler)`. The core strategy is to minimise the number of expensive DFT calculations by first performing a vast, cheap exploration using a pre-trained universal surrogate model, MACE-MP-0. The first task is to integrate the `mace-torch` library and implement the functionality to run large-scale molecular dynamics simulations using the MACE model as the ASE calculator. These simulations will generate massive trajectory files containing millions of candidate configurations. The second, more complex task is to develop the DIRECT sampling algorithm to intelligently select a small, diverse, and informative subset of these configurations for labelling. This involves implementing or integrating a high-performance atomic environment descriptor, specifically the Smooth Overlap of Atomic Positions (SOAP). A key part of this cycle is the optimisation of this descriptor calculation. A pure Python implementation would be too slow, so a significant effort will be dedicated to writing a custom, Numba-accelerated (`@jit`) version of the SOAP kernel that can process the trajectories efficiently. Once descriptors are calculated, an efficient clustering algorithm will be implemented to group them, followed by a stratified sampling logic to select representative frames from each cluster, ensuring both common and rare configurations are captured. The `Orchestrator`'s workflow will be updated to insert this exploration phase between the initial structure generation and the labelling phase. The end result of this cycle will be a smart, efficient system that uses its computational resources wisely, filtering a vast amount of low-value information to focus on a small kernel of high-value data points.

**CYCLE 04: The Active Learning Loop** (Words: 540)
This cycle implements the final, and most crucial, piece of the core technology: the closed-loop, on-the-fly (OTF) active learning capability. This transforms the pipeline from a static, feed-forward process into a dynamic, self-correcting system. The main deliverable is `Module E (Simulation Engine)`. This module will be responsible for taking a trained MLIP from Module D and using it to run production-style molecular dynamics simulations. The key innovation is the integration of an uncertainty quantification (UQ) mechanism directly into the simulation loop. For each MD step, the engine will calculate not only the forces but also a metric for the MLIP's uncertainty in that prediction (e.g., the variance between members of a committee of models). This uncertainty value will be compared against a dynamic threshold. A significant part of the work in this cycle is developing this dynamic threshold mechanism; it will be statistically determined from the distribution of uncertainties on the *existing* training set, allowing it to adapt as the model becomes more confident. If the threshold is exceeded, the simulation is automatically paused. At this point, another critical piece of technology is invoked: the periodic structure extraction algorithm. This is a non-trivial function that carefully cuts the high-uncertainty configuration out of the bulk simulation cell, creating a finite, charge-neutral cluster with a buffer region to avoid unphysical boundary artefacts. This extracted structure is then sent back to the `AseDB` to be picked up by the `Labelling Engine`. The `Orchestrator` logic will be substantially refactored to manage this cyclical workflow: train -> simulate -> detect uncertainty -> extract -> label -> retrain -> resume simulation. Upon completion, the MLIP-AutoPipe will be a fully autonomous learning machine, capable of iteratively and systematically improving its own robustness.

**CYCLE 05: Advanced Simulations & User Interface** (Words: 525)
The final development cycle focuses on polishing the system into a complete, professional, and powerful software package. This involves work on two parallel fronts: enhancing the scientific capabilities and refining the user experience. The first major feature is the extension of `Module E (Simulation Engine)` to support advanced, long-timescale simulation methods, specifically Adaptive Kinetic Monte Carlo (kMC). This is essential for studying rare events like atomic diffusion, which are inaccessible to standard MD. The implementation will require developing or integrating algorithms for saddle-point searching (e.g., the Dimer method) using the MLIP as the energy landscape. To make this computationally feasible, a "Tiered Rate Calculation" strategy will be implemented, where a cheap, approximate method is used to screen for likely events before a more expensive Hessian calculation is performed to get accurate vibrational prefactors for only the most promising candidates. The second major front is the development of a polished, robust, and user-friendly Command-Line Interface (CLI) in `main.py` using the `click` library. This will replace any temporary scripts with a professional, self-documenting interface, providing commands like `mlip-pipe run` and `mlip-pipe run-kmc`. The `pyproject.toml` file will be configured to install this CLI as a system command. Finally, this cycle includes the crucial "last mile" tasks of software development: writing comprehensive user documentation, including a getting-started guide and tutorials, and packaging the entire project for distribution via PyPI. This cycle transforms the powerful backend developed in the previous cycles into a finished, accessible, and distributable product for the scientific community.

## 6. Test Strategy

A comprehensive, multi-layered testing strategy is essential to ensure the correctness, robustness, and physical accuracy of the MLIP-AutoPipe system. Testing will be a continuous process throughout all development cycles, with a strong emphasis on automation. The strategy is divided into distinct levels: unit tests for isolated component verification, integration tests for validating interactions between components, and end-to-end tests for verifying the entire workflow.

**CYCLE 01: The Core Engine** (Words: 530)
Testing for the foundational cycle is critical, as all subsequent cycles will build upon this core. The focus is on verifying the data processing backbone.
*   **Unit Tests**:
    *   `Labelling Engine`: This component will be tested extensively by mocking the `subprocess.run` call to the external `pw.x` executable. The tests will not run a real DFT calculation. Instead, they will verify that the `_generate_input_file` method produces a character-for-character correct Quantum Espresso input file for a variety of `ase.Atoms` objects (e.g., bulk crystals, molecules, systems with different symmetries). We will use pre-validated "golden" input files as the ground truth for these comparisons. The `_parse_output` method will be tested with a suite of sample QE output files, including those from successful runs, runs with various error messages, and runs from different QE versions, ensuring the parsing logic is robust. The error recovery logic will be tested by configuring the mocked subprocess to initially fail and then asserting that the engine attempts to take corrective action.
    *   `Training Engine`: The external ACE training library will be completely mocked. The primary goal is to test our internal logic, not the library itself. Unit tests will focus on the `_compute_delta` method, providing it with a known `DFTResult` and a mock baseline calculator, and then asserting that the calculated residual values for energy, forces, and stress are correct. We will also test that the data is formatted and passed to the mocked training library's API precisely as the library expects.
    *   `AseDB`: The database wrapper will be tested against a temporary, in-memory SQLite database (`":memory:"`). The tests will verify all CRUD (Create, Read, Update, Delete) operations. We will assert that adding an atom returns a correct ID, that writing a `DFTResult` correctly updates the status of that atom to 'labelled', and that methods like `get_atoms_to_label` only return atoms with the appropriate status.

**CYCLE 02: Initial Structure Generation & Configuration** (Words: 520)
Testing for this cycle focuses on the user-facing entry point and the automated data seeding.
*   **Unit Tests**:
    *   `Config Expander`: As this should be a pure function, it is highly suitable for unit testing. A wide range of `MinimalConfig` inputs will be created as test cases. For each case, we will run the `expand` method and perform deep assertions on the resulting `FullConfig` object. For example, for an input of `{"elements": ["H", "O"]}`, we will assert that `structure_type` is correctly identified as 'Molecule' and that the DFT settings are appropriate for those elements. We will test edge cases, like unsupported elements, to ensure the expander raises clear, informative errors.
    *   `Structure Generator`: The external libraries used for structure generation (e.g., `icet` for SQS) will be mocked. The tests will verify that our code correctly constructs the input for these libraries based on the `FullConfig` and correctly parses their output back into a list of `ase.Atoms` objects. For each generation method, we will provide a simple input and assert that the output is a list of structurally diverse and physically plausible configurations with the correct atom types and counts.
*   **Integration Tests**: A key integration test will cover the flow from the CLI input to the database. It will use a `click.testing.CliRunner` to invoke the main command with a path to a temporary `input.yaml`. It will then assert that a valid `exec_config_dump.yaml` was created and that the temporary database has been populated with new structures in the 'initial_generated' state.

**CYCLE 03: Efficient Exploration & Optimisation** (Words: 540)
Testing for this cycle must cover both functional correctness and performance.
*   **Unit Tests**:
    *   `Explorer & Sampler`: The core logic of the DIRECT sampler will be tested by providing it with a small, pre-computed toy trajectory. This trajectory will be designed to have distinct structural motifs (e.g., half solid-like, half liquid-like). We will then assert that the clustering and sampling logic correctly identifies these groups and selects a balanced number of representative samples from each, proving the algorithm's correctness. The MD simulation itself will be mocked.
    *   `utils.descriptors`: This is a performance-critical area requiring rigorous testing. We will maintain two implementations of the SOAP descriptor calculation: a simple, readable, pure-Python version and the highly optimised Numba version. The unit tests will run both versions on a wide range of atomic structures and use `np.testing.assert_allclose` to ensure their outputs are numerically identical to within a tight tolerance. This guarantees that the optimisation does not introduce bugs.
*   **Property-Based Tests**: We will use a library like `hypothesis` to generate thousands of random but physically valid atomic configurations. These will be fed into the Numba-optimised descriptor function to ensure it is robust, does not crash on unexpected inputs, and remains consistent with the reference implementation.
*   **Performance Tests**: A separate benchmarking suite, not run with the standard unit tests, will be created. It will measure the wall-clock time taken by the pure-Python and Numba versions of the descriptor calculation on a large, realistic trajectory. The test will assert that the speed-up is substantial (e.g., >50x), providing concrete validation that the cycle's performance goals have been met.

**CYCLE 04: The Active Learning Loop** (Words: 550)
Testing the dynamic, stateful active learning loop is the most complex challenge.
*   **Unit Tests**:
    *   `Simulation Engine`: The MLIP model will be mocked. To test the `run_otf_md` method, the mock MLIP will be programmed to return low uncertainty for N steps, and then a high uncertainty on step N+1. The test will assert that the MD loop runs for exactly N+1 steps and then correctly returns the final `Atoms` object. The `_extract_structure_for_labelling` method will be tested with various periodic structures, asserting that the resulting non-periodic cluster is correctly formed with the appropriate buffer region and PBC flags disabled. The `_update_dynamic_threshold` method will be tested by providing a list of mock training data with known uncertainties and asserting that the calculated percentile is correct.
    *   `Orchestrator`: The orchestrator's loop logic will be tested by mocking all the engine modules. The test will configure the mock `SimulationEngine` to find an uncertain structure on its first run, and finish successfully on its second. We will then assert that the `LabellingEngine` and `TrainingEngine` were called exactly once in between, proving that the retraining cycle was correctly triggered and the loop terminated as expected.
*   **Integration Tests**: A full integration test is too slow. Instead, a "toy problem" integration test will be built using a 2D analytical potential. The `LabellingEngine` will be mocked to return the true value of this potential. The system will then train a simple 2D ML model (like a Gaussian Process) and run the active learning loop. We can then visualise the sequence of points chosen for labelling, asserting that the system intelligently adds new training points in the regions where the MLIP's error is highest, thus validating the behaviour of the entire closed-loop system.

**CYCLE 05: Advanced Simulations & User Interface** (Words: 530)
The final cycle's testing focuses on the new scientific features and the user-facing interface.
*   **Unit Tests**:
    *   `Simulation Engine (kMC)`: The complex kMC logic will be tested against a simple, known energy landscape, likely a 2D analytical potential with well-defined minima and saddle points. We will place the "atom" in a minimum and test the `_find_saddle_points` method, asserting that it correctly identifies the location of the known transition state. The `_calculate_hessian` method will be tested by comparing its numerical finite-difference result with the analytical second derivatives of the toy potential. This validates the core algorithms in a controlled environment.
    *   `CLI`: The entire command-line interface will be tested using the `click.testing.CliRunner`. The `Orchestrator` and all underlying modules will be mocked. We will invoke the CLI with various commands and arguments, like `runner.invoke(cli, ['run', 'fake.yaml', '--max-iterations', '3'])`. The tests will assert that the correct methods on the mock `Orchestrator` were called with the correctly parsed arguments. We will also test for user errors, such as typos in commands or paths to non-existent files, and assert that the CLI provides the expected friendly and informative error messages.
*   **End-to-End (E2E) Tests**: Two final E2E tests will be created to validate the entire system.
    1.  The first test will run the main `mlip-pipe run` command on a very simple, real system (e.g., H2 or Si) with the real DFT code and real training, but configured for only one or two active learning iterations. The test will simply assert that the process completes with an exit code of 0 and produces a model file, proving that all components are integrated correctly.
    2.  The second test will use the model generated by the first test to run the `mlip-pipe run-kmc` command for a small number of steps, asserting that it completes successfully and produces a trajectory where at least one kMC event has occurred.
