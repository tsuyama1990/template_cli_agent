# SYSTEM_ARCHITECTURE.MD

## 1. Summary

The Machine Learning Interatomic Potential Automatic Pipeline (MLIP-AutoPipe) is an ambitious system designed to fundamentally automate the process of creating high-fidelity, first-principles-accurate interatomic potentials. The core philosophy driving this project is the explicit goal of "removing the human expert from the loop." For decades, the development of force fields and potentials has been a significant bottleneck in computational materials science, relying heavily on the intuition, experience, and manual effort of seasoned researchers. This traditional, artisanal process involves meticulously selecting appropriate training structures from a vast chemical space, ensuring the robust convergence of demanding Density Functional Theory (DFT) calculations, tuning a complex array of hyperparameters for the machine learning model, and rigorously validating the final potential. This process is not only time-consuming and computationally expensive but also inherently subjective and difficult to reproduce. MLIP-AutoPipe aims to replace this craft with a robust, intelligent, and fully autonomous computational workflow, transforming potential generation from a research project into a reliable, on-demand computational service.

By leveraging a synergistic combination of modern machine learning techniques—particularly active learning—and physically-grounded heuristics, the system will autonomously navigate the vast potential energy surface to generate optimal and computationally efficient training datasets. The user interaction is designed to be minimal and high-level; a domain scientist will only need to provide a chemical composition (e.g., "FePt") and perhaps a target operational range (e.g., temperature). From this simple input, the pipeline will orchestrate the entire end-to-end process. This includes generating an initial set of diverse and physically-sound structures using established algorithms, performing thousands of automated DFT calculations for high-fidelity labelling, training sophisticated and accurate MLIP models like the Atomic Cluster Expansion (ACE), and iteratively validating and refining the potential through on-the-fly (OTF) simulations that actively seek out the model's areas of weakness. The system is designed to be a "fire-and-forget" tool that empowers a broader range of scientists, who may not be experts in quantum chemistry or machine learning, to generate bespoke, high-quality potentials for their specific materials of interest. This democratization of potential development will unlock new research possibilities and accelerate discoveries in complex phenomena such as phase transitions, defect dynamics, catalysis, and battery degradation, all of which require simulations of a scale, length, and accuracy previously unattainable without significant, prolonged manual intervention. The architecture is deliberately modular, allowing for future extensions with new structure generation algorithms, different DFT codes, or next-generation MLIP models, thereby ensuring its long-term relevance and adaptability in this rapidly evolving scientific field.

## 2. System Design Objectives

The design of MLIP-AutoPipe is guided by a set of clear, stringent objectives aimed at creating a system that is not only powerful and accurate but also reliable, efficient, and user-friendly. These objectives form the foundational principles upon which all architectural decisions are made.

*   **Autonomy and Usability**: The paramount objective is to create a system that is, for all practical purposes, fully autonomous. The ideal user experience involves a scientist specifying a material and, after some computation time, receiving a thoroughly validated potential ready for use in production simulations. This necessitates the creation of a sophisticated "Config Expander" module. This module acts as an expert system, translating high-level user intent (e.g., `elements: [Fe, Pt]`) into a detailed, low-level execution plan. It will automatically determine hundreds of parameters—such as DFT convergence settings (k-points, energy cutoffs), simulation protocols (temperature ramps, pressure controls), and MLIP model hyperparameters—based on built-in physical heuristics and best-practice rules encoded from materials science literature. For instance, it will query established databases or protocols like the SSSP to select appropriate pseudopotentials and cutoffs for the given elements. This radical simplification of the user interface is critical to the project's goal of removing the human expert from the operational loop, ensuring that the system is accessible to the widest possible scientific audience.

*   **Accuracy and Robustness**: The ultimate utility of the system depends on the quality of the potentials it produces. They must be accurate enough to faithfully reproduce first-principles (DFT) energies, forces, and stresses across a wide range of atomic configurations. This objective is met through a multi-pronged strategy. First, the system will use standardized, high-precision DFT protocols to ensure the quality of the training labels. Second, it will implement a "Delta Learning" approach, where the MLIP is trained to predict the *difference* between DFT results and a simpler, computationally cheap baseline physical potential (like Lennard-Jones). This allows the MLIP to focus its expressive power on capturing the complex, short-range quantum mechanical interactions that the baseline model misses. Robustness is equally critical for a long-running autonomous system. The pipeline must gracefully handle the inevitable failures of complex scientific codes. The DFT Labeling Engine will therefore incorporate a sophisticated, multi-stage error recovery mechanism, automatically attempting to fix convergence failures by adjusting parameters like electronic mixing schemes before flagging a structure as problematic.

*   **Computational Efficiency**: While accuracy is crucial, the system must be computationally efficient to be a practical tool rather than a theoretical curiosity. The brute-force approach of generating training data via extensive Ab Initio Molecular Dynamics (AIMD) is explicitly avoided due to its prohibitive cost. Instead, the system employs an intelligent, multi-tiered data generation strategy centered on active learning. A pre-trained, universal surrogate model (specifically, the MACE architecture) is used for rapid, large-scale exploration of the potential energy surface. This allows the system to sample millions of configurations and identify regions of high structural diversity or high model uncertainty at a tiny fraction of the cost of a DFT calculation. The subsequent active learning loop ensures that expensive DFT resources are allocated surgically, targeting only those atomic configurations that provide the most new information to the model. This information-theoretic approach maximizes the "knowledge gain" per DFT calculation, drastically reducing the number of labels required to achieve a target accuracy compared to random or grid-based sampling.

*   **Modularity and Extensibility**: The fields of materials informatics and machine learning are exceptionally dynamic. New algorithms, models, and simulation codes emerge on a continuous basis. To avoid obsolescence, the system is designed from the ground up with a highly modular architecture. Each core component (Structure Generator, Labeling Engine, Training Engine, etc.) is a distinct, swappable module with a well-defined and stable interface. This architectural principle ensures that the system can be easily extended and upgraded in the future. For example, if a user wishes to use a different DFT code like VASP, a new `VaspLabelingEngine` can be implemented that adheres to the same interface as the default `QuantumEspressoLabelingEngine`, and the orchestrator can use it without any other changes. Similarly, alternative MLIP frameworks like NequIP or Allegro, or novel active learning strategies, can be integrated by creating new plugin-like modules, ensuring the long-term viability and flexibility of the pipeline.

*   **Provenance and Reproducibility**: In science, results that are not reproducible are of little value. Every computational step, from the initial structure generation to the final simulation, must be meticulously tracked to ensure full scientific provenance. The system will use a centralized database (an ASE Database) to store not just the raw data (atomic structures, energies, forces), but also the complete metadata of its creation. This includes the precise configuration parameters used for every calculation, the versions of all external software (DFT code, training library), and even the version control hash of the pipeline's own source code at the time of the run. This comprehensive record-keeping ensures that any potential generated by the system is fully reproducible, providing a clear and unimpeachable audit trail for scientific validation, publication, and collaboration.

## 3. System Architecture

The MLIP-AutoPipe is architected as a modular, stateful pipeline orchestrated by a central workflow manager. This design promotes a clear separation of concerns, enhances testability, and ensures the system's extensibility. The components are designed to be loosely coupled, communicating primarily through a shared, persistent database and a unified configuration object, rather than through direct inter-module calls. This asynchronous-style communication, mediated by the database, makes the system more resilient to the failure of individual components.

```mermaid
graph TD
    subgraph User Interaction
        A[input.yaml] --> B{Config Expander};
    end

    subgraph Core Pipeline
        B --> C[exec_config_dump.yaml];
        C --> D{Workflow Orchestrator};
    end

    subgraph Data Flow & Persistence
        E[ASE Database (SQLite)]
    end

    subgraph Pluggable Modules
        M_A[Module A: Structure Generator]
        M_B[Module B: Explorer & Sampler]
        M_C[Module C: Labeling Engine]
        M_D[Module D: Training Engine]
        M_E[Module E: Simulation Engine]
    end

    D -- Manages & Schedules --> M_A;
    D -- Manages & Schedules --> M_B;
    D -- Manages & Schedules --> M_C;
    D -- Manages & Schedules --> M_D;
    D -- Manages & Schedules --> M_E;

    M_A -- "Writes Initial Structures" --> E;
    M_B -- "Writes Candidate Structures" --> E;
    M_C -- "Reads Unlabeled, Writes Labeled" --> E;
    M_D -- "Reads Labeled Data, Writes MLIP Artifact" --> E;
    M_E -- "Uses MLIP, Discovers & Writes Uncertain Structures" --> M_B;
    M_B -- "Selects & Queues Uncertain Structures" --> M_C;

    style User Interaction fill:#f9f,stroke:#333,stroke-width:2px
    style Core Pipeline fill:#ccf,stroke:#333,stroke-width:2px
    style Modules fill:#cfc,stroke:#333,stroke-width:2px
```

**Component Deep Dive:**

*   **Config Expander**: This crucial module serves as the primary user-facing component and the system's "expert translator." It ingests the user's high-level, minimal `input.yaml` and applies a sophisticated layer of physical heuristics and rule-based logic to generate the comprehensive `exec_config_dump.yaml`. This exhaustive file explicitly defines every parameter for the entire workflow, from DFT k-point meshes to MD simulation timesteps, effectively transforming a simple user request into a detailed, reproducible scientific protocol.

*   **Workflow Orchestrator**: This is the "brain" of the entire system. It operates as a state machine, reading the full configuration and directing the execution flow. It is responsible for managing the complex logic of the active learning loop: invoking the structure generator, dispatching batches of structures to the labeling engine, monitoring for completion, triggering the training engine, and initiating the simulation engine. It decides when the potential has converged based on pre-defined criteria (e.g., number of loops, or when the rate of finding new, uncertain structures drops below a threshold) and gracefully shuts down the pipeline.

*   **ASE Database**: The database is the central nervous system and shared memory of the application. It holds all atomic structures, their current state (e.g., `unlabeled`, `labeling_in_progress`, `labeled`, `failed`), their calculated DFT labels (energy, forces, stress), the trained MLIP models for each generation, and the full provenance metadata. By using the database as the primary communication channel, the modules are effectively decoupled, allowing them to operate independently and making the overall system more robust and easier to debug.

*   **Module A (Structure Generator)**: This module is responsible for creating the initial, physically diverse set of structures to "seed" the active learning process. It avoids expensive AIMD by employing a variety of established, cost-effective techniques tailored to different material types: Special Quasirandom Structures (SQS) for disordered alloys, Normal Mode Sampling (NMS) for molecular systems, and various crystal distortion methods (e.g., "rattling" and volume deformation) to capture elastic properties.

*   **Module B (Explorer & Sampler)**: This module is the heart of the system's computational efficiency. It leverages a fast, general-purpose surrogate model (MACE) to perform massive-scale molecular dynamics simulations, exploring a vast swath of the potential energy surface at a trivial cost. It then analyzes this huge trajectory, using geometric descriptors like SOAP to cluster the structures. From these clusters, it employs stratified sampling to select a small, highly informative subset of configurations to be passed to the expensive DFT labeling engine. It plays a second role during active learning, where it helps prioritize uncertain structures discovered by Module E.

*   **Module C (Labeling Engine)**: This module acts as a robust, automated wrapper around the chosen DFT code (Quantum Espresso). It takes a structure from the database, automatically generates the correct input file based on the system configuration, runs the DFT calculation, handles a wide range of common runtime errors through its recovery logic, and parses the output file, storing the resulting labels and metadata back in the database.

*   **Module D (Training Engine)**: This module implements the Delta Learning strategy. It queries the database for all available labeled data, calculates the contribution from a simple baseline potential, and then trains the high-fidelity ACE potential on the residual (the difference between DFT and the baseline). It saves the trained model artifact back into the database, ready for use by the simulation engine.

*   **Module E (Simulation Engine)**: This is the final consumer of the trained potential and the component that "closes" the active learning loop. It runs large-scale, long-time simulations (MD/kMC) to validate the potential and explore physical phenomena. During these simulations, it continuously monitors the MLIP model's uncertainty. If the simulation enters a configuration where the model is "unsure" of its prediction (i.e., the uncertainty exceeds a dynamic threshold), the engine pauses, extracts that structure, and sends it back to the pipeline to be labeled and incorporated into the next generation of the training set.

## 4. Design Architecture

The project will be structured as a standard, modern Python package, installable via `uv` and adhering to best practices for code organization and dependency management. The design emphasizes a strict separation of concerns between data structures (the "what"), business logic (the "how"), and infrastructure interaction (the "where"). This layered architecture is essential for creating a system that is both maintainable and highly testable.

**File Structure:**
```
.
├── pyproject.toml
├── README.md
├── src
│   └── mlip_autopipec
│       ├── __init__.py
│       ├── cli.py                # User-facing command-line interface (Typer/Click)
│       ├── config.py             # Pydantic models for configuration (Minimal and Full)
│       ├── orchestrator.py       # The main Workflow Orchestrator and state machine logic
│       ├── database.py           # High-level wrapper for ASE DB interaction
│       └── modules
│           ├── __init__.py
│           ├── a_structure_generator.py
│           ├── b_explorer_sampler.py
│           ├── c_labeling_engine.py
│           ├── d_training_engine.py
│           └── e_simulation_engine.py
└── tests
    ├── test_config.py
    ├── test_orchestrator.py
    └── modules
        ├── test_a_structure_generator.py
        ├── test_b_explorer_sampler.py
        ├── test_c_labeling_engine.py
        ├── test_d_training_engine.py
        └── e_simulation_engine.py
```

**Schema-First Design with Pydantic:**
The system's entire configuration will be rigidly defined by Pydantic models in `config.py`. This "Schema-First" approach is a cornerstone of the design. It ensures type safety, provides automatic and transparent validation, and serves as a form of executable documentation for all parameters. This prevents a wide class of common runtime errors caused by misconfigured or missing parameters.

*   **`MinimalConfig`**: This model represents the user's public-facing `input.yaml`. It will be very simple, containing only a few required fields like `elements` and `composition`, with sensible defaults for everything else.
*   **`FullConfig`**: This is the comprehensive, internal configuration model that is generated by the `ConfigExpander`. It will be a deeply nested model containing explicit, validated settings for every single parameter used by every module. This includes DFT parameters (`ecutwfc`, `kpoints_density`), MLIP hyperparameters (`body_order`, `r_cut`), and active learning settings (`max_generations`, `uncertainty_threshold`). Having this single, immutable `FullConfig` object passed through the system ensures that every part of a workflow is executed with a consistent and reproducible set of parameters.

**Class and Interface Design:**
The architecture will follow the principles of dependency inversion to maximize testability and modularity.

*   **`cli.py`**: Will contain the main user entry point, `mlip-pipe run <input.yaml>`, using the `Typer` library for a clean and modern command-line experience. It is responsible for parsing the input file, initiating the `ConfigExpander`, and handing off control to the `WorkflowOrchestrator`.
*   **`orchestrator.py`**: The `WorkflowOrchestrator` class will be the heart of the application's logic. It will be initialized *with* instances of the database wrapper and all the engine modules (a form of dependency injection). Its methods, such as `run_full_pipeline()`, will contain the high-level state machine logic for the active learning loop, but not the low-level implementation details of any specific task.
*   **`database.py`**: The `AseDBWrapper` class will provide a high-level, business-logic-aware API for all database operations. Instead of raw SQL-like queries, it will expose methods like `add_structures_to_labeling_queue()`, `get_next_unlabeled_structure()`, and `mark_structure_as_labeled(id, results)`. This abstracts the underlying database technology and centralizes data access logic.
*   **Modules**: Each module will be implemented as a class (e.g., `LabelingEngine`) with a primary public method that conforms to a consistent interface (e.g., `run()`). These modules will be stateless and will depend only on the `FullConfig` object and the `AseDBWrapper`, not on each other. This strict separation ensures that, for example, the `LabelingEngine` can be tested in complete isolation without needing to run a real `TrainingEngine`.

## 5. Implementation Plan

The project will be developed in two sequential, strategic cycles. This phased approach allows us to build and validate a foundational core before adding the more complex automation layers, reducing risk and ensuring a stable base to build upon.

**Cycle 1: Core Data Pipeline & Manual Workflow**
This foundational cycle focuses on creating a functional, albeit manually driven, pipeline. The primary goal is to establish the core data flow from a user-provided input structure to a trained potential, proving that the main components can interact correctly through the centralized database. The intelligent, autonomous features of the system are intentionally deferred to the next cycle. At the end of this cycle, a knowledgeable user will be able to manually add structures to the database and then use the command-line interface to trigger the labeling and training processes to generate a valid MLIP. This cycle essentially builds the robust chassis and powertrain of the car, ensuring all the fundamental mechanics are in place before we add the autopilot system. This approach allows us to tackle the most significant external dependencies—the DFT code and the MLIP training library—in a controlled environment.
*   **Features**:
    1.  **Project Scaffolding**: Set up the `pyproject.toml` using `uv`, define the core dependencies (`ase`, `pydantic`, `typer`), and establish the `src/mlip_autopipec` directory structure.
    2.  **Configuration Models (`FullConfig`)**: Implement the comprehensive `FullConfig` Pydantic model and its nested sub-models. The system will initially only work with this detailed configuration file, which the user must create manually.
    3.  **Database Wrapper**: Create the `AseDBWrapper` class to provide a clean, high-level API for all interactions with the project's SQLite database, abstracting away the raw `ase.db` calls.
    4.  **Module C (Labeling Engine)**: Implement the core functionality to run a single-point DFT calculation using Quantum Espresso. This includes robust input file generation from an `ase.Atoms` object, executing `pw.x` via a subprocess, and parsing the output file to extract energy, forces, and stress. Basic error detection (e.g., non-zero exit code) will be included.
    5.  **Module D (Training Engine)**: Implement the ACE model training logic, including the delta learning strategy. This module will query the database for all available labeled data, prepare it for the training library, execute the training, and save the resulting potential file.
    6.  **Basic Orchestrator**: A simplified `WorkflowOrchestrator` that can execute the labeling and training steps in a fixed, sequential order.
    7.  **Basic CLI**: A functional CLI with distinct commands for each step, such as `mlip-pipe label` and `mlip-pipe train`, allowing a user to manually drive the pipeline.

**Cycle 2: Automation, Active Learning & Advanced Simulation**
This cycle builds directly upon the stable foundation laid in Cycle 1 to deliver the full "human-out-of-the-loop" promise. It introduces the intelligent modules that automate the creative and decision-making aspects of potential generation. The active learning loop, which is the core of the system's efficiency and power, is the central feature of this cycle. This cycle transforms the manually-operated tool from Cycle 1 into a truly autonomous scientific instrument.
*   **Features**:
    1.  **`MinimalConfig` and `ConfigExpander`**: Develop the simple, user-facing `MinimalConfig` and the heuristic engine that intelligently expands it into the `FullConfig`, abstracting away the complexity from the end-user.
    2.  **Module A (Structure Generator)**: Implement the various initial structure generation techniques (SQS, NMS, rattling) to provide a diverse and physically meaningful starting point for the workflow.
    3.  **Module B (Explorer & Sampler)**: Integrate the MACE universal potential. Implement the logic to perform large-scale MD simulations, calculate geometric descriptors, and apply clustering and sampling algorithms to select the most informative structures for labeling.
    4.  **Module E (Simulation Engine)**: Integrate the newly trained ACE potential with a simulation engine (e.g., ASE's MD). The key feature is the implementation of the on-the-fly uncertainty monitoring to detect when the simulation has entered an unknown region of chemical space.
    5.  **Full Orchestrator**: Heavily enhance the orchestrator to manage the complete, asynchronous active learning loop: `Generate -> Explore -> Label -> Train -> Simulate -> Detect Uncertainty -> Loop`. This involves significant state management logic.
    6.  **Advanced Error Recovery**: Augment the `LabelingEngine` (Module C) with the sophisticated, multi-stage error recovery mechanisms (e.g., automatically retrying failed DFT calculations with safer parameters).
    7.  **Full CLI**: Finalize the user-facing CLI, making the simple `mlip-pipe run input.yaml` command the primary and default entry point for the entire autonomous workflow.

## 6. Test Strategy

A rigorous, multi-layered testing strategy is essential for ensuring the reliability of such a complex, autonomous system. Testing will be a critical part of both cycles, with a focus on unit tests for component-level correctness, integration tests for inter-component communication, and end-to-end tests for workflow validation.

**Cycle 1 Test Strategy:**
The testing focus in the first cycle will be on verifying the correctness of the core data modules and ensuring the integrity of the manual data pipeline from end to end.
*   **Unit Tests**: Each class and its methods will be tested in strict isolation using `pytest` and `pytest-mock`. For the `LabelingEngine`, this is particularly critical: we will thoroughly mock the `subprocess.run` call to Quantum Espresso. One test will simulate a successful run, providing a sample output file and asserting that the parser extracts the correct floating-point values for energy, forces, and stress. Another test will simulate a DFT failure by having the mock return a non-zero exit code, and we will assert that our engine catches this and raises a specific, informative exception. The Pydantic configuration models will also have dedicated tests to ensure their validation logic (e.g., rejecting negative cutoff energies) is working as expected.
*   **Integration Tests**: The primary integration test for Cycle 1 will be a "pipeline" test that verifies the correct interaction between the orchestrator, the database, and the core modules. This test will start with a known atomic structure (e.g., a 2-atom Si cell), and run it through the orchestrator's `label` and `train` processes. To keep the test fast and independent of an actual DFT installation, we will mock the `subprocess.run` call at the integration level, having it return a pre-computed, valid output string. This allows us to test the entire chain of custody—from file parsing, to database insertion, to data retrieval, to training preparation—verifying that the components are wired together correctly and that data flows between them as expected.

**Cycle 2 Test Strategy:**
In the second cycle, the testing strategy will expand significantly to cover the complex logic of the new automation modules and the state transitions of the active learning loop.
*   **Unit Tests**: All new modules (A, B, E) will receive comprehensive unit tests. For the `StructureGenerator`, we will assert that the output structures have the correct composition and properties. For the `Explorer & Sampler`, we will test the descriptor calculations and sampling algorithms on known data distributions to ensure they select diverse candidates. The most important new unit tests will be for the `SimulationEngine`. We will create a mock MLIP calculator that is programmed to return a high uncertainty value at a specific, predetermined timestep. We will then assert that the `SimulationEngine` runs for exactly that many steps and correctly extracts the right structure, thus validating the core trigger mechanism of the active learning loop.
*   **Integration Tests**: We will add several new integration tests. A key test will verify the interaction between the `SimulationEngine` and the `Explorer & Sampler`, ensuring that an "uncertainty" flag raised by the engine correctly triggers the selection and queuing of a new structure for labeling in the database. Another test will focus on the `ConfigExpander`, feeding it a `MinimalConfig` and asserting that the generated `FullConfig` contains a complete and physically sensible set of parameters.
*   **End-to-End (E2E) Test**: The capstone of our test strategy will be a full E2E test of the autonomous system. This test will start from a `minimal.yaml` file and execute a simplified, fast version of the entire active learning loop for one or two generations. This will use the same mocking strategy as the Cycle 1 integration test (mocking the DFT and training calls) to ensure it can run quickly in a CI environment. The purpose of this test is not to produce a good potential, but to validate that the `WorkflowOrchestrator`'s complex state machine is correct and that the data flows properly through the entire, recursive, closed-loop system. This is the ultimate validation of the "human-out-of-the-loop" objective.
