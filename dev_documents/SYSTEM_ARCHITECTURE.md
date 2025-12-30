# System Architecture: MLIP-AutoPipe

## 1. Summary

The MLIP-AutoPipe (Machine Learning Interatomic Potential - Automated Pipeline) system is an advanced, autonomous platform designed to revolutionise materials science research by automating the entire workflow for generating and validating bespoke Machine Learning Interatomic Potentials (MLIPs). The core philosophy of the project is to "remove the human expert from the loop," addressing the most significant bottleneck in computational materials science: the reliance on expert intuition and manual intervention for generating high-quality training data for MLIPs. The creation of a reliable potential has traditionally been an artisanal process, demanding a profound understanding of quantum mechanics, extensive hands-on experience with Density Functional Theory (DFT) calculations, and a laborious, iterative process of selecting appropriate atomic structures and simulation conditions. An expert researcher must make dozens of critical decisions: which crystal phases to include, what types of defects are relevant, the temperature and pressure ranges for molecular dynamics, and how to sample the resulting trajectories. This manual curation is not only inefficient and expensive, consuming countless hours of researcher time and supercomputer resources, but it is also inherently subjective and suffers from a lack of reproducibility. A potential developed in one research group can be difficult to replicate by another, hindering scientific progress. MLIP-AutoPipe is engineered to replace this manual, intuition-driven craft with a systematic, robust, and fully automated scientific pipeline.

The system accomplishes this through a sophisticated synthesis of physics-based heuristics, advanced sampling algorithms, and a closed-loop active learning framework. The user's interaction is intentionally minimised to the highest level of abstraction: providing the chemical composition of a material, for example, "FePt" or "SiGe". From this single piece of information, the system autonomously initiates a multi-stage workflow. It begins by intelligently generating a diverse portfolio of initial structures, employing algorithms specifically chosen for the material's bonding character. It then proceeds to a high-throughput exploration phase, using a universal surrogate potential (MACE) to run large-scale molecular dynamics simulations, scanning a vast landscape of atomic configurations at a tiny fraction of the cost of first-principles methods. Following this exploration, the system executes high-precision DFT calculations via Quantum Espresso, but only for a small, intelligently selected subset of structures deemed most informative. This labelled data is then used to train a specialised, high-fidelity MLIP using the Atomic Cluster Expansion (ACE) framework. The final, and most crucial, stage is the deployment of this newly trained potential within the system's own simulation engine. An on-the-fly uncertainty quantification mechanism continuously monitors this simulation, detects when the model encounters novel atomic configurations (its "unknown unknowns"), and automatically triggers a retraining loop. This feedback mechanism sends the new, challenging structure back to the DFT engine for labelling, augments the dataset, and refines the model. This active learning cycle ensures that computational resources are spent judiciously, drastically reducing the cost and time compared to traditional methods that rely on brute-force Ab Initio Molecular Dynamics (AIMD). The final output is not just a potential, but a complete, self-consistent package containing the high-fidelity MLIP, the full provenance of its training data, and the configurations used to generate it. This enables researchers to investigate complex phenomena such as phase transitions, defect migration, and chemical reactions over previously inaccessible time and length scales with unprecedented confidence and reproducibility.

## 2. System Design Objectives

The design of the MLIP-AutoPipe system is underpinned by a set of fundamental objectives aimed at forging a more powerful, reliable, and accessible paradigm for computational materials science. These principles guide the architecture and functionality of the entire framework.

**1. Maximise Automation and Eliminate Human Bottlenecks:** The primary objective is to engineer a truly "zero-touch" pipeline. The system must automate every conceivable step of the workflow, from interpreting the user's high-level input to generating the final, validated potential and its associated report. This automation transcends simple scripting of existing codes. It involves encapsulating the cognitive processes and decision-making of a human expert. This includes the automated classification of the material's chemical bonding to select the most physically appropriate structure generation algorithm (e.g., SQS for disordered alloys, NMS for molecules). It also involves the complete automation of the DFT calculation process, which traditionally requires expert supervision to select pseudopotentials, test for k-point and energy cutoff convergence, and troubleshoot self-consistent field (SCF) convergence failures. The system will implement robust error-recovery protocols for DFT calculations, automatically adjusting parameters like the mixing beta or smearing temperature if a calculation fails. The ultimate goal is the complete removal of manual data curation, parameter tuning, and convergence checking, which are collectively the most time-consuming and expertise-dependent stages in the current MLIP creation process. By achieving this, the system democratises the creation of high-quality potentials, making them accessible to researchers who may not be world experts in electronic structure theory.

**2. Guarantee Physical Plausibility and Robustness:** While full automation is the primary goal, the generated potential must be physically meaningful, robust, and transferable. The system is designed to embed physical principles at every stage of the workflow to constrain the machine learning model and prevent unphysical behaviour. This philosophy is first embodied in the "Two-Tier Configuration" strategy, where a minimal user input is expanded into a comprehensive execution plan by a physics-based heuristic engine. This engine infers appropriate simulation temperatures from melting point estimations and selects DFT parameters based on established protocols like SSSP. The core of the physical robustness comes from the "Delta Learning" approach. Instead of learning the entire potential energy surface from scratch, the MLIP is trained to learn the *difference* between a simple, analytical baseline potential (like Lennard-Jones or ZBL) and the true DFT energy landscape. This simple baseline correctly describes the strong repulsive interactions at very short interatomic distances, preventing atoms from unphysically overlapping. By learning only the complex residual, the MLIP is free to focus its capacity on the chemically relevant bonding region, leading to a more accurate and stable model. Furthermore, during the active learning phase, advanced boundary condition treatments will be applied to structures extracted from periodic simulations, ensuring that artefacts from artificial periodicity are meticulously handled to produce a potential that is accurate for both bulk and defect systems.

**3. Optimise Computational Efficiency:** The generation of training data via DFT calculations is, by a significant margin, the most computationally expensive part of the workflow. A core design objective is therefore to minimise the number of required DFT calculations to the greatest extent possible without sacrificing the accuracy or transferability of the final potential. The system attacks this challenge with a multi-pronged strategy. The first line of attack is the use of a pre-trained, universal surrogate potential (MACE). This model, trained on a massive database of existing materials, can predict energies and forces with reasonable accuracy at a computational cost that is orders of magnitude lower than DFT. The system leverages this to rapidly explore the potential energy surface, generating millions of candidate structures. The second part of the strategy is intelligent down-selection. Instead of running DFT on all of these candidates, the DIRECT sampling algorithm is employed. It uses structural descriptors to identify and cluster similar configurations, and then selects only a few diverse, representative structures from each cluster. This avoids wasting computational effort on redundant structures and ensures the initial training set is not biased towards a single region of phase space. The final and most powerful efficiency optimisation is the on-the-fly active learning loop. This ensures that subsequent DFT calculations are *only* performed when the model itself signals that it is uncertain, thereby focusing expensive resources exclusively on the task of patching the most significant gaps in the model's knowledge.

**4. Ensure Modularity, Extensibility, and Reproducibility:** The fields of machine learning and materials modelling are evolving at a breathtaking pace. To avoid obsolescence, the MLIP-AutoPipe is architected as a collection of modular, interchangeable components. The five principal modules (Structure Generator, Explorer, Labeling Engine, Training Engine, Simulation Engine) are designed with clean, well-defined interfaces. This modularity allows for future extensions and adaptations. For example, the Quantum Espresso `LabelingEngine` could be swapped for a different DFT code like VASP or SIESTA by simply writing a new class that conforms to the engine's interface. Similarly, a new and improved MLIP architecture, like an equivariant graph neural network, could be slotted into the `TrainingEngine` to replace ACE. This "pluggable" architecture ensures the long-term viability and relevance of the platform. Crucially, the system is designed to be fully reproducible. Every action, from the initial expansion of the configuration file to every DFT calculation and every retraining step, is meticulously logged. The final output is a comprehensive bundle that includes not just the trained potential, but also the complete training dataset stored in the ASE database, and the final `exec_config_dump.yaml` file that specifies every parameter used in the workflow. This guarantees that any result generated by the pipeline is transparent, auditable, and can be precisely replicated by other researchers.

## 3. System Architecture

The MLIP-AutoPipe system is designed as a modular pipeline, orchestrated by a central workflow manager. The five core modules interact sequentially and, in the case of the active learning loop, cyclically. Data is persisted and passed between modules via a central ASE (Atomic Simulation Environment) database, which acts as the single source of truth for all generated structures and their calculated properties. This database-centric design decouples the modules, allowing them to operate asynchronously and making the overall system more robust.

```mermaid
graph TD
    subgraph User Input
        A[input.yaml];
    end

    subgraph Orchestrator
        B[Workflow Manager];
    end

    subgraph Core Modules
        C[Module A: Structure Generator];
        D[Module B: Explorer & Sampler];
        E[Module C: Labeling Engine];
        F[Module D: Training Engine];
        G[Module E: Simulation Engine];
    end

    subgraph Data Store
        H[(ASE Database)];
    end

    subgraph Output
        I[Final MLIP];
        J[exec_config_dump.yaml];
    end

    A --> B;
    B --> C;
    C -- Initial Structures --> H;
    B --> D;
    H -- Initial Structures --> D;
    D -- Candidate Structures --> H;
    B --> E;
    H -- Candidate Structures --> E;
    E -- Labeled Data --> H;
    B --> F;
    H -- Labeled Data --> F;
    F -- Trained MLIP --> G;
    B --> G;
    G -- High-Uncertainty Structures --> H;
    G -- Simulation Trajectories --> I;
    B --> J;

    %% Active Learning Loop
    E_loop_start(( )) --> F_loop(( )) --> G_loop(( )) --> E_loop_start;
    linkStyle 11 stroke:red,stroke-width:2px,stroke-dasharray: 5 5;
    linkStyle 12 stroke:red,stroke-width:2px,stroke-dasharray: 5 5;
    linkStyle 13 stroke:red,stroke-width:2px,stroke-dasharray: 5 5;

    note right of G
        Active Learning Loop:
        The Simulation Engine (G) identifies
        high-uncertainty structures, which
        are fed back to the Labeling Engine (E)
        for DFT calculation. The new data is
        used by the Training Engine (F) to
        retrain and improve the MLIP.
    end
```

**1. Workflow Orchestrator:**
This is the central controller of the pipeline, implemented as a stateful script that reads the configuration files and invokes the modules in the correct sequence. It is responsible for the overall logic of the workflow, managing the transitions between different stages of the potential's development. Its key responsibilities include:
-   **Configuration Expansion:** Its first task is to read the user's minimal `input.yaml` and invoke the internal physics-based heuristic engine. This engine populates all the necessary default parameters, from DFT convergence settings to MD simulation temperatures, generating the complete `exec_config_dump.yaml`. This complete configuration is then used for the rest of the session, ensuring reproducibility.
-   **Module Invocation:** The orchestrator calls the entry point of each of the five core modules in the correct order, passing the necessary configuration and database connection information. It acts as the "glue" that connects the independent modules into a coherent pipeline.
-   **Active Learning Loop Management:** This is the orchestrator's most complex task. It manages the iterative cycle between the Simulation, Labeling, and Training engines. After an initial model is trained, the orchestrator passes control to the Simulation Engine. If the engine identifies a new high-uncertainty structure, the orchestrator directs this structure to the Labeling Engine, waits for it to finish, and then triggers the Training Engine to produce a refined model. It monitors the number of these retraining generations and terminates the loop when a convergence criterion is met, such as reaching the maximum number of generations or when a full simulation run completes without finding any new uncertain regions.

**2. Module A: Structure Generator:**
This module is responsible for creating a diverse, physically plausible set of initial atomic structures without performing any expensive DFT calculations. Its goal is to provide a rich starting point for the exploration phase. It first determines the material's bond type (alloy, molecular, ionic, covalent) and then employs the most suitable algorithm:
-   **Alloys:** Special Quasirandom Structures (SQS) are used to model the atomic correlations in random solid solutions. The module also systematically applies volumetric and shear strains to these structures to provide initial data on the material's elastic properties.
-   **Molecules:** For molecular systems, it uses Normal Mode Sampling (NMS) to explore the vibrational landscape around an equilibrium geometry, providing data on bond stretches and angle bends.
-   **Ionic Crystals:** An approach equivalent to Ab Initio Random Structure Searching (AIRSS) is used to discover stable and metastable polymorphs, which are crucial for materials that exhibit phase transitions.
-   **Covalent Materials:** For materials like silica, it employs "Deep Rattling" (large random displacements) and simulated melt-quench protocols using the surrogate model to generate topologically diverse amorphous and defected structures.

**3. Module B: Explorer & Sampler:**
This module performs a high-throughput, computationally inexpensive exploration of the potential energy surface.
-   **Surrogate Model:** It uses a pre-trained, universal MACE potential to run large-scale molecular dynamics simulations on the initial structures. This allows it to sample millions of distinct atomic configurations in a tiny fraction of the time a first-principles simulation would take.
-   **DIRECT Sampling:** From the millions of frames generated, it uses the DIRECT algorithm. This involves calculating a structural descriptor (like SOAP or ACE) for each configuration, performing a clustering analysis on the descriptor data, and then selecting a few hundred or thousand representative configurations. This ensures the initial training set is diverse, unbiased, and information-rich.

**4. Module C: Labeling Engine:**
This is the robust interface to the external DFT code (Quantum Espresso). It takes candidate structures from the database and calculates their energies, forces, and stresses.
-   **Parameter Automation:** It automatically sets robust and high-precision DFT parameters based on the SSSP protocol, including pseudopotentials, plane-wave cutoffs, and k-point meshes, removing the need for manual convergence studies.
-   **Automated Execution & Recovery:** It manages the execution of `pw.x` and includes sophisticated error recovery logic. If an SCF calculation fails to converge, it will automatically try a sequence of remedies, such as reducing the mixing beta or changing the diagonalization algorithm, before flagging the structure as problematic.

**5. Module D: Training Engine:**
This module takes the DFT-labeled data from the database and trains the MLIP.
-   **Delta Learning:** It implements the "delta learning" scheme, where the model learns the residual between a simple physical baseline potential and the DFT data. This improves robustness and ensures correct physical behaviour at short range.
-   **MLIP Framework:** It uses the Atomic Cluster Expansion (ACE) model, a powerful and systematically improvable descriptor, as its primary training framework.
-   **Hyperparameter Optimisation:** Before the main training run, it performs an automated, lightweight search to find optimal hyperparameters (e.g., body order, polynomial degree) for the ACE model, ensuring the model has the optimal balance of flexibility and efficiency for the target system.

**6. Module E: Simulation Engine:**
This module uses the newly trained MLIP to run large-scale simulations and drives the active learning process.
-   **On-the-Fly (OTF) Inference:** It runs molecular dynamics or kinetic Monte Carlo simulations using the MLIP. During the simulation, it continuously quantifies the model's uncertainty for the current atomic configuration.
-   **Dynamic Uncertainty Threshold:** When the uncertainty exceeds a dynamically adjusted threshold, the simulation is paused. The threshold adapts as the model becomes more accurate, starting high and gradually decreasing.
-   **Structure Extraction:** The high-uncertainty structure is extracted, using advanced periodic boundary treatments to create a physically meaningful cluster for retraining. This new structure is then passed back to the Labeling Engine, closing the autonomous, self-improving learning loop.

## 4. Design Architecture

The software is designed as a modern Python package, `mlip_autopipec`, with a strong emphasis on clean separation of concerns, testability, and maintainability. The design segregates the core scientific logic from the infrastructure code (like process execution and database interaction) and the user interface.

**1. File Structure:**
-   `src/mlip_autopipec/`: The main package directory.
    -   `main.py`: The main CLI entry point, using a library like `click` or `typer`. This is the sole interface for the user.
    -   `orchestrator.py`: Contains the primary workflow logic, managing the state of the pipeline and calling the various modules.
    -   `config/`: Modules for handling the "Two-Tier Configuration." This includes the `ConfigExpander` heuristic engine and Pydantic models for validating both the user's input and the fully-expanded configuration.
    -   `data/`:
        -   `database.py`: Defines a custom `AseDB` wrapper class. This class provides a clean, high-level API for all database interactions (e.g., `add_structure`, `get_unlabeled_candidates`). It completely abstracts the underlying `ase.db` implementation, ensuring that the rest of the application is agnostic to the database schema.
        -   `models.py`: Contains Pydantic models for ensuring structured data transfer between components, such as `DFTResult` or `TrainingConfig`.
    -   `modules/`: Each of the five core scientific modules resides in its own file, encapsulating its specific logic.
        -   `structure_generator.py`: Contains the logic for SQS, NMS, etc.
        -   `explorer.py`: Implements the MACE-based exploration and DIRECT sampling.
        -   `labeling_engine.py`: A wrapper for Quantum Espresso, handling input file generation, execution, and output parsing.
        -   `training_engine.py`: Contains the logic for training the ACE model, including the delta learning implementation.
        -   `simulation_engine.py`: Implements the OTF uncertainty-driven MD/kMC simulations.
    -   `utils/`:
        -   `dft_utils.py`: Contains helper functions for DFT parameter automation, cleanly separating this domain knowledge from the main engine logic.
        -   `numba_kernels.py`: Contains performance-critical functions (e.g., descriptor calculations, kMC rate calculations) which are accelerated with Numba's JIT compilation. This isolates the highly optimised code.
-   `tests/`: A parallel directory for unit and integration tests, mirroring the structure of the `src` directory.
-   `pyproject.toml`: The single source of truth for all project metadata, dependencies, and tool configurations.

**2. Core Classes and Design Patterns:**
-   **`ConfigExpander`:** This class is a key part of the user experience. It will implement the heuristic logic for transforming the minimal `input.yaml` into the comprehensive `exec_config_dump.yaml`. It acts as the "expert in a box" for setting up a high-quality calculation.
-   **`AseDB` Wrapper:** This class encapsulates all `ase.db.connect` calls, providing methods like `get_atoms_by_state('unlabeled')` or `update_atoms_with_dft_results(id, dft_result)`. This abstraction is crucial for testability, as it allows the entire database to be mocked in unit tests, and for maintainability, as any changes to the database schema only need to be made in this one place.
-   **Engine Abstractions (Strategy Pattern):** While the initial implementation targets Quantum Espresso and ACE, the `LabelingEngine` and `TrainingEngine` will be designed around abstract base classes or protocols. This will allow for future extensions to support other DFT codes (like VASP) or MLIP models (like NequIP). The orchestrator will work with the abstract engine interface, and the concrete implementation will be chosen at runtime based on the configuration file. This is an application of the Strategy Pattern.
-   **Decoupling from `subprocess` (Dependency Inversion):** Domain logic, particularly the `LabelingEngine`, must be decoupled from infrastructure concerns like executing external processes. An infrastructure-level `ProcessRunner` class will be created to handle the execution of shell commands. The `LabelingEngine` will depend on an *abstraction* of this runner, which will be injected into it (Dependency Injection). This design follows the Dependency Inversion Principle and makes unit testing vastly simpler, as a test can inject a mock runner that returns canned output, completely avoiding any actual `subprocess` calls.

## 5. Implementation Plan

The project is decomposed into five development cycles. This staged approach allows for incremental development and testing, ensuring that each new layer of functionality is built upon a solid and validated foundation. Each cycle delivers a significant, demonstrable enhancement to the pipeline.

**Cycle 01: Core Engine - DFT Labeling and MLIP Training**
This foundational cycle focuses on building the absolute heart of the pipeline: the ability to perform a single DFT calculation and train a basic model from that single data point. This validates the most critical components of the scientific workflow.
-   **Features:**
    -   Implement the `LabelingEngine` for Quantum Espresso. This includes robust generation of input files based on an `ase.Atoms` object and, crucially, a reliable parser for the text-based output of `pw.x` to extract energy, forces, and stress.
    -   Implement the `TrainingEngine` for the ACE model. This involves writing the code to take a dataset and correctly interface with the chosen ACE library to produce a potential file.
    -   Implement the "Delta Learning" scheme. A simple, fixed baseline potential (e.g., Lennard-Jones for a noble gas) will be implemented to prove the concept.
    -   Create the `AseDB` wrapper class for persisting and querying the state of structures.
    -   A minimal orchestrator script will be developed to chain these components together: take a single structure file, run it through the `LabelingEngine`, save the result to the database, and then train a model with the `TrainingEngine`.

**Cycle 02: Structure Generation and Configuration**
This cycle adds the intelligent "front-end" to the pipeline, automating the creation of initial training data and drastically simplifying the user's input requirements.
-   **Features:**
    -   Implement the `StructureGenerator` module. At least two distinct generation methods will be built to demonstrate the system's context-awareness (e.g., SQS for alloys and NMS for molecules).
    -   Develop the "Two-Tier Configuration" system. This involves creating the Pydantic models for the input/full configurations and implementing the `ConfigExpander` heuristic engine. The engine will be capable of inferring the material type and populating the configuration with sensible defaults.
    -   Integrate the `StructureGenerator` and `ConfigExpander` into the main workflow. The pipeline will now be able to start from a simple composition string in `input.yaml` and autonomously generate its own initial dataset.

**Cycle 03: High-Throughput Exploration and Sampling**
This cycle focuses on dramatically improving the efficiency and intelligence of the data generation process by integrating the surrogate-based exploration model. This moves the system from using static initial guesses to actively exploring the energy landscape.
-   **Features:**
    -   Integrate the pre-trained MACE model as the universal surrogate potential within the new `Explorer` module. This involves handling the loading of the model and attaching it as an ASE calculator.
    -   Implement the DIRECT sampling algorithm. This requires integrating or writing a component for calculating structural descriptors (e.g., SOAP) and then implementing the clustering and stratified sampling logic.
    -   Identify performance bottlenecks in the descriptor calculation and apply Numba's `@jit` decorator for acceleration. The code will be profiled to ensure the JIT compilation is applied to the most computationally intensive loops.
    -   Refine the main workflow to insert this exploration step between the initial structure generation and the final DFT labeling.

**Cycle 04: Active Learning and Advanced Simulations**
This cycle "closes the loop," transforming the linear pipeline into a self-correcting, autonomous system that intelligently refines itself.
-   **Features:**
    -   Implement the `SimulationEngine` capable of running MD simulations with the trained ACE potential.
    -   Integrate on-the-fly uncertainty quantification into the `SimulationEngine`. This is a highly model-dependent task that requires careful implementation based on the chosen MLIP's architecture.
    -   Implement the dynamic uncertainty threshold logic, which adjusts the system's sensitivity based on the maturity of the model.
    -   Develop the advanced periodic boundary condition treatment for extracting physically meaningful, non-periodic clusters from bulk simulations for retraining.
    -   Modify the orchestrator to manage the full active learning loop, cycling between simulation, labeling, and retraining until a convergence criterion is met.
    -   (Stretch Goal) Implement a basic version of the tiered-rate kinetic Monte Carlo (kMC) to demonstrate the potential for exploring rare events.

**Cycle 05: UI/UX, Finalisation, and Documentation**
This final cycle polishes the system, focusing on the user experience to make the tool accessible, intuitive, and easy to use for the target audience of materials scientists.
-   **Features:**
    -   Develop a professional command-line interface (CLI) using `typer`. This will provide a single, consistent entry point with clear commands (e.g., `mlip-pipe run`), options, and helpful error messages.
    -   Implement a centralised logging system and add rich progress bars for long-running tasks to provide clear feedback to the user.
    -   Generate comprehensive user documentation, including a "Getting Started" guide, a detailed specification of the `input.yaml` file, and several tutorials for different classes of materials.
    -   Perform a final code cleanup, ensuring full type hint coverage (enforced by `mypy`), adding docstrings to all public functions, and enforcing a strict code style with `ruff`.

## 6. Test Strategy

The test strategy is designed to ensure the correctness, reliability, and physical plausibility of this complex scientific workflow. It combines rigorous unit testing for individual components, integration testing for module interactions, and a small number of critical end-to-end tests to validate the entire pipeline.

**Unit Testing:**
The focus of unit testing is to verify the correctness of each individual component in complete isolation from the rest of the system. This is achieved through extensive use of mocking for all external dependencies (filesystem, database, subprocesses).
-   **`ConfigExpander`:** This class is highly testable. It will be tested with a wide variety of `input.yaml` snippets to ensure the heuristic engine correctly determines bond types, selects appropriate defaults, and produces a fully-formed, valid configuration dictionary.
-   **DFT Parsers:** The functions within the `LabelingEngine` that parse Quantum Espresso output files are critical. They will be tested against a comprehensive suite of saved, real output files. This suite will include examples of successful runs, various common failure modes (e.g., SCF non-convergence, crashes), and different versions of Quantum Espresso to ensure the parser is robust.
-   **`AseDB` Wrapper:** The database wrapper will be tested with a mocked `ase.db.connect` object. Tests will assert that high-level API calls like `get_atoms_by_state` are translated into the correct underlying database queries with the correct parameters.
-   **Numba Kernels:** The JIT-compiled functions will be tested as pure numerical functions. Their output for simple, known atomic configurations will be compared against the output of a well-established, trusted reference implementation to verify their numerical correctness to within a tight tolerance.

**Integration Testing:**
Integration tests are designed to verify that different components of the system work together correctly, focusing on the interfaces and data contracts between them. These tests will use a real (but temporary) ASE database but will still mock the most expensive external processes.
-   **Generation -> Labeling:** A key test will run the `StructureGenerator`, save its output to a temporary database, and then have the `LabelingEngine` read from that database and successfully generate a valid Quantum Espresso input file. This verifies that the `Atoms` objects created by the generator are correctly formatted and can be consumed by the labeling module.
-   **Labeling -> Training:** Another test will pre-populate a database with mock DFT results (e.g., read from a text file). It will then run the `TrainingEngine` to ensure it can correctly read the data, process it through the delta learning scheme, and execute a (mocked) training run. This validates the data flow required for training.
-   **Surrogate Workflow (Module B -> C):** A test will run the MACE exploration and DIRECT sampling, saving candidate structures to the database. It will then verify that the `LabelingEngine` can correctly process these candidates, ensuring compatibility between the exploration and labeling stages.

**End-to-End (E2E) Testing:**
E2E tests will validate the entire pipeline for a small, well-defined system where the expected outcome is known. These tests are the final check that all components work together as intended in a production-like scenario.
-   **"Dry Run" Test:** An E2E test will be created that runs the entire pipeline with the `dft_compute.command` mocked to simply `echo "success"`. This test checks no numerical results but verifies that the entire workflow can execute from `input.yaml` to a final "trained" model file without crashing, ensuring all file paths, configurations, and module interactions are correctly wired.
-   **Simple System Validation (e.g., Argon Dimer):** A full E2E test will be run for a very simple physical system like an Argon dimer. The `StructureGenerator` will create stretched and compressed bonds. The `LabelingEngine` will call a real (but very fast) DFT calculation. The `TrainingEngine` will train a simple model. The test will then assert that the resulting potential energy curve has the correct qualitative shape (e.g., a single minimum at the known equilibrium distance). This provides confidence that the entire, real pipeline is physically sound.
-   **Active Learning Loop Test:** A dedicated E2E test will be designed to verify the core logic of the self-improving loop. It will start with a model trained on a single data point. The `SimulationEngine` will then be run, and the test will mock the uncertainty calculation to report high uncertainty immediately. It will then assert that the system correctly extracts the new structure, calls the `LabelingEngine` again, and successfully retrains the model. This validates the orchestrator's control flow for the most complex part of the system.