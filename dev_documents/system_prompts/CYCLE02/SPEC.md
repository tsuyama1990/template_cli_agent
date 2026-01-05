# Specification: CYCLE02 - Advanced Exploration Engine and Intelligent Sampling

## 1. Summary

This document provides the detailed technical specification for the second and final development cycle of the MLIP-AutoPipe project. Building upon the solid architectural foundation and the fully functional, albeit simplified, pipeline established in Cycle 1, this cycle is dedicated to implementing the sophisticated, scientific core of the application. The primary objective is to replace the placeholder components with their full, computationally intensive implementations, thereby transforming the tool from a basic structure generator into an intelligent and powerful framework for creating high-quality, diverse training datasets for Machine Learning Interatomic Potentials (MLIPs). This cycle represents the culmination of the project's vision, delivering the key features that provide its unique scientific value and distinguish it from simpler, less robust data generation scripts. The work in this cycle is more complex and computationally focused, directly tackling the challenges of simulating realistic material behaviour and making intelligent, data-driven decisions about which atomic configurations are most valuable for training a machine learning model.

The central focus of this cycle will be the implementation of the `MDExplorer`, the computational engine that uses Molecular Dynamics (MD) and hybrid Monte Carlo (MC) simulations to explore the potential energy surface of a material. This is the most critical and complex component of the entire system. Its development will involve the careful integration of external MLIP models (such as MACE or SevenNet) to calculate energies and forces, the management of the simulation dynamics using the robust integrators provided by the Atomic Simulation Environment (ASE), and the implementation of advanced, physics-aware features. These features include the logic for automatically switching between different thermodynamic ensembles based on the physical properties of the system and the integration of the ZBL potential to handle the high-energy atomic interactions that are common in high-temperature simulations. Furthermore, this cycle will replace the simple random sampler from Cycle 1 with an intelligent Farthest Point Sampling (FPS) algorithm. This advanced sampler will use structural descriptors to select the most geometrically diverse configurations for the final dataset, ensuring maximum information content and minimal redundancy. To enhance usability and make the tool more accessible, especially to new users, we will also introduce a simple, interactive web-based user interface for configuration and visualization. By the end of this cycle, MLIP-AutoPipe will be a feature-complete, production-ready, and scientifically rigorous tool that fully realises the project's ambitious goal of automating the generation of superior training data for next-generation MLIPs.

## 2. System Architecture

The architecture in Cycle 2 is a direct and planned extension of the foundation laid in Cycle 1. The modular design established in the first cycle allows us to seamlessly replace the placeholder components (`MDExplorer`, `RandomSampler`) with their full implementations without requiring significant changes to the core orchestration logic. This demonstrates the value of the initial architectural planning. In addition to these replacements, we will be adding new modules to support the advanced physics calculations and the new user interface. The changes are targeted and designed to fit cleanly into the existing structure, preserving the principles of high cohesion and low coupling. The focus is on enriching the capabilities of the system while maintaining its structural integrity and clarity. The blueprint below highlights the specific files that will be the focus of development in this cycle.

**File Structure for Cycle 2:**

The following ASCII tree shows the files that will be created or significantly modified in this cycle. Files marked in **bold** represent the primary deliverables for this phase of the project. This structure shows the evolution of the codebase as it matures from a basic pipeline into a fully-featured scientific application.

```
.
├── src
│   └── mlip_autopipec
│       ├── **web_ui.py**
│       ├── core
│       │   └── **physics.py**
│       ├── exploration
│       │   └── **engine.py**
│       ├── generators
│       │   └── **ionic.py**
│       └── sampling
│           └── **samplers.py**
└── tests
    ├── **test_physics.py**
    ├── **test_exploration.py**
    └── **test_integration.py**
```

**Component Blueprint:**

*   **`exploration/engine.py`**: This file, which contained a simple placeholder in Cycle 1, will be completely overhauled to become the computational heart of the application. The placeholder `MDExplorer` will be replaced with a full implementation that can run actual MD and hybrid MD/MC simulations.
    *   It will dynamically load a user-specified MLIP calculator (e.g., MACE) using ASE's `get_calculator` interface.
    *   It will use ASE's `Langevin` or other dynamics integrators to run the simulation for a specified number of steps at a given temperature.
    *   It will contain the logic for performing hybrid MD/MC moves, specifically the `SwapMove`, which is essential for exploring chemical ordering in alloys.
    *   It will implement the crucial logic for automatic NVT/NPT ensemble switching by calling the `detect_vacuum` function from the `core.physics` module.
    *   It will include the advanced logic to mix in the ZBL potential for handling short-range atomic repulsions, which is vital for preventing simulation failures at high temperatures.

*   **`core/physics.py`**: This is a new module that will be created to house sophisticated, reusable physics-based utility functions that are too complex for the main orchestrator module.
    *   Its primary initial component will be the `detect_vacuum` function. This function will implement a robust, grid-based algorithm to determine if an atomic structure contains a vacuum layer (and is therefore a slab or surface) or if it is a fully periodic bulk material. This information is critical for selecting the correct thermodynamic ensemble.

*   **`sampling/samplers.py`**: This file will be significantly updated to include the new, intelligent `FPSSampler`, which will replace the basic `RandomSampler` as the default choice.
    *   The `FPSSampler` will be implemented. This will involve integrating a third-party library (e.g., `dscribe`) to compute SOAP (Smooth Overlap of Atomic Positions) descriptors for each atomic environment in every structure generated during the exploration phase.
    *   The core FPS algorithm will then be implemented to iteratively select a subset of structures that are farthest apart in this high-dimensional SOAP feature space, thus maximizing the diversity of the final dataset.

*   **`generators/ionic.py`**: A new generator will be added to the `generators` package to extend the system's capabilities to a new class of materials.
    *   The `IonicGenerator` will extend the `BaseStructureGenerator` abstract class. It will include specific logic to handle the oxidation states of the elements and to ensure that the final generated crystal structures are charge-neutral, a fundamental physical requirement.

*   **`web_ui.py`**: This new file will contain a simple, user-friendly web-based graphical user interface.
    *   It will be built using a simple, lightweight web framework like Streamlit or Gradio, chosen for its ease of development and deployment.
    *   It will provide a graphical way for users to define the `SystemConfig` parameters, start a pipeline run in the background, and view visualizations of the generated structures from the output database. This will make the tool much more accessible to novice users.

*   **`tests/test_exploration.py`**: A new test file dedicated to unit testing the complex setup and decision-making logic within the `MDExplorer`.
*   **`tests/test_physics.py`**: A new test file for unit testing the `detect_vacuum` algorithm with a variety of known input structures.
*   **`tests/test_integration.py`**: A new and critically important test file for the end-to-end integration tests, which are a key deliverable of this cycle, proving that all the new, complex components work together correctly.

## 3. Design Architecture

The design in Cycle 2 is sharply focused on ensuring robustness, performance, and scientific correctness, particularly within the computationally intensive `MDExplorer` and `FPSSampler` components. The key design principle is the careful management of state, external dependencies (like ML models), and computational resources. The architecture is designed to be both powerful enough for expert users and simple enough to be used by newcomers, with complexity hidden behind clean, well-defined interfaces. The design choices prioritize long-term maintainability and the ability to withstand the rigours of real-world scientific research, where simulations can be unpredictable and computationally demanding.

**MDExplorer Design:**

*   **Late Binding of Calculator:** The `MDExplorer` will not instantiate the computationally expensive MLIP calculator in its `__init__` method. Instead, the calculator will be loaded lazily, on-demand, inside the method that actually runs the simulation. This is a critical design choice for performance and memory management, especially when the explorer is run in a parallel processing environment. It avoids the need to pickle and transfer the large ML model (which can be hundreds of megabytes) between processes and ensures that GPU resources, if used, are managed correctly within each worker process.
*   **Stateful Dynamics and Trajectory Handling:** The class will encapsulate an ASE `Langevin` (or other) dynamics object. It will be responsible for attaching the calculator and the atomic structure to the dynamics object and running the simulation for the specified number of steps. To capture the results, it will attach an `Observer` to the dynamics object. This observer will be a callable that is invoked at regular intervals during the simulation, and its role will be to save the state of the system (positions, forces, energy) to a trajectory file.
*   **Robust Error Handling:** The simulation runner method will be wrapped in a comprehensive `try...except` block to gracefully handle potential `PhysicsViolationError` exceptions that can be raised by ASE (e.g., if atoms get too close during a high-temperature simulation, causing forces to become unstable). In case of such an error, the explorer will log the issue clearly, discard the failed structure, and continue its work on the other structures, ensuring that the failure of a single simulation does not cause the entire pipeline to crash.

**Intelligent Sampling (FPS) Design:**

*   **Decoupling of Descriptors and Sampling:** The `FPSSampler` will be designed in two logical parts. The first part will be responsible for the computationally intensive task of calculating the feature vectors (the SOAP descriptors) for all structures generated by the exploration stage. The second part will be the implementation of the core FPS algorithm itself, which will operate on these pre-calculated vectors. This separation makes the code cleaner and easier to test.
*   **Dependency Management:** The implementation will require adding new, specialized scientific dependencies to the `pyproject.toml` file for computing the SOAP descriptors (e.g., `dscribe`). Care will be taken to choose a library that is well-maintained, performs efficiently, and integrates well with the existing ASE-based workflow.
*   **Consistent API:** The `FPSSampler` will expose a simple `sample(structures: list[ase.Atoms])` method, which is consistent with the sampler interface defined by the placeholder in Cycle 1. Internally, it will manage the complex process of computing the descriptors and executing the sampling algorithm, hiding this complexity from the `PipelineRunner`.

**Vacuum Detection Algorithm Design:**

*   **Grid-Based Physical Approach:** The `detect_vacuum` function in `core/physics.py` will work by constructing a 3D grid of points over the simulation cell. This approach is more robust than simply checking the cell dimensions, as it can handle non-orthogonal cells and complex geometries.
*   **Distance Calculation and Thresholding:** For each point in this grid, the algorithm will calculate the distance to the nearest atom in the structure. If a contiguous region of grid points is found where all points are farther than a certain threshold distance from any atom, this region is considered to be a vacuum. The function will then analyse the dimensionality and connectivity of this vacuum region to distinguish between a 1D vacuum (e.g., a void), a 2D vacuum layer (indicating a slab geometry), or a fully periodic system (indicating a bulk material). This classification is then used to make the crucial decision about the thermodynamic ensemble.

## 4. Implementation Approach

The implementation of Cycle 2 will prioritise the most complex and highest-risk components first. This strategy ensures that the most significant challenges are tackled early in the cycle, leaving ample time for testing and refinement. The development process will continue to follow a TDD-like methodology, with new unit and integration tests being developed in parallel with the new features.

1.  **Vacuum Detection Algorithm:** The implementation will begin with the `detect_vacuum` function in `core/physics.py`. This is a self-contained and purely algorithmic piece of logic that can be developed and unit-tested in complete isolation from the rest of the system. This makes it an ideal starting point. It will be thoroughly tested with a variety of known structures before being integrated.
2.  **`MDExplorer` Implementation:** This is the largest and most central task of the cycle.
    *   First, the basic logic for loading a calculator (using a simple, fast potential like EMT for initial testing) and running a standard MD simulation using ASE's `Langevin` dynamics will be implemented.
    *   Next, the automatic ensemble switching feature will be integrated. This will involve calling the newly created `detect_vacuum` function to make the decision between using an NVT or NPT ensemble for each input structure.
    *   The hybrid MD/MC `SwapMove` will be added. This will likely involve creating a custom ASE dynamics move class that can be attached to the main dynamics object to perform atom swaps periodically.
    *   Finally, the ZBL potential mixing logic will be implemented. This is an advanced feature that may require creating a custom ASE calculator class that internally holds both the MLIP and ZBL calculators and returns a combined energy and force.
3.  **`FPSSampler` Implementation:**
    *   A suitable library for SOAP descriptor calculation will be researched, chosen, and added as a dependency to `pyproject.toml`.
    *   The `FPSSampler` class will be implemented in `sampling/samplers.py`. The logic for calculating the SOAP descriptors for all input structures will be written first. This will be followed by the implementation of the core FPS algorithm, which iteratively selects new structures based on their distance in the SOAP feature space.
4.  **New `IonicGenerator`:** The `generators/ionic.py` module will be created. The `IonicGenerator` will be implemented, extending the `BaseStructureGenerator`. This will include adding logic to look up the common oxidation states of the specified elements and ensuring that the final generated structure is charge-neutral.
5.  **Web UI:** The `web_ui.py` file will be created. A simple UI will be built using Streamlit. It will provide widgets for setting the key configuration parameters (e.g., elements, temperature). A "Run" button will be implemented to trigger the execution of the `PipelineRunner` in a separate thread or subprocess to avoid blocking the UI. An area will be designated for displaying a 3D visualization of one of the generated structures using a library like `py3Dmol`.
6.  **Updating the Orchestrator:** The `PipelineRunner` in `core/orchestrator.py` and the `cli.py` will require minor updates to accommodate the new samplers (`fps`) and generators (`ionic`), for example, by adding them to a factory or to the conditional logic that selects which component to use based on the configuration.
7.  **Testing:** Throughout the process, corresponding unit tests will be developed in `test_exploration.py` and `test_physics.py`. Once all the individual components are implemented, the final and most crucial step will be to write the end-to-end integration tests in `test_integration.py`, which will verify that the entire, complex pipeline works as expected.

## 5. Test Strategy

The testing strategy for Cycle 2 shifts its focus from individual unit tests to more comprehensive **integration tests**, while still ensuring that new, complex algorithms are thoroughly unit-tested in isolation. The goal of this phase of testing is to ensure not only that each component works correctly on its own, but that they work together seamlessly to produce a scientifically valid result. Given the computational expense of the new components, the test strategy must be carefully designed to be both thorough and efficient.

**Unit Testing Approach:**

*   **`test_exploration.py`**:
    *   **Goal:** To verify the complex configuration and setup logic of the `MDExplorer` without running expensive and time-consuming simulations. The focus is on the "control plane" of the explorer, not its full data plane.
    *   **Test Cases:**
        1.  **Ensemble Selection Logic:** A key test will be written to verify the automatic ensemble switching. It will provide the `MDExplorer` with a pre-defined "bulk" structure and another "slab" structure. By mocking the `detect_vacuum` function to return known, deterministic values, the test will assert that the explorer correctly selects the NPT ensemble for the bulk system and the NVT ensemble for the slab system.
        2.  **Calculator Loading and Configuration:** A test will be written to ensure the explorer correctly loads and configures the specified MLIP model. This can be done by mocking ASE's `get_calculator` function and asserting that it was called with the correct model name and parameters from the configuration file.
        3.  **Graceful Error Handling:** A test will be created where the mocked ASE dynamics object is forced to raise a `PhysicsViolationError`. The test will assert that the `MDExplorer` correctly catches this exception, logs an informative message, and continues its execution without crashing the entire pipeline.

*   **`core/physics.py` Tests (`test_physics.py`)**:
    *   **Goal:** To rigorously verify the correctness of the `detect_vacuum` algorithm.
    *   **Test Cases:** A series of tests will be created, each with a different, pre-constructed `ase.Atoms` object representing a known physical system (e.g., a bulk FCC crystal, a slab with a 10 Angstrom vacuum layer, a single molecule in a large box, a porous material). A unit test will be written for each of these cases, asserting that the `detect_vacuum` function returns the correct classification.

**Integration Testing Approach:**

*   **`test_integration.py`**:
    *   **Goal:** To verify that the entire, feature-complete pipeline works from end-to-end, including the computationally intensive steps. A critical aspect of these tests is that they will use a fast, simple, and deterministic potential (like ASE's built-in EMT potential) instead of a real MLIP. This allows the tests to run efficiently and reproducibly within a CI/CD environment.
    *   **Test Cases:**
        1.  **Full MD Run Verification:** An integration test will be created for a small (e.g., 8-atom) copper system. The test will configure the pipeline to run a short MD simulation (e.g., 50 steps). The assertions will be multi-faceted:
            *   The pipeline must complete successfully without any errors.
            *   The final database must contain the correct number of structures.
            *   Crucially, the average of the final atomic positions in the database must be measurably different from the average of the initial positions, providing definitive proof that the MD simulation actually ran and modified the geometry.
        2.  **FPS Sampling Logic Verification:** A test will be designed to verify the core logic of FPS. It will generate a few very distinct initial structures (e.g., a compressed lattice, an expanded lattice, and a sheared lattice). It will then run the pipeline with the `FPSSampler` configured to select only two structures. The test will assert that the final database contains the two structures that are most different from each other (the compressed and expanded ones), thus verifying that the FPS logic is working as intended.
        3.  **Hybrid MC Swap Move Verification:** An integration test will be created for a small binary alloy (e.g., a 2-atom CuAu cell) and will configure the pipeline to use the hybrid MD/MC `SwapMove`. After the run, the test will load the final structure and assert that the chemical identities of the atoms at positions 0 and 1 have been swapped, confirming that the custom MC move was successfully executed during the dynamics.
