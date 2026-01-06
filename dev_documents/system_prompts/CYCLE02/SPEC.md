# Specification: MLIP-AutoPipe Cycle 2

## 1. Summary

This document provides the comprehensive technical specification for the second implementation cycle of the MLIP-AutoPipe project. Building upon the robust foundational command-line tool established in Cycle 1, this cycle is strategically dedicated to implementing the advanced, scientifically sophisticated features that will define MLIP-AutoPipe as a premier, intelligent data generation framework in the materials science domain. The primary objectives for this cycle are threefold: to dramatically enhance the physical realism and exploratory power of the simulations, to significantly increase the diversity and information content of the generated datasets, and to vastly improve the overall usability and accessibility of the tool by introducing an intuitive graphical user interface. This cycle transforms the tool from a functional utility into a powerful research instrument.

The key enhancements are strategically distributed across the pipeline's core services. The **Generation** service will be substantially expanded beyond simple alloys to encompass a wider variety of complex physical systems. We will implement new, specialized generators, including an `IonicGenerator` that programmatically enforces charge neutrality for materials like oxides and salts, and an `InterfaceGenerator` capable of constructing heterostructures and slab models for studying surface phenomena. This expansion will significantly broaden the applicability of the tool to a much wider and more diverse range of materials science research problems, from catalysis to battery materials.

The **Exploration** service will undergo the most profound transformation. The pure Molecular Dynamics (MD) engine from Cycle 1 will be evolved into a state-of-the-art hybrid MD/MC (Monte Carlo) engine. This hybrid approach will introduce powerful MC moves, such as atomic swaps and vacancy hops, which are absolutely crucial for efficiently exploring the vast compositional and configurational phase space, especially in multi-component alloys and systems with point defects. Furthermore, a cornerstone of this cycle is the integration of support for modern Machine Learning Interatomic Potential (MLIP) models, with an initial focus on the highly accurate MACE architecture. This integration includes a critical, non-negotiable feature: the on-the-fly mixing of the MLIP with a classical ZBL potential. This technique correctly handles the high-energy repulsive interactions at very short interatomic distances, a common failure point for MLIPs, thereby preventing simulation crashes and enhancing physical realism. The explorer will also be endowed with "intelligence" to automatically switch between NVT and NPT thermodynamic ensembles based on a physical analysis of the simulation cell, ensuring that bulk crystals and surface slab models are treated with the correct thermodynamics without requiring manual user specification.

The **Sampling** service will be upgraded from a simple random selector to an intelligent data curator with the implementation of a `FPSSampler` (Farthest Point Sampling). This advanced technique utilizes atomic environment descriptors (specifically, SOAP vectors) to select a maximally diverse and non-redundant subset of structures from the vast number of configurations generated during the exploration trajectories. This promises a significant improvement in the quality, efficiency, and information density of the resulting training datasets compared to the naive random sampling of Cycle 1. Finally, a major deliverable of this cycle is the development of a **Web-based User Interface (Web UI)**. This graphical interface, built with Streamlit, will provide users with an intuitive, interactive environment to build configuration files, launch pipeline runs, monitor their progress, and visualize the generated structures, dramatically lowering the barrier to entry and making the tool accessible to a broader scientific audience.

## 2. System Architecture

The architecture in Cycle 2 strategically expands upon the solid and scalable foundation laid in Cycle 1. It incorporates new modules for the advanced features and adds a distinct new entry point for the Web UI, all while preserving the core service-oriented and decoupled structure. This ensures that the new, complex components are integrated in a modular, maintainable, and testable fashion.

The `services/generators.py` file will be significantly updated with the new `IonicGenerator` and `InterfaceGenerator` classes, both inheriting from the `BaseStructureGenerator` abstract base class to ensure seamless integration with the existing factory mechanism. The most substantial changes will occur in `services/exploration.py`, which will be carefully refactored to accommodate the hybrid MD/MC logic. This includes modifying the main simulation loop and introducing a new `MixedCalculator` class to handle the MLIP/ZBL potential blending. To support the new `FPSSampler` in `services/sampling.py`, a new project dependency on a library for computing SOAP descriptors (e.g., DScribe) will be introduced and managed via `pyproject.toml`. A completely new, top-level file, `web_ui.py`, will be created to serve as the entry point for the graphical user interface. This application, built using the Streamlit framework for rapid and interactive development, will interface with the core pipeline by programmatically generating valid Pydantic configuration objects and invoking the `PipelineRunner` service in a subprocess. This demonstrates the power and flexibility of the decoupled architecture, where the core logic is entirely independent of its client.

The `PipelineRunner` itself in `orchestration.py` will require surprisingly minimal changes. The existing factory patterns for selecting generators and samplers mean that the orchestrator only needs to be able to receive a configuration object that specifies, for example, `"sampler_type": "fps"`. It will then automatically instantiate and invoke the correct new service without needing to be aware of its internal implementation details. This highlights the extensibility of the initial design. The Pydantic models in `config/models.py` will be extended to include all the new configuration options, with careful use of `Optional` fields to maintain backward compatibility with Cycle 1 configuration files.

**File Structure and Code Blueprints (Cycle 2):**

New or significantly modified files and directories in this cycle are marked in **bold**.

```
.
├── pyproject.toml
└── src
    └── mlip_autopipec
        ├── __init__.py
        ├── cli.py
        ├── **web_ui.py**             # Streamlit application entry point.
        │   └── `def main():`
        │       # 1. Build UI layout (sidebar for config, main area for results).
        │       # 2. On "Run" button click, construct FullConfig from UI state.
        │       # 3. Launch PipelineRunner in a subprocess.
        │       # 4. Display logs and visualize results from output DB.
        ├── config
        │   ├── __init__.py
        │   └── **models.py**         # Pydantic models updated with new features (MCConfig, FPSConfig, etc.).
        ├── domain
        │   ├── __init__.py
        │   └── models.py
        ├── services
        │   ├── __init__.py
        │   ├── orchestration.py
        │   ├── **generators.py**       # Add IonicGenerator, InterfaceGenerator classes.
        │   │   └── `class IonicGenerator(BaseStructureGenerator): ...`
        │   ├── **exploration.py**      # Refactor to add Hybrid MD/MC, MLIP+ZBL support.
        │   │   └── `class MixedCalculator(Calculator): ...`
        │   │   └── `def _run_single_md(...):`
        │   │       # Add logic for MC moves and Metropolis check inside MD loop.
        │   ├── **sampling.py**         # Add FPSSampler class.
        │   │   └── `class FPSSampler(BaseSampler):`
        │   │       `def _calculate_soap_vectors(...): ...`
        │   │       `def sample(...): ...`
        │   └── storage.py
        └── shared
            ├── __init__.py
            ├── ase_utils.py
            └── **physics.py**          # Add vacuum detection logic.
                └── `def detect_vacuum(atoms: ase.Atoms) -> bool:`
```

This updated structure maintains the clear separation of concerns. The Web UI acts as an alternative client to the core services, just like the CLI. The new scientific features are encapsulated within their respective service modules, minimizing the impact on the overall orchestration logic.

## 3. Design Architecture

The design in Cycle 2 meticulously extends the Pydantic schema to accommodate the full suite of new configuration options, ensuring that the new services adhere to the established abstract base class interfaces and that the entire system remains robust and validated.

**Pydantic-based Schema Design (Extensions):**

The central `FullConfig` model in `config/models.py` will be thoughtfully updated. The `Literal` types for `generator_type` and `sampler_type` will be expanded to include the new options, and entirely new Pydantic sub-models will be introduced to manage the configuration of the new, complex features.

*   `SystemConfig` (updated):
    *   `generator_type`: `Literal['alloy', 'ionic', 'interface']` - The enumeration is expanded to include the new generator types, providing compile-time safety against typos.
    *   `charge_balance`: `Optional[dict[str, int]]` - An optional dictionary mapping elements to their expected oxidation states, to be used by the `IonicGenerator`.
*   `ExplorationConfig` (updated):
    *   `use_mlip`: `bool = False` - A flag to enable the use of an MLIP model. Defaults to `False` for backward compatibility.
    *   `mlip_model_path`: `Optional[Path] = None` - The file path to the trained MACE model. This is only required if `use_mlip` is `True`.
    *   `use_zbl_mixing`: `bool = True` - A flag to enable the ZBL potential mixing, defaulting to `True` as a safe default when using MLIPs.
    *   `mc_moves`: `Optional[MCConfig] = None` - An optional sub-model to contain the configuration for all Monte Carlo moves.
*   **New `MCConfig` model:**
    *   `move_frequency`: `int = 10` - Attempt an MC move every N MD steps.
    *   `swap_probability`: `float = 0.0` - The probability of attempting an atom swap move, validated with `Field(ge=0, le=1)`.
    *   `hop_probability`: `float = 0.0` - The probability of attempting a vacancy hop move.
*   `SamplingConfig` (updated):
    *   `sampler_type`: `Literal['random', 'fps']` - The enumeration is expanded to include the 'fps' sampler.
    *   `fps_config`: `Optional[FPSConfig] = None` - A sub-model for configuring the FPS sampler, required only if `sampler_type` is `'fps'`.
*   **New `FPSConfig` model:**
    *   `n_max`: `int` - The number of radial basis functions for the SOAP descriptor.
    *   `l_max`: `int` - The maximum degree of spherical harmonics for the SOAP descriptor.
    *   `r_cut`: `float` - The cutoff radius for the local atomic environment.
    *   `sigma`: `float` - The width of the Gaussian smearing for the SOAP descriptor.

The key producer of this configuration will now be twofold: the user, who can author a more complex YAML file for the CLI, and the new Web UI, which will programmatically construct a `FullConfig` object from the user's interactive inputs in the browser. The primary consumer remains the `PipelineRunner` and its subsidiary services. The extensive use of `Optional` fields and default values is a deliberate design choice to ensure perfect backward compatibility; a user can still run a simple MD simulation using a Cycle 1-style configuration file without encountering validation errors. This schema-driven extension strategy ensures that all the powerful new features are integrated into the system in a structured, validated, and maintainable manner.

## 4. Implementation Approach

The implementation of Cycle 2 will be executed in a layered, sequential fashion, starting with the core scientific enhancements that are dependencies for other features, and concluding with the user-facing Web UI. This ensures a logical build order and allows for incremental testing.

**Step 1: Advanced Generators (`services/generators.py`)**
The first step is to broaden the tool's applicability. We will implement the `IonicGenerator` and `InterfaceGenerator` classes, ensuring they strictly adhere to the `BaseStructureGenerator` interface defined in Cycle 1. The `IonicGenerator` will require careful logic to balance the charges based on the provided oxidation states, ensuring the generated structures are electrostatically neutral. The `InterfaceGenerator` will leverage powerful functions from the Pymatgen library to create surface slabs, determine optimal orientations, and combine them into interface heterostructures.

**Step 2: SOAP Descriptors and FPSSampler (`services/sampling.py`)**
Next, we will implement the advanced sampling. This begins with integrating the chosen SOAP descriptor library (e.g., DScribe) and adding it as a dependency in `pyproject.toml`. Then, we will implement the `FPSSampler` class. This class will contain a private helper method, `_calculate_soap_vectors`, which will efficiently compute the descriptors for all frames in a trajectory. The main `sample` method will then implement the core FPS algorithm, which is an iterative process: it first selects a random structure, then repeatedly selects the structure that is "farthest" (in the Euclidean distance of the SOAP feature space) from the set of already-selected structures. This component will be developed before the explorer is modified, allowing it to be unit-tested with pre-existing trajectory files.

**Step 3: Enhance Exploration Service (`services/exploration.py` and `shared/physics.py`)**
This is the most complex and critical step of the cycle. We will carefully refactor the `services/exploration.py` module.
*   **MLIP Integration:** We will add logic to the `_run_single_md` worker function to dynamically load a MACE model from a given file path. We will then implement a new ASE-compatible `MixedCalculator` class. This class will wrap both the MLIP calculator and a ZBL calculator, delegating calls to the MLIP by default, but smoothly switching to the ZBL potential for atom pairs that are closer than a certain cutoff distance.
*   **Hybrid MD/MC:** The main MD simulation loop within the worker function will be modified. At a specified frequency (e.g., every 10 steps), the loop will be paused to attempt an MC move. A random choice will be made between different move types (e.g., atom swap) based on their respective probabilities. The move will be provisionally applied, the change in potential energy will be calculated, and the move will be accepted or rejected based on the standard Metropolis criterion, which involves the energy change and the system temperature.
*   **Auto-Ensemble Switching:** We will implement the `detect_vacuum` function in `shared/physics.py`. This function will use a grid-based algorithm: it will discretize the simulation cell into a 3D grid and check if any grid points are sufficiently far from all atoms, which robustly identifies a vacuum slab. Before starting the MD run, the `_run_single_md` worker will call this function. Based on the boolean result, it will intelligently select either an `NPT` (for bulk systems) or `NVT` (for systems with a vacuum) ASE dynamics object for the simulation.

**Step 4: Update Configuration and CLI (`config/models.py`, `cli.py`)**
With the backend logic in place, we will update the Pydantic models in `config/models.py` to include all the new parameters defined in the Design Architecture. We will also update the `cli.py` file to add a new `ui` command (`mlip-autopipec ui`) that launches the Streamlit application.

**Step 5: Develop Web UI (`web_ui.py`)**
The final step is to build the user interface in `web_ui.py` using Streamlit.
*   **Configuration Form:** The UI will be structured with a sidebar containing a series of widgets (sliders for temperature, text inputs for elements, dropdowns for generator type) that directly correspond to the parameters in the `FullConfig` Pydantic model. We will use Streamlit's session state to manage the configuration as the user interacts with the form.
*   **Pipeline Execution:** A prominent "Run Pipeline" button will be the main call to action. When clicked, this will trigger a function that instantiates the `PipelineRunner` and calls its `execute` method with the configuration object built from the UI's state. To prevent the UI from becoming unresponsive during the run, the pipeline will be launched in a separate thread or subprocess, and logs will be streamed back to a text area in the UI in real-time.
*   **Output Visualization:** After a run completes, the UI will provide a results tab. This tab will connect to the output database, display a summary of the results, and show a paginated table of the generated structures. When a user clicks on a row in the table, it will use a library like `nglview` to render an interactive 3D visualization of the selected atomic structure.

## 5. Test Strategy

The testing strategy for Cycle 2 is necessarily more complex, designed to cover the sophisticated new scientific features and the new interactive user interface. We will expand our use of mocking for external dependencies (like MLIP models) and introduce a new class of end-to-end tests for the UI.

**Unit Testing Approach:**
Unit tests for Cycle 2 will be highly targeted, aiming to validate the correctness of the new, complex algorithms in isolation.
For the **Advanced Generators**, the `IonicGenerator` will be tested with various chemical formulas (e.g., CeO2, NaCl) to ensure that all generated structures are perfectly charge-neutral. The `InterfaceGenerator` tests will verify that the created structures have the correct slab orientations, interface separation distances, and produce no atomic overlaps at the interface.
For the **FPSSampler**, the key test will involve creating a small, artificial dataset of structures where some are nearly identical and one is distinctly different. We will assert that the `FPSSampler` always picks the dissimilar structure first, providing a deterministic validation of the core logic of the diversity selection algorithm. The calculation of the SOAP descriptors themselves will be mocked to isolate the FPS algorithm for this unit test.
The **Exploration Service** enhancements require the most careful and detailed testing. We will test the hybrid MD/MC logic on a simple 2-atom system (e.g., a Cu and an Au atom). The test will attempt a swap move, and we will assert that the atomic numbers of the two atoms are correctly exchanged. We will pre-calculate the energy change for this swap and test both the acceptance case (for a favorable energy change) and the rejection case (for an unfavorable one) of the Metropolis criterion. For MLIP integration, we will create a mock MACE calculator object that returns predictable, hardcoded energies and forces for specific configurations. We will then assert that the `Explorer` correctly calls this mock calculator and uses its output. This allows us to test the integration logic without needing a real, computationally expensive MLIP model. The auto-ensemble switching will be tested by creating a known bulk structure (e.g., a filled conventional cell) and a surface slab structure, and asserting that the `detect_vacuum` function in `physics.py` correctly returns `False` for the former and `True` for the latter, leading to the selection of NPT and NVT dynamics, respectively.

**Integration Testing Approach:**
The integration test suite will be expanded to validate the new end-to-end workflows that are now possible. A new, comprehensive test case will be created that runs the pipeline with a configuration that enables the hybrid MD/MC explorer and the `FPSSampler`. This test will use a mock MLIP model (as developed for the unit tests) to keep the runtime manageable within a CI environment. The assertions for this test will be scientifically focused. For example, in a CuAu alloy run with atom swaps enabled, we will connect to the generated database, calculate a chemical short-range order parameter for every structure, and assert that the distribution of this parameter shows a wider range of values compared to a run with pure MD. This provides quantitative proof that the swap moves are effectively exploring the configurational space. Another part of the assertion will involve calculating the average diversity of the dataset (e.g., the mean of the lower triangle of the SOAP kernel distance matrix) and confirming that it is significantly higher for the FPS-sampled run compared to a random-sampled run. This provides concrete evidence that the advanced sampling method is delivering its intended value. A separate, completely new class of integration test will be created for the Web UI. We will use the Playwright framework to write a script that automates interaction with the Streamlit application. The test script will perform actions like a real user: navigate to the UI's URL, interact with widgets to fill in the configuration form (e.g., `page.locator(...).select_option(...)`), click the "Run" button, and then wait for the run to complete by polling for a "success" message in the UI's log panel. Finally, it will check the UI's output display to verify that the results table is populated and that a 3D visualization widget is present. This end-to-end test validates that the entire UI is correctly wired to the backend pipeline.
