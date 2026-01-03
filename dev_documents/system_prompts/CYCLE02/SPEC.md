# Specification: CYCLE02 - Advanced Exploration and User Interface

## 1. Summary

This document provides the detailed technical specification for the second development cycle of the MLIP-AutoPipe project. CYCLE02 represents a significant leap forward, building directly upon the foundational pipeline established in CYCLE01. The primary objective of this cycle is to evolve the tool from a functional core into a sophisticated, intelligent, and versatile framework for MLIP dataset generation. This will be achieved by focusing on three key themes: implementing advanced, physics-aware exploration techniques; broadening the scope of supported material systems; and dramatically improving user accessibility through a graphical interface. Where CYCLE01 was about building a solid chassis and engine, CYCLE02 is about adding the high-performance transmission, the intelligent suspension, and the user-friendly dashboard.

The centerpiece of CYCLE02, and its most complex technical challenge, is the implementation of a **Hybrid MD/MC Exploration Engine**. This advanced engine will augment the standard Molecular Dynamics simulations from CYCLE01 with periodic Monte Carlo (MC) moves, specifically focusing on atomic swaps. This is a critical feature for efficiently exploring the vast and complex configuration space of materials like multi-component alloys. Pure MD can often get trapped in local energy minima, but the stochastic nature of MC swaps allows the simulation to overcome energy barriers and discover entirely new structural arrangements, leading to a far more diverse and comprehensive dataset. To complement this, we will also implement **Automatic Ensemble Switching**. This feature adds a layer of intelligence to the simulation setup, enabling the engine to automatically detect whether a system is a bulk material or a surface slab (by checking for vacuum regions) and then dynamically apply the correct thermodynamic ensemble—NPT for constant pressure simulations of bulk systems, and NVT for constant volume simulations of surfaces. This is a crucial feature that prevents unphysical simulation artifacts and removes a complex decision-making step from the user.

To further enhance the quality and efficiency of the generated datasets, this cycle will introduce **Farthest Point Sampling (FPS)** as a new, advanced sampling method. Unlike the naive random sampling from CYCLE01, which can inadvertently select many similar structures, FPS is an intelligent algorithm that uses structural descriptors to select a subset of configurations that are maximally different from each other. This ensures that the final dataset is information-rich and minimally redundant. The project's versatility will be expanded by adding a new generator for **Ionic Crystals**, which requires careful handling of constraints like charge neutrality. Finally, to make the tool more accessible to a wider range of users, particularly those less comfortable with the command line, a simple but effective **Web-based User Interface (UI)** will be developed. This GUI will provide an intuitive, interactive way to configure all pipeline parameters, launch runs, and visualize the resulting atomic structures, offering a powerful and user-friendly alternative to the CLI.

## 2. System Architecture

Development in CYCLE02 will seamlessly extend the modular architecture established in CYCLE01. The key to this smooth evolution is the interface-driven design. New services will be added, and existing ones will be modified to incorporate the advanced features, but the core `PipelineOrchestrator` will remain largely unchanged. This demonstrates the power and flexibility of the initial design, which was created with precisely this kind of extension in mind.

**File Structure for CYCLE02:**

The following tree shows the new and modified files for this cycle. The structure is a direct continuation of the work in CYCLE01. Bolded files are the primary focus of implementation.

```
src/mlip_autopipec/
├── cli/
│   └── **main.py**              # Add a new CLI command to launch the Web UI.
├── core/
│   ├── **factories.py**         # Update factories to recognize and dispatch new components (Ionic, FPS).
│   └── pipeline_orchestrator.py # No major changes needed, as it relies on abstract interfaces.
├── domain/
│   └── **configuration.py**     # Add and extend Pydantic models for all new components and features.
├── services/
│   ├── generation/
│   │   ├── __init__.py
│   │   ├── alloy.py
│   │   ├── base.py
│   │   └── **ionic.py**         # New concrete implementation for IonicGenerator.
│   ├── exploration/
│   │   ├── __init__.py
│   │   └── **md_engine.py**     # Heavily modify to become the HybridMDMC_Engine, adding new logic.
│   ├── sampling/
│   │   ├── __init__.py
│   │   ├── **fps.py**           # New concrete implementation for FarthestPointSampler.
│   │   └── random.py
│   └── storage/
│       ├── __init__.py
│       └── ase_db_wrapper.py
├── utils/
│   ├── __init__.py
│   └── **physics.py**           # Add new utility functions for vacuum detection.
└── gui/
    └── **main_gui.py**          # New Web UI application, likely using Streamlit.
```

**Component Blueprint:**

This section details the specific design and responsibilities of the new and modified components for CYCLE02.

*   **`services/exploration/md_engine.py`:** This file will undergo the most significant modification in this cycle. The existing `MDEngine` will be refactored and expanded into a more powerful and intelligent `HybridMDMC_Engine`. This is not a simple addition but a deep enhancement of its core logic.
    *   The simulation loop will be updated to periodically pause the MD integration and attempt a series of Monte Carlo moves, with the frequency and type of moves (e.g., which elements to swap) being controlled by the Pydantic configuration.
    *   At the start of each worker simulation, it will incorporate a call to the new `detect_vacuum` utility. The result of this check will be used in a conditional block to select and attach the correct ASE `Dynamics` object and associated thermostat/barostat—either `NPT` for bulk systems to allow cell relaxation or `Langevin` (NVT) for surface systems to maintain the cell volume. This automates a critical expert decision.
    *   It will also integrate a ZBL potential mixer. This involves using ASE's `MixedPotential` to combine the primary MLIP with a classical short-range repulsive potential (ZBL). This is a vital feature for robustness, as it prevents atoms from getting unphysically close at high simulation temperatures, a common cause of simulation crashes.

*   **`services/sampling/fps.py`:** This new file will contain the `FarthestPointSampler` class, which will be a new implementation of the `ISampler` interface.
    *   Its `sample` method will first read and aggregate all atomic configurations from the provided trajectory files.
    *   It will then use an external library (e.g., `dscribe`) to compute a high-dimensional structural descriptor for each configuration. The SOAP (Smooth Overlap of Atomic Positions) descriptor is a prime candidate.
    *   Finally, it will implement the iterative FPS algorithm. This starts by selecting a random structure, then iteratively finds the structure in the remaining set that is furthest away (in descriptor space) from the already selected set, and adds it to the selection. This is repeated until the desired number of samples is reached.

*   **`services/generation/ionic.py`:** This new file will contain the `IonicGenerator`.
    *   It will implement the `IGenerator` interface, making it a plug-and-play component.
    *   Its configuration will include the expected oxidation states for each element. The core logic of the generator will ensure that the final generated structure is perfectly charge-neutral, a fundamental physical constraint for ionic materials.
    *   It will leverage known crystallographic prototypes (like rock salt, cesium chloride, or perovskite structures) as templates, which it will then populate with the specified elements while respecting the charge neutrality constraint.

*   **`gui/main_gui.py`:** This new file will contain the entire Web UI application.
    *   It will be built using a simple yet powerful web framework like Streamlit or Flask, chosen for rapid development.
    *   The UI will dynamically generate input fields (sliders for temperature, text boxes for elements, dropdowns for lattice types) by introspecting the Pydantic configuration models. This creates a tight coupling between the configuration schema and the UI, ensuring they can never go out of sync.
    *   It will feature a prominent "Run Pipeline" button. When clicked, this will trigger a function that serializes the current state of the UI's input fields into a valid Hydra YAML configuration file on disk and then launches the main `mlip_autopipec run-pipeline` command in a detached subprocess.
    *   It will include a results section to display the status of the run (e.g., by tailing the log file) and, upon completion, to visualize some of the final generated structures using a library like `ase.visualize`.

*   **`domain/configuration.py`:** This central file will be updated to define the schemas for all the new features, ensuring they are also strongly typed, validated, and self-documenting. The `Union` type from Python's `typing` module will be used extensively to allow the user to select which component implementation they want to use in a type-safe manner (e.g., `sampler: Union[RandomSamplerConfig, FPSSamplerConfig]`).

## 3. Design Architecture

The design for CYCLE02 is a direct and logical extension of the schema-first, interface-driven architecture established in CYCLE01. The use of Pydantic schemas remains the cornerstone of the design, ensuring that the new, more complex configuration options are robust, easy to understand, and validated at the system's boundary.

**Pydantic Schema Design (`domain/configuration.py` additions):**

We will extend the existing set of Pydantic models to accommodate the new features. This is done by adding new models for the new components and using `Union` types to allow for a choice between implementations.

```python
# --- Additions and Modifications for CYCLE02 ---
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Literal, Union

# 1. New Ionic Generator Config
class IonicGeneratorConfig(BaseModel):
    elements: List[str] = Field(..., description="List of element symbols.")
    oxidation_states: Dict[str, float] = Field(..., description="Mapping of element to its oxidation state.")
    composition: Dict[str, int] = Field(..., description="Mapping of element to the number of atoms in the formula unit.")
    lattice_prototype: Literal["rocksalt", "perovskite", "cesiumchloride"] = Field(..., description="The crystal prototype to use.")

    @validator('composition')
    def cell_must_be_charge_neutral(cls, v, values):
        if 'oxidation_states' not in values:
            return v # Another validator will catch the missing field
        total_charge = sum(values['oxidation_states'][el] * num for el, num in v.items())
        if abs(total_charge) > 1e-6:
            raise ValueError("The formula unit must be charge-neutral.")
        return v

# Update GeneratorConfig to be a Union, discriminated by the 'name' field
class GeneratorConfig(BaseModel):
    name: Literal["alloy", "ionic"]
    params: Union[AlloyGeneratorConfig, IonicGeneratorConfig]

# 2. Updated Explorer Config for Hybrid MD/MC
class MonteCarloConfig(BaseModel):
    swap_elements: List[str] = Field(..., description="List of two element symbols to swap.")
    swap_frequency: int = Field(10, gt=0, description="Attempt a swap every N steps.")

class HybridMDMC_EngineConfig(MDEngineConfig): # Inherits from the old CYCLE01 config
    mc_moves: MonteCarloConfig
    use_auto_ensemble: bool = Field(True, description="Automatically switch between NVT/NPT based on vacuum detection.")
    use_zbl_repulsion: bool = Field(True, description="Use ZBL potential for short-range repulsion.")

# Update ExplorerConfig to be a Union
class ExplorerConfig(BaseModel):
    name: Literal["md", "hybrid_md_mc"]
    params: Union[MDEngineConfig, HybridMDMC_EngineConfig]
    max_workers: int = Field(4, gt=0)

# 3. New FPS Sampler Config
class FPSSamplerConfig(BaseModel):
    num_samples: int = Field(100, gt=0)
    descriptor: Literal["soap"] = Field("soap", description="The structural descriptor to use for FPS.")

# Update SamplerConfig to be a Union
class SamplerConfig(BaseModel):
    name: Literal["random", "fps"]
    params: Union[RandomSamplerConfig, FPSSamplerConfig]

# The FullConfig model is updated to use the new Union-based configs
class FullConfig(BaseModel):
    generator: GeneratorConfig
    explorer: ExplorerConfig
    sampler: SamplerConfig
    storage: StorageConfig
    class Config:
        extra = "forbid"
```

**Key Invariants and Constraints:**
*   **Charge Neutrality:** The new `IonicGeneratorConfig` will have a custom Pydantic validator that cross-references the `composition` and `oxidation_states` dictionaries to enforce the fundamental physical constraint that the defined chemical formula unit must be charge-neutral. This prevents a user from attempting to generate a physically impossible material.
*   **Discriminated Unions:** The use of `Union` types, critically combined with a `Literal` `name` field, creates a "discriminated union." This is a powerful pattern that allows Pydantic (and the factories) to unambiguously determine which schema to use for validation and which component to instantiate. It makes the configuration both flexible and type-safe.

**Extensibility:**
This design beautifully demonstrates the power of the Open/Closed Principle (open for extension, closed for modification), which was a core goal of the architecture established in CYCLE01. To add these major new features, we are primarily *adding* new code: new service classes (`ionic.py`, `fps.py`), new Pydantic models, and new logic to the `md_engine.py`. We are updating the factories to make them aware of the new components. Crucially, the `PipelineOrchestrator`, the very core of the application, requires almost no changes. It continues to operate on the abstract `IExplorer` and `ISampler` interfaces, blissfully unaware of the complex new logic happening inside the concrete implementations. This is the hallmark of a well-designed, extensible system.

## 4. Implementation Approach

The implementation of CYCLE02 will be an incremental process of integrating these new, advanced components into the existing, stable pipeline structure. The modular architecture allows us to develop and test each new feature in a relatively isolated manner before integrating it into the whole.

1.  **Schema Extension:** The first and most important step is to update the Pydantic schemas in `domain/configuration.py`. This involves adding the new models (`IonicGeneratorConfig`, `HybridMDMC_EngineConfig`, `FPSSamplerConfig`) and refactoring the top-level configuration models (`GeneratorConfig`, etc.) to use the `Union` type. This defines the "shape" of all the new features and provides the data contracts for the rest of the development.
2.  **Core Utilities First:** We will then implement the new, self-contained utility functions. The `detect_vacuum` logic will be added to `utils/physics.py`. This is a critical dependency for the new explorer engine, and it can be developed and unit-tested in complete isolation from the rest of the system.
3.  **Refactor and Enhance the Explorer:** This is the most significant implementation task of the cycle. We will heavily refactor the existing `MDEngine` in `services/exploration/md_engine.py` to become the `HybridMDMC_Engine`. This will involve:
    *   Modifying the main simulation loop to include a conditional block that, every `swap_frequency` steps, pauses the MD and executes a series of Monte Carlo swap moves.
    *   Adding logic at the very beginning of the worker function to call the newly implemented `detect_vacuum` utility and then select the appropriate ASE `Dynamics` object (`Langevin` for NVT or `NPT` for NPT).
    *   Implementing the ZBL potential mixing logic using ASE's `MixedPotential` class.
4.  **Implement New Services:** We will then implement the new, self-contained services.
    *   **`IonicGenerator`:** Create the new generator in `services/generation/ionic.py`. This will involve writing logic to construct a charge-neutral lattice based on the chosen crystallographic prototype.
    *   **`FarthestPointSampler`:** Create the new sampler in `services/sampling/fps.py`. This will require adding a new dependency (e.g., `dscribe`) to the `pyproject.toml` for calculating the SOAP descriptors. The core of the implementation will be the iterative FPS selection algorithm itself.
5.  **Update the Factories:** With the new services implemented, we will modify the factory classes in `core/factories.py`. The logic will be updated to handle the discriminated unions in the configuration. For example, the `SamplerFactory` will now inspect `config.sampler.name` and return an instance of `FarthestPointSampler` if the value is `"fps"`, or `RandomSampler` if it is `"random"`.
6.  **Develop the Web UI:** This is a major user-facing feature. We will create the `gui/main_gui.py` application using Streamlit. The development will focus on two parts:
    *   Building the UI components (widgets) to allow the user to set all the configuration parameters defined in our Pydantic models.
    *   Writing the backend logic for the "Run Pipeline" button, which will involve taking the current state of the UI widgets, programmatically generating a temporary Hydra YAML file, and then using Python's `subprocess` module to launch the main CLI command in a non-blocking way.
7.  **Update the CLI:** We will add a new, simple command to `cli/main.py`, such as `mlip_autopipec gui`, whose sole purpose is to launch the Streamlit application.
8.  **Comprehensive Testing:** Throughout this entire process, we will write a comprehensive suite of new unit and integration tests for all the new features to ensure they are not only correct in isolation but also function as expected within the full pipeline, and that they do not introduce any regressions in the existing CYCLE01 functionality.

## 5. Test Strategy

The testing strategy for CYCLE02 is designed to be comprehensive, ensuring that the new, complex algorithms are not only functionally correct but also effective in achieving their scientific goals. The strategy expands upon the foundation of tests from CYCLE01, adding new unit tests for the new logic and new integration tests for the new end-to-end workflows.

**Unit Testing Approach (Min 300 words):**
Unit tests are absolutely crucial in this cycle for verifying the correctness of the complex, algorithm-heavy new components in isolation. This allows us to debug and validate the core logic without the overhead of running a full pipeline.

*   **`HybridMDMC_Engine`:** We will write a suite of highly targeted unit tests for the new functionalities. The automatic ensemble switching logic will be tested by creating mock ASE `Atoms` objects (one representing a tightly-packed bulk material, another with a large vacuum slab) and passing them to the internal switching function. We will then assert that this function returns the correct ASE `Dynamics` class (`NPT` for bulk, `Langevin` for the slab). The MC swap move logic will be meticulously tested. A test will create a simple structure, count the number of atoms of each element, run the swap function, and then assert that the total composition of the cell remains identical, but that the specific identities of two atoms have been exchanged.

*   **`IonicGenerator`:** The generator will be tested to verify that the structures it produces adhere to fundamental physical laws. A key test will be to iterate through all supported `lattice_prototype` options, generate a structure for each, and assert that every single generated structure is perfectly charge-neutral based on the provided oxidation states. We will also write tests to ensure that for a given prototype like "rocksalt," the generated structure has the correct coordination number for each atom.

*   **`FarthestPointSampler`:** This component's selection algorithm is critical to test for correctness. A purely mathematical unit test will be created. Instead of using real atomic structures and slow descriptors, we will create a small, dummy set of 2D vectors (e.g., `[[0,0], [1,0], [10,0], [10,1]]`) where the most diverse subset is known by inspection. The test will run the FPS algorithm on this dummy data and assert that it correctly identifies and returns the known optimal subset (e.g., `[[0,0], [10,1]]`). This verifies that the core selection algorithm is mathematically sound without any external dependencies.

*   **`detect_vacuum` Utility:** The vacuum detection function in `utils/physics.py` will be tested with a variety of edge cases to ensure its robustness. We will test it with a tightly packed bulk structure (expect `False`), a cell with a clear vacuum slab along one axis (expect `True`), an completely empty cell (expect `True`), and a cell containing a molecule with a spherical void in the middle (expect `False`, as it's not a slab).

**Integration Testing Approach (Min 300 words):**
The end-to-end integration tests will be expanded to include new, sophisticated scenarios that specifically trigger and validate the advanced features introduced in this cycle. These tests are designed to ensure that the new components work correctly within the full pipeline and deliver their intended benefits.

*   **Hybrid MD/MC End-to-End Test:** A new integration test case will be added that uses a configuration file specifically designed to activate the hybrid engine with a high MC swap frequency. The test will run the full pipeline on a suitable binary alloy. After completion, the assertions will go beyond simply checking the final database. The test will programmatically load the intermediate trajectory files and perform an analysis to find direct evidence of successful swaps. For example, it could track the initial neighbors of a specific atom (e.g., a Germanium atom surrounded by Silicon) and then scan the trajectory to assert that at a later timestep, its list of neighbors has changed to include other Germanium atoms, an event that is statistically impossible in a short, low-temperature MD run but is expected with MC swaps.

*   **FPS vs. Random Sampling Comparative Test:** This is a crucial test to provide quantitative validation of the FPS sampler's effectiveness. The test will first run the pipeline up to the end of the exploration stage to generate a single, common set of trajectories. Then, it will execute only the sampling and storage stages twice. The first run will use the `random` sampler. The second run will use the `fps` sampler on the exact same trajectory files. The test will then load the two resulting databases. The primary assertion will be a quantitative comparison of the diversity of the two datasets. We will calculate a structural descriptor (like SOAP) for every structure in each database and then compute a diversity metric, such as the standard deviation of the pairwise descriptor distances. The test will assert that the diversity metric for the FPS-generated dataset is significantly higher than for the one generated by random sampling, providing clear proof of the feature's value.

*   **Ionic System End-to-End Test:** A new end-to-end test will be created that uses the `IonicGenerator` by setting `generator.name` to `"ionic"` in the configuration. It will run the full pipeline for a simple ionic compound (e.g., NaCl). The assertions will check that the final database contains structures that are all charge-neutral and have the expected elemental composition. This validates the new generator and ensures that the rest of the pipeline (explorer, sampler, storage) can correctly handle this different class of materials.
