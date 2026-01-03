# CYCLE02 SPECIFICATION: Advanced Exploration and User Interface

## 1. Summary

This document provides the detailed technical specification for the second development cycle of the MLIP-AutoPipe project. Building upon the robust foundation established in Cycle 1, this cycle introduces a suite of advanced features designed to significantly enhance the scientific capabilities, versatility, and usability of the tool. The primary focus is on upgrading the exploration engine to a sophisticated hybrid Molecular Dynamics/Monte Carlo (MD/MC) model, which is essential for efficiently exploring the complex potential energy surfaces of many materials.

The key deliverables for Cycle 2 are:
1.  **Advanced Exploration:** The simple MD engine from Cycle 1 will be replaced with a powerful hybrid engine capable of performing MC atom swaps and other moves. This will be complemented by logic for automatic ensemble switching (NPT for bulk, NVT for slabs) and the integration of the ZBL potential to handle high-energy atomic interactions gracefully.
2.  **Expanded Material Support:** The tool's capabilities will be extended beyond simple alloys. New structure generators for ionic crystals and surface adsorption systems will be implemented.
3.  **Intelligent Sampling:** The basic random sampler will be augmented with a state-of-the-art Farthest Point Sampling (FPS) algorithm, enabling the selection of maximally diverse and informative structures for the final dataset.
4.  **Web User Interface:** A graphical user interface will be developed to provide an intuitive, interactive way to configure and run the pipeline, as well as to visualize the generated structures.

By the end of this cycle, MLIP-AutoPipe will transform from a foundational tool into a powerful and flexible framework suitable for advanced materials science research. It will be capable of generating high-quality datasets for a much wider range of physical systems and will be accessible to a broader audience through its interactive web interface.

## 2. System Architecture

The architecture in Cycle 2 expands upon the existing structure, adding new service implementations and a new top-level directory for the web-based GUI. The core orchestration logic and interfaces defined in Cycle 1 will remain largely unchanged, demonstrating the extensibility of the initial design.

```
src/
└── mlip_autopipec/
    ├── __init__.py
    ├── cli/
    │   └── ... (no changes)
    ├── core/
    │   └── ... (no changes)
    ├── domain/
    │   └── configuration.py     # **Add Pydantic models for new features**
    ├── services/
    │   ├── __init__.py
    │   ├── generation/
    │   │   ├── __init__.py
    │   │   ├── alloy.py
    │   │   ├── base.py
    │   │   ├── ionic.py         # **IonicGenerator implementation**
    │   │   └── adsorption.py    # **AdsorptionGenerator implementation**
    │   ├── exploration/
    │   │   ├── __init__.py
    │   │   └── md_engine.py     # **Upgrade to Hybrid MD/MC engine**
    │   ├── sampling/
    │   │   ├── __init__.py
    │   │   ├── fps.py           # **FarthestPointSampler implementation**
    │   │   └── random.py
    │   └── storage/
    │       └── ... (no changes)
    ├── utils/
    │   ├── __init__.py
    │   └── physics.py           # **Add vacuum detection logic**
    └── gui/                     # **New directory for Web UI**
        ├── __init__.py
        └── main_gui.py          # **Web UI application entrypoint (e.g., Streamlit)**

tests/
├── __init__.py
├── integration/
│   └── ... (add new integration tests for advanced features)
└── unit/
    ├── __init__.py
    ├── services/
    │   ├── generation/
    │   │   ├── test_ionic_generator.py # **Tests for IonicGenerator**
    │   │   └── test_adsorption_generator.py # **Tests for AdsorptionGenerator**
    │   └── sampling/
    │       └── test_fps.py               # **Tests for FarthestPointSampler**
    └── utils/
        └── test_physics.py         # **Tests for vacuum detection**
```

The factory classes in `core/factories.py` will be updated to recognize and instantiate the new implementations (e.g., `IonicGenerator`, `FarthestPointSampler`) based on the user's configuration.

## 3. Design Architecture

The design in Cycle 2 continues the schema-first approach, extending the Pydantic configuration models to support the new features. This ensures that the new, complex functionality is still governed by strict, validated data contracts.

**Extended Pydantic Configuration (`domain/configuration.py`):**
The existing configuration models will be extended using unions and discriminators to allow the user to select between different implementations.

```python
# Example of extending Pydantic models
from pydantic import BaseModel, Field, conlist
from typing import Union, Literal

# --- Exploration ---
class HybridMDMCConfig(BaseModel):
    # Inherits basic MD settings...
    mc_swap_probability: float = Field(0.1, ge=0, le=1, description="Probability of attempting an MC atom swap move")
    auto_ensemble_switching: bool = Field(True, description="Enable automatic NVT/NPT ensemble switching")

class ExplorerConfig(BaseModel):
    # Use a discriminated union to select the engine
    engine: Union[MDConfig, HybridMDMCConfig]

# --- Generation ---
class IonicGeneratorConfig(BaseModel):
    type: Literal["ionic"]
    cation: str
    anion: str
    # ... other ionic crystal parameters

class AdsorptionGeneratorConfig(BaseModel):
    type: Literal["adsorption"]
    slab_elements: conlist(str, min_length=1)
    adsorbate_molecule: str
    # ... other surface parameters

class GeneratorConfig(BaseModel):
    # A discriminated union allows selecting the generator type
    generator_type: Union[AlloyGeneratorConfig, IonicGeneratorConfig, AdsorptionGeneratorConfig] = Field(..., discriminator='type')

# --- Sampling ---
class FPSSamplerConfig(BaseModel):
    type: Literal["fps"]
    num_samples: int = Field(100, gt=0)

class RandomSamplerConfig(BaseModel):
    type: Literal["random"]
    num_samples: int = Field(100, gt=0)

class SamplerConfig(BaseModel):
    sampler_type: Union[FPSSamplerConfig, RandomSamplerConfig] = Field(..., discriminator='type')

```
This design allows for a highly flexible and self-documenting configuration. A user can simply change the `type` field in their YAML file to switch from a `random` sampler to an `fps` sampler, and Pydantic will automatically validate the corresponding block of parameters.

**Key Design Elements:**
-   **Hybrid MD/MC Engine:** The `MDEngine` service will be refactored to check its configuration. If MC moves are enabled, it will, at certain intervals during the MD run, pause the dynamics and attempt a Monte Carlo move (e.g., swapping the positions of two different atom types). This requires careful management of the simulation state.
-   **Vacuum Detection:** The `utils/physics.py` module will contain a new function, `detect_vacuum`. This function will analyze an `Atoms` object by creating a 3D grid over the simulation cell, identifying grid points that are far from any atom, and checking if these "empty" points form a contiguous layer. The `MDEngine` will call this function at the start of a run to decide whether to use an NVT (if vacuum is detected) or NPT (if not) ensemble.
-   **Web UI (`gui/main_gui.py`):** A simple web application will be built using a framework like Streamlit. It will dynamically generate input forms (sliders, text boxes) based on the Pydantic configuration models. When the user clicks "Run", the GUI will serialize the inputs into the same YAML format the CLI uses, and then it will call the `PipelineOrchestrator` using a `subprocess` to launch the calculation in the background. It will also provide a simple interface for visualizing the final structures from the output database.

## 4. Implementation Approach

The implementation of Cycle 2 features will be done in parallel where possible, focusing on maintaining the modularity of the system.

1.  **Extend Configuration:** First, update all the Pydantic models in `domain/configuration.py` to include the new options for hybrid MD/MC, new generators, and FPS sampling, as designed above.
2.  **Implement New Generators:** Create `ionic.py` and `adsorption.py` in `services/generation/`. Each will contain a new class that implements the `IStructureGenerator` interface. The `IonicGenerator` will use crystallographic logic to place atoms while ensuring charge neutrality. The `AdsorptionGenerator` will first build a surface slab and then place molecules on it.
3.  **Implement FPS Sampler:** Create `fps.py` in `services/sampling/`. The `FarthestPointSampler` class will implement the `ISampler` interface. It will use a library like `scikit-learn` or a custom implementation to compute SOAP descriptors for all structures in a trajectory and then iteratively select the structure that is farthest from the already-selected set.
4.  **Upgrade the Exploration Engine:** This is the most complex step. Refactor `services/exploration/md_engine.py`.
    -   Add the `detect_vacuum` utility in `utils/physics.py`.
    -   In the engine's main loop, add a call to this utility to select the correct ASE dynamics ensemble (NPT or NVT).
    -   Incorporate a probabilistic check inside the MD simulation loop to trigger MC swap moves.
    -   Add the logic to mix the MLIP calculator with ASE's ZBL potential for improved stability.
5.  **Update Factories:** Modify the factory classes in `core/factories.py` to recognize the new configuration options and return instances of the new generator and sampler classes.
6.  **Develop Web UI:** Create the `gui/main_gui.py` application.
    -   Build the UI components to reflect the Pydantic configuration structure.
    -   Implement the callback logic for the "Run" button to generate a config file and launch the pipeline as a background process.
    -   Add a results page to browse and visualize structures from an `ase.db` file.

## 5. Test Strategy

The testing strategy for Cycle 2 focuses on verifying the correctness of the new, scientifically complex components and ensuring they integrate smoothly into the existing pipeline.

**Unit Testing Approach (Min 300 words):**
New, highly specific unit tests will be developed for each new feature. For the `IonicGenerator`, unit tests will be created to assert that the generated structures are always perfectly charge-neutral. For example, when generating a structure for a material like NaCl, the test will assert that the number of Na atoms equals the number of Cl atoms. For the `AdsorptionGenerator`, tests will verify that the adsorbate molecule is placed at a reasonable distance from the slab surface and that the underlying slab structure remains intact.

The `FarthestPointSampler` will be tested with a simple, well-defined set of "structures" (e.g., 2D points) to verify that the algorithm correctly selects the points that are maximally far apart. This isolates the core logic of the FPS algorithm from the complexities of atomic structures and SOAP descriptors. The `detect_vacuum` function in `utils/physics.py` will be tested with several pre-defined `Atoms` objects: one representing a bulk crystal, one with a clear vacuum slab, and one with internal pores. The tests will assert that the function correctly returns `True` only for the slab case. These targeted unit tests are crucial for ensuring the correctness of the complex scientific logic being introduced in this cycle.

**Integration Testing Approach (Min 300 words):**
The end-to-end integration test suite will be expanded to cover the new workflows. A new test file, `test_advanced_pipelines.py`, will be added. One key test will configure and run the full pipeline for an ionic system (e.g., NaCl) using the new `IonicGenerator`. This test will validate that the entire system works with the new generator type.

Another critical integration test will focus on the hybrid MD/MC engine. It will run a short pipeline for a binary alloy (e.g., CuAu) with the hybrid engine enabled. After the run completes, the test will analyze the output trajectory file to programmatically verify that atom swap moves have actually occurred. This could be done by tracking the indices of specific atom types over time.

A third integration test will compare the output of the `RandomSampler` versus the `FarthestPointSampler`. It will run the pipeline twice on the same initial structure and exploration run, once with each sampler. The test will then calculate a diversity metric (e.g., the variance of energies or structures) on the two resulting databases and assert that the diversity score for the FPS-generated database is significantly higher than that of the randomly sampled one. These integration tests provide confidence that the new advanced features not only work correctly in isolation but also deliver their intended scientific benefits when plugged into the full pipeline.
