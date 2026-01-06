# Specification: CYCLE 02 - Advanced Exploration, Sampling, and UI

## 1. Summary

This document provides the detailed specifications for the second development cycle of the MLIP-AutoPipe project. Building upon the robust foundation established in Cycle 01, this cycle is focused on introducing the advanced, scientifically sophisticated features that distinguish MLIP-AutoPipe as a powerful and intelligent data generation framework. The primary goals of this cycle are to enhance the exploration capabilities of the simulation engine, implement intelligent data sampling techniques to improve the quality of the final dataset, and provide a user-friendly web-based interface for interactive configuration and visualization.

The core of this cycle is the enhancement of the exploration engine. We will move beyond the basic NVT molecular dynamics of Cycle 01 and implement a more powerful hybrid MD/MC (Monte Carlo) engine. This will include the implementation of MC moves such as atom swaps, which are critical for efficiently exploring the compositional space in alloys. Furthermore, the engine will be upgraded with an automatic ensemble switching capability. By programmatically detecting the presence of a vacuum layer in a simulation cell, the engine will intelligently switch between NVT (for slab/surface systems) and NPT (for bulk systems) ensembles, a feature that significantly improves the physical realism of the simulations and removes a complex decision-making step from the user.

A key deliverable of this cycle is the implementation of a dedicated sampling stage in the pipeline. While Cycle 01 stored all generated trajectory data, this is often inefficient and redundant. Cycle 02 will introduce a `SamplingEngine` that can intelligently select a small, diverse subset of structures from the vast number of configurations generated during exploration. The flagship feature of this engine will be the implementation of Farthest Point Sampling (FPS), a powerful technique that uses structural descriptors (SOAP) to maximize the diversity of the selected data points. This ensures that the final dataset is rich in information and highly efficient for training MLIPs. Finally, to improve usability and accessibility, this cycle includes the development of a simple, intuitive Web UI. This interface will provide a graphical way for users to create and validate configuration files and visualize the atomic structures, making the tool more approachable for a broader audience.

## 2. System Architecture

The architecture in Cycle 02 is an extension of the foundation built in Cycle 01. The existing modular structure allows for the seamless integration of new components, such as the `SamplingEngine`, and the enhancement of existing ones, like the `ExplorerEngine`, without requiring major refactoring.

**File Structure (Cycle 02):**

The file tree below highlights the new and modified files for this cycle. Files in **bold** are the primary targets for creation or significant modification.

```
.
└── src/
    └── mlip_autopipec/
        ├── **cli.py**                # Modified to include Web UI command
        ├── config/
        │   └── **models.py**         # Extended with new config options
        ├── generators/
        │   ├── base.py
        │   ├── alloy.py
        │   └── ionic.py
        ├── explorers/
        │   └── **engine.py**         # Enhanced with hybrid MD/MC and auto-ensemble
        ├── sampling/
        │   ├── __init__.py
        │   ├── **base.py**           # Abstract BaseSamplingEngine class
        │   ├── **fps.py**            # Farthest Point Sampling implementation
        │   └── **random.py**         # Simple random sampler
        ├── storage/
        │   └── database.py
        ├── **pipeline.py**           # Modified to include the sampling stage
        └── ui/
            ├── __init__.py
            └── **app.py**            # Web UI application (Streamlit/Flask)
```

**Code Blueprints:**

This section provides the detailed blueprints for the key modifications and new components in Cycle 02.

**`src/mlip_autopipec/explorers/engine.py` (Enhancement)**
The `ExplorerEngine` will be significantly upgraded. The core change will be in the worker function that runs the simulation, which will now incorporate the hybrid MD/MC logic and other advanced features.

```python
# (Inside the MD worker function)

# ... setup calculator and dynamics ...

def run_single_md_process(atoms: Atoms, config: ExplorationConfig) -> List[Atoms]:
    # ...

    # 1. Automatic Ensemble Switching
    if detect_vacuum(atoms): # detect_vacuum is a new utility function
        # Use NVT for slabs/surfaces
        dyn = Langevin(atoms, timestep=config.time_step_fs * units.fs, ...)
    else:
        # Use NPT for bulk materials
        dyn = NPT(atoms, timestep=config.time_step_fs * units.fs, ...)

    # 2. Hybrid MD/MC Logic
    trajectory = []
    for step in range(config.num_steps):
        dyn.run(1) # Run one MD step

        if step % config.mc_interval == 0:
            # Attempt a Monte Carlo move
            if config.enable_atom_swap:
                # new_atoms = attempt_atom_swap(atoms, temperature_k)
                # if accepted, atoms = new_atoms
                pass

        trajectory.append(atoms.copy())

    return trajectory

# Utility function to be implemented (e.g., in a new `physics.py` module)
def detect_vacuum(atoms: Atoms) -> bool:
    # Implementation using grid-based detection of empty space
    pass
```

**`src/mlip_autopipec/sampling/base.py` (New File)**
This new file defines the abstract interface for all sampling strategies, ensuring they can be used interchangeably by the pipeline.

```python
from abc import ABC, abstractmethod
from typing import List
from ase import Atoms

class BaseSamplingEngine(ABC):
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def sample(self, trajectories: List[List[Atoms]]) -> List[Atoms]:
        """
        Takes a list of trajectories (list of lists of Atoms) and returns
        a single list of down-selected Atoms objects.
        """
        pass
```

**`src/mlip_autopipec/sampling/fps.py` (New File)**
This file will contain the implementation of the Farthest Point Sampling algorithm. It will require dependencies like `dscribe` for calculating the SOAP descriptors.

```python
from typing import List
from ase import Atoms
# from dscribe.descriptors import SOAP
# from sklearn_extra.cluster import KMedoids (or a manual implementation of FPS)

from .base import BaseSamplingEngine

class FPSSamplingEngine(BaseSamplingEngine):
    def sample(self, trajectories: List[List[Atoms]]) -> List[Atoms]:
        # 1. Flatten all trajectories into a single list of structures
        all_structures = [frame for traj in trajectories for frame in traj]

        # 2. Calculate SOAP descriptors for all structures
        # soap_desc = self._calculate_soap(all_structures)

        # 3. Implement FPS algorithm to select N diverse structures
        # selected_indices = self._run_fps(soap_desc, self.config.num_samples)

        # 4. Return the selected structures
        # return [all_structures[i] for i in selected_indices]
        pass # Placeholder for the complex implementation
```

**`src/mlip_autopipec/pipeline.py` (Enhancement)**
The `PipelineRunner` will be updated to include the new sampling stage, which will be executed after exploration and before storage.

```python
# (Inside PipelineRunner class)

    def _get_sampler(self) -> BaseSamplingEngine:
        # Factory logic for samplers
        if self.config.sampling.method == "fps":
            from .sampling.fps import FPSSamplingEngine
            return FPSSamplingEngine(self.config.sampling)
        elif self.config.sampling.method == "random":
            from .sampling.random import RandomSamplingEngine
            return RandomSamplingEngine(self.config.sampling)
        else:
            raise NotImplementedError(...)

    def run(self):
        # Stage 1: Generation (as before)
        # ...

        # Stage 2: Exploration (as before)
        # ...
        explored_trajectories = explorer.run_parallel(initial_structures)

        # Stage 3: Sampling (New Stage)
        sampler = self._get_sampler()
        sampled_structures = sampler.sample(explored_trajectories)

        # Stage 4: Storage (Modified to store sampled structures)
        self.db_wrapper.write_atoms_list(sampled_structures, "final_dataset")
```

This extended architecture seamlessly integrates the new, advanced functionalities while maintaining the modular and decoupled design established in the first cycle. The introduction of an abstract base class for sampling engines mirrors the successful pattern used for generators, ensuring the system remains extensible and easy to maintain.

## 3. Design Architecture

The design for Cycle 02 focuses on enriching the application's domain logic and providing a more interactive user experience. The Pydantic configuration schema will be extended to control the new features, and a new service layer for the Web UI will be introduced.

**Pydantic-based Schema Design (Extensions):**

The existing configuration models in `src/mlip_autopipec/config/models.py` will be updated.

```python
from pydantic import BaseModel, Field
from typing import List, Dict, Literal

# ... Existing SystemConfig ...

class ExplorationConfig(BaseModel):
    # ... Existing fields ...
    enable_atom_swap: bool = Field(False, description="Enable Monte Carlo atom swaps.")
    mc_interval: int = Field(10, gt=0, description="Perform MC move every N MD steps.")
    # Add fields for ZBL potential, etc.

class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""
    method: Literal["fps", "random"] = Field("random", description="Sampling method to use.")
    num_samples: int = Field(100, gt=0, description="Number of structures to select for the final dataset.")
    # Add FPS-specific parameters if needed

class FullConfig(BaseModel):
    # ... Existing fields ...
    sampling: SamplingConfig # Add the new sampling configuration
```

**Domain Concepts and Invariants:**
*   **`ExplorationConfig`**: This model is extended to manage the more complex simulation logic.
    *   **Invariants**: The `mc_interval` must be a positive integer. If `enable_atom_swap` is true, the system must be multi-component.
    *   **Consumers**: The enhanced `ExplorerEngine` will use these new fields to control the hybrid MD/MC simulation loop.
*   **`SamplingConfig`**: This new model defines the parameters for the data selection stage.
    *   **Invariants**: `num_samples` must be positive and should ideally be less than the total number of frames generated, although this is a soft constraint.
    *   **Producers**: User-defined in the YAML configuration.
    *   **Consumers**: The `PipelineRunner` uses this to select the appropriate sampling engine, which in turn uses the parameters to perform the sampling.

**Web UI Design (`src/mlip_autopipec/ui/app.py`):**
The Web UI will be a simple, single-page application built with a framework like Streamlit for rapid development. Its purpose is twofold: configuration and visualization.
1.  **Configuration Builder**: The UI will present all the configuration options defined in the Pydantic models (`FullConfig`) using interactive widgets (sliders for numbers, text boxes for strings, dropdowns for choices). As the user modifies these widgets, the UI will dynamically build a corresponding YAML configuration in a text area. This provides a user-friendly way to create valid configuration files without needing to remember all the specific keys and data types.
2.  **Structure Visualizer**: The UI will include a file uploader that allows the user to upload an ASE-readable file (like `.xyz` or a database file). It will then use a library like `nglview` to render the 3D atomic structure, allowing for interactive inspection (rotating, zooming). This is invaluable for quickly verifying the output of the generation or exploration stages.
The UI will be a separate concern and will not contain any core pipeline logic. It will be launched via a new CLI command (e.g., `mlip-autopipec ui`).

## 4. Implementation Approach

The implementation of Cycle 02 will be phased to tackle the new features in a logical order, starting with the backend enhancements and concluding with the user-facing UI.

**Step 1: Extend Configuration Models**
The first step is to update the Pydantic models in `config/models.py` as detailed in the Design Architecture. The `ExplorationConfig` will be extended, and the new `SamplingConfig` will be added. The `FullConfig` model will be updated to include the new sampling section.

**Step 2: Enhance the Exploration Engine**
This is the most complex backend task. The `ExplorerEngine`'s worker function in `engine.py` will be refactored.
*   A new utility module, perhaps `src/mlip_autopipec/physics/validation.py`, will be created to house the `detect_vacuum` function.
*   The worker function will be modified to call `detect_vacuum` and select the appropriate ASE dynamics object (NVT or NPT).
*   The main simulation loop will be changed from a single `dyn.run()` call to a `for` loop that runs one step at a time, allowing for the periodic introduction of MC moves.
*   The logic for `attempt_atom_swap` will be implemented. This involves randomly selecting two atoms of different species and accepting or rejecting the swap based on the Metropolis criterion (change in energy).

**Step 3: Implement the Sampling Framework**
*   The abstract base class `BaseSamplingEngine` will be created in `sampling/base.py`.
*   A simple `RandomSamplingEngine` will be implemented first in `sampling/random.py`. This is straightforward and will allow for the integration of the sampling stage into the pipeline before the more complex FPS is ready.
*   The `FPSSamplingEngine` will then be implemented in `sampling/fps.py`. This is a significant sub-task that will require integrating `dscribe` to compute SOAP descriptors and then implementing the iterative FPS algorithm to select the most diverse structures.

**Step 4: Integrate Sampling into the Pipeline**
The `PipelineRunner` in `pipeline.py` will be updated. The factory method `_get_sampler` will be added. The main `run` method will be modified to include the third "Sampling" stage. It will call the selected sampler after the exploration stage and pass the resulting sampled structures to the storage stage.

**Step 5: Develop the Web UI**
A new directory `src/mlip_autopipec/ui` will be created. Inside `app.py`, a Streamlit application will be built.
*   The UI layout will be created with a sidebar for configuration options and a main panel for visualization and the generated YAML output.
*   Streamlit widgets will be mapped to the fields in the Pydantic `FullConfig` model.
*   An `nglview` component will be used for 3D visualization.
*   A new command, `ui`, will be added to `cli.py` to launch the Streamlit app.

**Step 6: Final Integration and Testing**
All new features will be tested together. An end-to-end integration test will be developed that uses a configuration enabling all the new features: hybrid MD/MC, FPS sampling, etc. The test will verify that the final database contains a smaller, sampled set of structures and that the simulation ran without errors.

## 5. Test Strategy

The testing for Cycle 02 will be layered on top of the existing test suite, with a strong focus on the new, complex functionalities.

**Unit Testing Approach (Min 300 words):**
The unit tests for Cycle 02 will be crucial for validating the complex scientific logic being introduced.
*   **Explorer Engine (`tests/unit/explorers/test_engine.py`):** The enhanced `ExplorerEngine` requires targeted unit tests. The `detect_vacuum` function will be tested with known bulk and slab `Atoms` objects to ensure it correctly classifies them. A separate test will focus on the hybrid MD/MC logic. We can create a mock `Atoms` object and a mock dynamics object and verify that the simulation loop correctly attempts a Monte Carlo swap at the specified `mc_interval`. The acceptance/rejection logic of the swap can also be tested by mocking the energy change.
*   **Sampling Engines (`tests/unit/sampling/test_fps.py`, `test_random.py`):** Testing the sampling engines is vital. For the `RandomSamplingEngine`, we can provide it with a dummy trajectory of 100 frames and assert that it correctly returns the `num_samples` specified in the configuration. Testing the `FPSSamplingEngine` is more involved. We will create a small, artificial dataset of structures where some are very similar and one is very different. The test will assert that the FPS algorithm correctly picks the unique structure first, thus validating the core of its diversity-maximization logic. We will mock the expensive SOAP descriptor calculation to keep the test fast.

**Integration Testing Approach (Min 300 words):**
The integration tests will ensure that the new components work correctly within the larger pipeline and that the data flows correctly between all four stages.
*   **Full End-to-End Pipeline (`tests/integration/test_full_pipeline_cycle02.py`):** The primary integration test for this cycle will be a comprehensive test of the entire, four-stage workflow. A specific configuration file will be created for this test, enabling all the new features: a hybrid MD/MC run on an alloy, followed by FPS sampling. The test will invoke the `mlip-autopipec run` command. After execution, the test will perform a series of assertions on the output database (`mlip_data.db`):
    1.  It will verify that the final dataset stored under the key "final_dataset" exists.
    2.  It will assert that the number of structures in this final dataset is exactly equal to the `num_samples` specified in the configuration, proving that the sampling stage was executed correctly.
    3.  It will check that the number of structures in the final dataset is significantly smaller than the total number of frames that would have been generated in the exploration trajectories, confirming that down-sampling occurred. This single test provides high confidence that all the new backend components are correctly integrated and working in sequence.
*   **Web UI Test (Manual or using Selenium/Playwright):** While harder to automate, a test plan for the UI will be created. This will involve launching the UI with the `mlip-autopipec ui` command and manually verifying that changing the widgets correctly updates the YAML output. It will also involve uploading a sample `.xyz` file to ensure the visualizer renders the structure correctly.
