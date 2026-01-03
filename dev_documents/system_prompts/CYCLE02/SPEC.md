# SPEC.md for CYCLE02: Automation, Active Learning & Advanced Simulation

## 1. Summary

This document provides the detailed technical specification for CYCLE02 of the MLIP-AutoPipe project. This cycle represents the culmination of the project's vision, building directly upon the foundational manual pipeline established in CYCLE01. The singular, ambitious goal of this cycle is to introduce and integrate the intelligent automation that is the core promise of the project. The primary objective is to fully "remove the human expert from the loop" by implementing a complete, end-to-end autonomous workflow that can generate a high-quality, validated Machine Learning Interatomic Potential (MLIP) from a bare-minimum user specification. This cycle transforms the tool from a manually-operated instrument into a fully autonomous, "fire-and-forget" scientific discovery platform.

To achieve this, CYCLE02 will deliver the remaining critical modules of the system architecture: the `StructureGenerator` (Module A) for intelligent and diverse initial data seeding without resorting to expensive AIMD; the `ExplorerSampler` (Module B) for computationally efficient exploration of the vast potential energy surface using a universal surrogate MACE model; and the `SimulationEngine` (Module E), which closes the loop by providing the on-the-fly (OTF) active learning capability. The `WorkflowOrchestrator`, which was a simple sequential runner in the previous cycle, will be significantly enhanced to become a sophisticated state machine, capable of managing the complex, cyclical, and asynchronous state transitions of the active learning loop. Furthermore, this cycle will introduce the crucial user-facing `ConfigExpander`, an expert system that translates a simple, high-level user request (e.g., the elements of a material) into the detailed, exhaustive `FullConfig` that the pipeline requires for execution. Upon the successful completion of CYCLE02, the MLIP-AutoPipe will be a feature-complete, autonomous system. A user will be able to provide a simple `input.yaml` file containing just the chemical elements and composition, and the system will autonomously execute the entire multi-stage process of data generation, DFT labeling, model training, and active refinement to produce a validated, high-quality potential, ready for production simulations.

## 2. System Architecture

CYCLE02 involves implementing the remaining automation-focused modules and significantly enhancing the orchestrator and CLI to support the full, autonomous workflow. The architecture, as defined in the main `SYSTEM_ARCHITECTURE.md`, will now be feature-complete, with all modules present and fully integrated.

**File Structure for CYCLE02:**
(Bold files/directories are to be newly created in this cycle. *Italicized* files are those from CYCLE01 that will be significantly modified and enhanced.)

```
.
├── pyproject.toml
├── README.md
├── src
│   └── mlip_autopipec
│       ├── __init__.py
│       ├── *cli.py*
│       ├── *config.py*
│       ├── *orchestrator.py*
│       ├── database.py
│       └── modules
│           ├── __init__.py
│           ├── **a_structure_generator.py**
│           ├── **b_explorer_sampler.py**
│           ├── *c_labeling_engine.py*
│           ├── d_training_engine.py
│           └── **e_simulation_engine.py**
└── tests
    ├── __init__.py
    ├── *test_config.py*
    ├── *test_orchestrator.py*
    └── modules
        ├── __init__.py
        ├── **test_a_structure_generator.py**
        ├── **test_b_explorer_sampler.py**
        ├── *test_c_labeling_engine.py*
        ├── test_d_training_engine.py
        └── **test_e_simulation_engine.py**
```

**Component Blueprint:**

*   **`config.py`**: This file will be substantially updated to support the new autonomous workflow. It will introduce the `MinimalConfig` Pantic model, which represents the simple, user-facing `input.yaml`. A new and critical component, the `ConfigExpander` class, will be implemented here. This class's core logic will be responsible for programmatically converting a `MinimalConfig` instance into a `FullConfig` instance. It will achieve this by applying a codified set of physical heuristics, such as looking up recommended DFT cutoff energies for the specified elements from a data file based on the SSSP protocol, or setting default temperature ranges for simulations based on tabulated melting points.

*   **`modules/a_structure_generator.py`**: This new file will contain the `StructureGenerator` class. It will implement a suite of different strategies for creating an initial, physically diverse set of atomic structures without relying on expensive first-principles calculations. It will include methods for generating Special Quasirandom Structures (SQS) for modeling disordered alloys, performing Normal Mode Sampling (NMS) to explore the vibrational space of molecules, and creating various rattled and strained structures to capture the elastic behavior of crystalline materials. The primary output of this module will be a list of `ase.Atoms` objects, which will be added to the database to seed the learning process.

*   **`modules/b_explorer_sampler.py`**: This new file will house the `ExplorerSampler` class, a component that is absolutely critical for the computational efficiency of the entire pipeline. It will use a pre-trained, universal machine learning potential (specifically, a MACE model from the Materials Platform) to run fast, large-scale molecular dynamics simulations. From the resulting multi-million-frame trajectory, it will use a sophisticated workflow of descriptor calculation (e.g., using SOAP to characterize local atomic environments), dimensionality reduction, and clustering to identify the unique structural motifs present. It will then use a stratified sampling algorithm to select a small but highly representative set of geometrically diverse structures to be passed to the computationally expensive DFT labeling engine.

*   **`modules/e_simulation_engine.py`**: This new file will contain the `SimulationEngine` class. This class takes the MLIP that was trained by Module D and uses it to run production-level simulations (e.g., molecular dynamics to study thermal properties or kinetic Monte Carlo to study rare events). Its most important feature is the implementation of on-the-fly (OTF) uncertainty quantification. It will monitor the predictions of the MLIP (or an ensemble of MLIPs) during the simulation. If the uncertainty for a given atomic configuration exceeds a dynamically adjusted threshold, the engine will pause, extract that "uncertain" structure, and flag it for re-training, thereby closing the active learning loop.

*   **`modules/c_labeling_engine.py`**: The existing `LabelingEngine` from CYCLE01 will be enhanced with the advanced, multi-stage error recovery logic that was specified in the original requirements. This will make it far more robust. For instance, upon an SCF convergence failure, it will now automatically try a sequence of recovery strategies, such as reducing the electronic mixing parameter (`mixing_beta`), switching to a more conservative mixing mode, or temporarily increasing the smearing temperature, before finally giving up on a structure.

*   **`orchestrator.py`**: The `WorkflowOrchestrator` will undergo a complete transformation. It will be refactored from a simple sequential runner into a sophisticated, stateful class that manages the entire asynchronous, cyclical logic of the active learning loop. It will have a main `run_autonomous()` method that executes the full workflow:
    1.  Invoke the `ConfigExpander`.
    2.  Call Module A to `run_initial_seeding()`.
    3.  Begin the main loop:
        a. Call Module B to `run_exploration_phase()`.
        b. Call Module C to `run_labeling_phase()` on all currently unlabeled structures.
        c. Call Module D to `run_training_phase()` to produce a new, improved MLIP.
        d. Call Module E to `run_simulation_phase()`.
        e. If uncertain structures are found, add them to the database and continue the loop.
        f. If the simulation completes without finding uncertain structures, or if the maximum number of generations is reached, terminate the workflow.

*   **`cli.py`**: The CLI will be greatly simplified from the user's perspective to reflect the new autonomous nature of the system. The `label` and `train` commands from the previous cycle will be either removed or hidden as internal debugging commands. The primary, and likely only, user-facing entry point will become `mlip-pipe run <input.yaml>`, which takes the minimal user configuration file and triggers the entire, unattended workflow managed by the enhanced orchestrator.

## 3. Design Architecture

The architecture in CYCLE02 evolves significantly to support a fully autonomous, intelligent, and closed-loop system. The Pydantic schemas are expanded to include the user-facing minimal configuration, and the interactions between the modules become more complex and dynamic, managed entirely by the state machine within the `WorkflowOrchestrator`.

**Pydantic Schema Design (`config.py`):**

The `config.py` file will now define two distinct top-level models, representing the two levels of user interaction with the system. The `ConfigExpander` will be the bridge between them.

```python
#
# Pseudocode for additions and modifications to config.py
#
from pydantic import BaseModel, Field, PositiveFloat, PositiveInt
from typing import List, Optional, Dict

# The NEW user-facing configuration model
class MinimalConfig(BaseModel):
    """A simple, high-level configuration provided by the user."""
    elements: List[str] = Field(..., min_items=1, description="A list of the chemical elements in the system.")
    composition: str = Field(..., description="A string representing the composition, e.g., 'FePt' or 'Si2O4'.")
    temperature_range_k: Optional[List[int]] = Field(default=[300, 1000], description="Target temperature range in Kelvin.")

# The ConfigExpander will be a class that performs this conversion.
class ConfigExpander:
    def __init__(self, minimal_config: MinimalConfig):
        self.minimal = minimal_config
        self.full = FullConfig(...) # Initialize a default FullConfig

    def expand(self) -> FullConfig:
        # Heuristics are applied here.
        self._expand_dft_settings()
        self._expand_simulation_settings()
        return self.full

    def _expand_dft_settings(self):
        # Example: lookup SSSP recommendations for elements to find ecutwfc
        # Example: if 'Fe' in elements, set magnetism to 'ferromagnetic'
        pass

    def _expand_simulation_settings(self):
        # Example: set MD temperature steps based on temperature_range_k
        pass

# FullConfig will be significantly expanded to include settings for the new modules.
class ActiveLearningSettings(BaseModel):
    """Settings controlling the active learning loop."""
    strategy: str = "direct_then_active"
    surrogate_model: str = "mace_mp"
    max_generations: PositiveInt = 10
    uncertainty_threshold_strategy: str = "dynamic_95percentile"

class StructureGenerationSettings(BaseModel):
    """Settings for initial structure generation."""
    sqs_size: List[int] = [4, 4, 4]
    rattle_stdev: PositiveFloat = 0.1

class FullConfig(BaseModel):
    """The complete, internal configuration, now with settings for all modules."""
    # ... existing fields from CYCLE01 (SystemSettings, DFTSettings, MLIPSettings) ...
    active_learning: ActiveLearningSettings
    structure_generation: StructureGenerationSettings
    # ... other new fields for modules B and E ...
```

**Key Architectural Invariants and Constraints:**
*   **Orchestrator as the Single Point of Truth for State**: No module will ever call another module directly. All interactions are mediated by the `WorkflowOrchestrator`. This is a critical invariant for managing the complex, non-linear state of the active learning loop. The orchestrator is the only component that knows, for instance, which "generation" of the loop is currently active.
*   **Idempotency of Data Labeling**: The labeling step must be designed to be idempotent. If an "unlabeled" structure is submitted to the `LabelingEngine` and the process is interrupted for any reason (e.g., a machine reboot), re-running the labeling step should either gracefully skip the already-completed calculation or correctly overwrite it without creating duplicate or inconsistent data in the database.
*   **Explicit Convergence Criteria**: The `WorkflowOrchestrator` must have a clear and explicit set of convergence criteria to gracefully terminate the active learning loop. This will be defined in the `ActiveLearningSettings` and will likely include a maximum number of generations (`max_generations`) as a hard stop, potentially supplemented by a more sophisticated criterion, such as the loop terminating when a long simulation runs to completion without triggering the uncertainty threshold.

**Data Consumers and Producers:**
*   **`ConfigExpander`**: Consumes a `MinimalConfig` object, produces a `FullConfig` object.
*   **`StructureGenerator` (Module A)**: Consumes the `FullConfig` object, produces a list of `ase.Atoms` objects that are written to the database.
*   **`ExplorerSampler` (Module B)**: Consumes the `FullConfig` object, produces a list of high-value `ase.Atoms` objects that are queued for DFT labeling.
*   **`SimulationEngine` (Module E)**: Consumes a trained MLIP model file and the `FullConfig`. Produces simulation trajectory data for the user and, most importantly, produces new "uncertain" `ase.Atoms` objects that are fed back into the labeling queue, forming the core of the active learning process.

## 4. Implementation Approach

The implementation for CYCLE02 will focus on building the new automation modules first as standalone components, then integrating them into the significantly enhanced orchestrator, and finally exposing the simplified workflow through the CLI.

1.  **Update `config.py`**: The first and most important step is to implement the `MinimalConfig` Pantic model and the logic for the `ConfigExpander`. This step defines the new, simplified user experience and provides the necessary configuration infrastructure for all the other new modules.
2.  **Implement `modules/a_structure_generator.py`**: We will build the `StructureGenerator` class. This may involve identifying and integrating external libraries like `pymatgen` or `icet` for robust SQS generation. We will write unit tests to verify that the generated structures have the correct properties (e.g., composition, cell size, diversity).
3.  **Implement `modules/b_explorer_sampler.py`**: This is a major implementation step.
    *   First, we will integrate the MACE model. This will require adding `mace-torch` or a similar library as a dependency in `pyproject.toml`.
    *   Next, we will implement a simple molecular dynamics loop using ASE's MD classes, driven by the MACE calculator.
    *   Then, we will implement the descriptor calculation pipeline, likely by integrating the `dscribe` library for SOAP descriptors, and a standard clustering algorithm, for instance `sklearn.cluster.KMeans`.
    *   **Testing**: Unit tests for this module are vital. We will create a test that runs a very short MD trajectory (perhaps only a few dozen steps) on a small, simple system and verifies that the sampling logic correctly selects at least one structure.
4.  **Implement `modules/e_simulation_engine.py`**:
    *   First, we will implement the main simulation loop (likely MD using ASE's `VelocityVerlet`).
    *   The core of this module is the uncertainty logic. This will involve implementing a mechanism to query the MLIP model (or an ensemble of models, if applicable) for the variance or standard deviation of its prediction.
    *   We will then implement the dynamic thresholding logic (e.g., by tracking the distribution of uncertainty values on the existing training set and flagging any new configurations that are significant outliers).
    *   **Testing**: We will mock the MLIP calculator object. We will write a test where the mock calculator is programmed to return a high uncertainty value after a predefined number of steps, and we will assert that the `SimulationEngine` correctly captures the corresponding atomic structure and exits gracefully.
5.  **Enhance `orchestrator.py`**: This is the main integration step of the entire project.
    *   We will refactor the `WorkflowOrchestrator` to manage the new state machine logic for the full active learning loop. This will likely involve breaking its `run` method into several private helper methods for each stage of the loop.
    *   We will then integrate the newly implemented modules (A, B, and E) into the orchestrator's workflow.
    *   The orchestrator's main `run_autonomous()` method will now be the top-level entry point that executes the entire unattended pipeline from start to finish.
6.  **Update `cli.py`**:
    *   We will refactor the main command to be `mlip-pipe run input.yaml`.
    *   This command will now perform the full, high-level sequence: instantiate the `MinimalConfig` from the YAML file, run the `ConfigExpander` to get the `FullConfig`, initialize the `WorkflowOrchestrator` with this full configuration, and finally, call the orchestrator's main `run_autonomous()` method.
7.  **Enhance `modules/c_labeling_engine.py`**: As a final step, we will add the advanced, multi-stage error recovery logic to the `LabelingEngine` to improve the overall robustness and reliability of the long-running, autonomous workflow.

## 5. Test Strategy

The testing strategy in CYCLE02 must necessarily be more sophisticated, focusing on the complex logic of the new automation modules and, most importantly, the intricate state transitions and feedback mechanisms of the active learning loop.

**Unit Testing Approach:**

Our unit testing approach will continue to follow the strict principle of isolation for each new component. For the `StructureGenerator` (Module A), we will write tests that call its public methods (e.g., `generate_sqs`) and assert that the returned list of `ase.Atoms` objects has the specified number of atoms, the correct chemical elements, and the correct overall stoichiometry. We will not be testing the underlying SQS algorithm itself (as we trust the third-party library), but rather that our wrapper correctly configures and uses that library.

For the `ExplorerSampler` (Module B), our tests will mock the MACE calculator to return predictable, deterministic energies and forces. This allows us to test the sampling logic in isolation. For instance, we can provide a small, pre-defined trajectory of `ase.Atoms` objects and assert that our SOAP-based clustering and sampling algorithm correctly selects the most geometrically distinct structures, thereby verifying the core of our "smart sampling" logic without needing a real ML model.

The unit tests for the `SimulationEngine` (Module E) are critical for validating the active learning trigger mechanism. We will create a mock ACE calculator using `pytest-mock`. This mock calculator will be cleverly programmed to return a low, stable uncertainty value for the first `N` simulation steps, and then a very high uncertainty value on step `N+1`. The test will then execute the `SimulationEngine` with this mock calculator and assert that the engine runs for exactly `N+1` steps and returns the specific atomic structure that caused the high uncertainty. This directly and robustly tests the uncertainty detection and structure extraction mechanism, which is the absolute heart of the on-the-fly learning loop. Finally, the `ConfigExpander` itself will be unit tested by providing it with a sample `MinimalConfig` object and asserting that the output `FullConfig` contains correctly derived, physically reasonable values.

**Integration Testing Approach:**

The integration testing for CYCLE02 culminates in a full end-to-end test of the complete autonomous workflow. This test will, by necessity, be a carefully scaled-down version of a real scientific run, specifically designed to execute quickly and reliably in an automated CI environment.

The main integration test, which we can call `test_full_autonomous_loop`, will perform the following steps:
1.  It will start with a very simple, minimal `input.yaml` file (e.g., for elemental Si, which is computationally cheap).
2.  It will invoke the main CLI entry point: `mlip-pipe run input.yaml`. This immediately tests the `ConfigExpander` and the correct initialization of the main `WorkflowOrchestrator`.
3.  The test will configure the orchestrator to run for a small, fixed number of active learning generations (e.g., `max_generations: 2`).
4.  Crucially, to make this test fast and deterministic, we will patch the most expensive external calls:
    *   The `LabelingEngine`'s `subprocess.run` call to Quantum Espresso will be patched to return pre-computed, valid DFT results for any given Si structure almost instantly.
    *   The `TrainingEngine`'s actual model fitting process will be patched to simply copy a pre-trained "dummy" potential file into the expected output location, bypassing the time-consuming training step.
5.  The test will then assert that the following sequence of events occurs in the correct order:
    *   The `StructureGenerator` runs and successfully adds initial structures to the database.
    *   The `ExplorerSampler` runs, using a real (but very fast, since the system is small) MACE model.
    *   The `LabelingEngine` is called for the structures selected by the explorer (and our mock instantly provides the results).
    *   The `TrainingEngine` is called to produce the first-generation potential (and our mock creates the file).
    *   The `SimulationEngine` starts, runs for a few steps, and (by careful design of our mocked potential/uncertainty) finds an "uncertain" structure.
    *   The orchestrator correctly loops, sending this new uncertain structure back to be labeled.
    *   The loop completes successfully after the configured number of generations.
6.  Finally, the test will assert that a final MLIP file exists in the output directory.

This comprehensive test provides a powerful validation of the `WorkflowOrchestrator`'s entire state machine and ensures that all the modules are correctly wired together to perform the autonomous data generation and retraining loop, proving that the system meets its primary architectural design objective.
