# CYCLE 03: SPEC.md - Surrogate-Based Exploration and Sampling

## 1. Summary

Cycle 03 introduces a critical optimization phase into the MLIP-AutoPipe workflow, designed to drastically reduce the number of expensive first-principles calculations required to build a high-quality potential. The core objective is to move from a "brute-force" labeling of all generated structures to an intelligent, "active" selection of only the most informative data points. This is achieved by implementing `Module B: Explorer & Sampler`, which leverages a pre-trained, general-purpose surrogate model (MACE) to perform rapid and extensive exploration of the material's potential energy surface.

Instead of immediately sending the initial structures from Cycle 02 for DFT labeling, the new module will first use the fast MACE potential to run large-scale molecular dynamics (MD) simulations. This will generate vast trajectories containing millions of atomic configurations, covering a much wider range of thermal and structural diversity than the initial static structures. The system will then apply the DIRECT (Diverse Representative SubseT) sampling algorithm. This involves calculating a structural descriptor (like SOAP) for each frame in the trajectory, clustering these frames in the high-dimensional descriptor space, and then selecting a small, diverse, and representative subset for the actual DFT labeling. Performance-critical parts of this process, such as custom descriptor calculations or clustering algorithms, will be optimized with Numba for maximum efficiency. By the end of this cycle, the pipeline will be significantly more cost-effective, focusing its expensive computational budget on unique and structurally important configurations.

## 2. System Architecture

This cycle activates `Module B` and integrates it into the main workflow. The primary new file is `b_explorer_sampler.py`, and the `WorkflowOrchestrator` will be updated to accommodate this new, intermediate step.

**File Structure (Cycle 03 Focus):**

```
.
├── src/
│   └── mlip_autopipec/
│       ├── modules/
│       │   ├── a_structure_generator.py
│       │   ├── **b_explorer_sampler.py** # NEW module implementation
│       │   ├── c_labeling_engine.py
│       │   └── d_training_engine.py
│       ├── **utils/**
│       │   └── **descriptor_utils.py** # NEW helper for descriptor calculations
│       └── **workflow.py**         # Modified to include the explorer/sampler step
└── pyproject.toml
```

The files marked in **bold** are the primary deliverables. `b_explorer_sampler.py` is the new module. We also introduce a new utility file, `descriptor_utils.py`, to encapsulate the logic for calculating structural descriptors, keeping the main module clean. The `workflow.py` will be modified to insert this new sampling stage into the pipeline.

## 3. Design Architecture

The design for this cycle centers on the creation of an efficient exploration and data reduction engine.

**Pydantic Schema (`configs/models.py`):**

The configuration will be updated to include a new section for this module.

*   `ExplorerSamplerConfig` (to be added to `MainConfig`):
    *   `surrogate_model: str` (e.g., "mace_mp")
    *   `md_temperature: float` (Temperature for the exploratory MD)
    *   `md_steps: int` (Number of steps in the MD simulation)
    *   `sampling_descriptor: str` (e.g., "soap")
    *   `n_samples: int` (The final number of structures to select)

**Class and Module Design:**

*   **`ExplorerSampler` (`modules/b_explorer_sampler.py`):**
    *   The constructor will accept an `ExplorerSamplerConfig` and an `AseDBWrapper` instance.
    *   Its main public method, `run()`, will execute the exploration and sampling logic.
    *   **Internal Logic:**
        1.  **Load Surrogate:** It will load the specified pre-trained MACE model.
        2.  **Run MD:** It will retrieve the initial structures from the database. For each structure, it will attach the MACE calculator and run an MD simulation using ASE's dynamics modules (e.g., `Langevin`). It will collect all frames from the resulting trajectory.
        3.  **Calculate Descriptors:** For each frame in the massive trajectory, it will call a helper function in `descriptor_utils.py` to compute the SOAP descriptor vector.
        4.  **Cluster and Sample:** It will perform clustering (e.g., using a fast algorithm like Mini-batch K-Means) on the descriptor vectors. After clustering, it will implement a stratified sampling strategy to pick a representative structure from each of the `n_samples` most populated or diverse clusters.
        5.  **Update Database:** The module will save the *selected* structures back to the database with a status like 'selected_for_labeling', ensuring the `LabelingEngine` only processes this curated subset. The original, generated structures can be marked as 'archived'.

*   **`descriptor_utils.py` (`utils/descriptor_utils.py`):**
    *   This file will contain functions for descriptor calculations.
    *   `calculate_soap_descriptors(frames)`: This function will take a list of `Atoms` objects and use a library like `dscribe` to compute SOAP vectors.
    *   **Performance:** Any custom logic or post-processing of descriptors that proves to be a bottleneck will be identified and accelerated using `@numba.jit`.

## 4. Implementation Approach

The implementation will focus on integrating the MACE model and building the sampling pipeline.

1.  **Add Dependencies:** Add `mace-torch` and `dscribe` to the `pyproject.toml` dependencies.
2.  **Update Configuration:** Add the `ExplorerSamplerConfig` to the Pydantic models in `configs/models.py` and update the `ConfigExpander` to populate it with reasonable defaults.
3.  **Implement MACE Integration:** In `b_explorer_sampler.py`, write the code to load a pre-trained MACE model and attach it as a calculator to an `ase.Atoms` object.
4.  **Implement MD Loop:** Write the logic to run the MD simulation using ASE's dynamics tools. Ensure that the trajectory frames are collected efficiently.
5.  **Implement Descriptor Calculation:** Create `descriptor_utils.py` and implement the SOAP calculation function, wrapping the `dscribe` library.
6.  **Implement Sampling:** In `ExplorerSampler`, implement the DIRECT algorithm. This involves using a library like `scikit-learn` for clustering and then writing the selection logic.
7.  **Profile and Optimize:** Use profiling tools to identify any performance bottlenecks in the descriptor calculation or clustering loops. Apply Numba decorators to custom Python functions that are identified as slow.
8.  **Update WorkflowOrchestrator:** Modify the `run` method in `WorkflowOrchestrator`. The new sequence will be:
    1.  Config Expansion (Cycle 02)
    2.  Structure Generation (Cycle 02)
    3.  **Explorer & Sampler (Cycle 03)**
    4.  Labeling Engine (Cycle 01)
    5.  Training Engine (Cycle 01)
    The orchestrator will now call `explorer_sampler.run()` after `structure_generator.run()` and before `labeling_engine.run()`.

## 5. Test Strategy

Testing will focus on the correctness of the sampling algorithm and the integration of the new module into the pipeline.

**Unit Testing Approach:**
(Located in `tests/unit/`)
*   **`modules/b_explorer_sampler.py` (`ExplorerSampler`):**
    *   This test will be crucial. We will create a small, artificial trajectory of `Atoms` objects where some structures are very similar and some are very different.
    *   Mock the MACE model and the database wrapper.
    *   Mock the descriptor calculation function in `descriptor_utils.py` to return pre-defined descriptor vectors that reflect the structural similarities in the artificial trajectory.
    *   Run the `explorer_sampler.run()` method.
    *   Assert that the sampling logic correctly selects the diverse structures and discards the redundant ones. For example, if there are 10 identical structures and 2 unique ones, the sampler should pick one of the 10 and both of the unique ones.
    *   Assert that the correct structures are saved to the mocked database with the correct status.
*   **`utils/descriptor_utils.py`:**
    *   Write a simple test to ensure the wrapper around the `dscribe` library is working correctly. Provide a simple `Atoms` object (e.g., H2O) and assert that the output descriptor has the expected shape and type.

**Integration Testing Approach:**
(Located in `tests/e2e/`)
*   **Surrogate Workflow Test:** This test will verify the new end-to-end data flow.
    *   **Setup:** Create a temporary directory and a minimal `input.yaml` for a system like Silicon.
    *   **Execution:** The test will use the `click.testing.CliRunner`. It will mock the external processes: the MACE MD run will be replaced by a very short run (e.g., a few steps), the `LabelingEngine`'s `subprocess.run` will be mocked, and the `TrainingEngine`'s training call will be mocked.
    *   **Verification:**
        1.  The test will check the database after the `ExplorerSampler` step. It should find a set of initial structures (marked 'archived') and a new, smaller set of structures (marked 'selected_for_labeling').
        2.  The number of 'selected_for_labeling' structures must match the `n_samples` parameter in the configuration.
        3.  The test will verify that the `LabelingEngine` was subsequently called only for this smaller, selected set of structures. This confirms that the data reduction is happening correctly within the pipeline.
