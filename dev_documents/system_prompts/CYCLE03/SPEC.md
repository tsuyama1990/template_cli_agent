# CYCLE03/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 3 of the MLIP-AutoPipe project. The previous cycles established a pipeline that can autonomously generate initial structures and train a baseline model. However, the quality of this initial model is limited by the diversity of the starting structures. This cycle aims to solve that problem by introducing a crucial optimization step: intelligent exploration and sampling. The primary objective is to make the data generation process vastly more efficient and effective by using a cheap, pre-trained universal surrogate model to explore a wide range of atomic configurations *before* committing to expensive DFT calculations. This embodies the system's core philosophy by replacing the brute-force (or human-guided) approach to finding important structures with a targeted, algorithm-driven search.

The core of this cycle is the implementation of **Module B: Explorer & Sampler**. This module will use the MACE-MP foundation model, a state-of-the-art universal potential, to run large-scale molecular dynamics (MD) simulations. These simulations will generate millions of atomic configurations, covering a much broader region of the material's phase space than the static structures from Cycle 2. From this massive dataset of "candidate" structures, the module will then use a sophisticated sampling technique called DIRECT (Descriptive, Informative, Representative, and Equidistributed Clustering-based Technique) to select a small, highly informative subset of structures. Only this curated subset will be passed to the labeling engine. Furthermore, performance is critical; key parts of the sampling algorithm, such as descriptor calculations, will be optimized with `Numba` to handle the large datasets efficiently. By the end of this cycle, MLIP-AutoPipe will be able to generate training datasets of significantly higher quality for the same computational budget, leading to far more accurate and robust MLIPs.

## 2. System Architecture

The architecture for Cycle 3 integrates the new Explorer & Sampler module into the main workflow. This module fits between the initial structure generation (Module A) and the DFT labeling (Module C), acting as an intelligent filter to improve the quality of data sent for expensive calculations.

**File Structure:**

The following ASCII tree highlights the new or significantly modified files for this cycle. The focus is on implementing Module B and the necessary utilities for descriptor calculation and sampling.

```
.
├── src/
│   └── mlip_autopipec/
│       ├── data/
│       │   └── models.py     # Modified to include new config sections
│       ├── modules/
│       │   ├── a_structure_generator.py
│       │   ├── **b_explorer_sampler.py**
│       │   ├── c_labeling_engine.py
│       │   └── d_training_engine.py
│       ├── orchestrator.py # Modified to include exploration step
│       └── utils/
│           └── **descriptors.py** # Numba-optimised descriptor calculations
│           └── dft_utils.py
├── tests/
│   ├── unit/
│   │   ├── modules/
│   │   │   └── **test_b_explorer_sampler.py**
│   │   └── utils/
│   │       └── **test_descriptors.py**
│   └── e2e/
│       └── **test_cycle03_workflow.py**
└── ...
```

**Component Blueprint:**

*   **`modules/b_explorer_sampler.py`**: This new file will contain the `ExplorerSampler` class. This class will orchestrate the entire exploration and sampling process. Its `execute` method will:
    1.  Take the initial static structures from the database.
    2.  Run high-temperature MD simulations for each structure using the MACE-MP potential via the `mace-torch` library. This will generate long trajectory files.
    3.  Iterate through the trajectory frames, calculating a structural descriptor (e.g., SOAP) for each frame using the highly optimized functions from `descriptors.py`.
    4.  Perform clustering on the descriptor data to group similar structures.
    5.  Execute the DIRECT sampling logic to select a diverse and representative set of structures from the clusters.
    6.  Save this information-rich subset of structures back into the `AseDB` with a status indicating they are ready for DFT labeling.
*   **`utils/descriptors.py`**: This new utility module will house functions related to calculating atomic environment descriptors. To handle the millions of frames from MD simulations, these functions (e.g., `calculate_soap_descriptors`) will be heavily optimized for performance. The implementation will use libraries like `dscribe` but will have key loops accelerated with `@numba.jit` to ensure maximum throughput.
*   **`orchestrator.py` (Modified)**: The `Orchestrator`'s main workflow method, `run_full_pipeline()`, will be updated to include this new stage. The execution sequence will now be: `Generate (A) -> Explore & Sample (B) -> Label (C) -> Train (D)`.
*   **`data/models.py` (Modified)**: The Pydantic configuration models will be extended to include a new section for the explorer and sampler, allowing the user to control parameters for this stage.

## 3. Design Architecture

The design for Cycle 3 focuses on efficient data handling and the seamless integration of a surrogate model into the pipeline. Pydantic schemas will ensure that the complex parameters for this new stage are well-defined and validated.

**Pydantic Schema Design:**

*   **`ExplorerSamplerConfig` (Model, part of FullConfig)**: This model will encapsulate all settings for `Module B`.
    *   `surrogate_model`: A string specifying the pre-trained model to use (e.g., `"mace_mp"`).
    *   `md_temperature`: A float representing the temperature (in Kelvin) for the MD simulation.
    *   `md_steps`: An integer for the number of steps in the exploration MD run.
    *   `sampling_strategy`: A string literal, fixed to `"direct"` for this cycle.
    *   `num_samples_to_select`: An integer specifying the target size of the final curated dataset.
    *   *Invariants*: Validation will ensure that `md_steps` and `num_samples_to_select` are positive integers.

*   **`FullConfig` (Top-level Model, modified)**: The `FullConfig` model will be updated to include the new configuration section.
    *   `system`: {...}
    *   `generation`: {...}
    *   **`exploration`**: An instance of the `ExplorerSamplerConfig` model.
    *   `dft_compute`: {...}
    *   `mlip_training`: {...}

**Data Flow and Consumers/Producers:**

1.  **Input Data**: The `ExplorerSampler` consumes the initial static structures produced by `Module A`, which it retrieves from the `AseDB`.
2.  **Surrogate Model (MACE)**: The module will dynamically load the specified pre-trained MACE model. This model is a key **producer** of data, generating millions of atomic configurations (positions, velocities, forces) during the MD simulation.
3.  **Descriptor Calculation**: The `descriptors.py` utility **consumes** the MD trajectory frames and **produces** a high-dimensional array of descriptor vectors.
4.  **DIRECT Sampling Algorithm**: The core logic within `ExplorerSampler` **consumes** the descriptor array and **produces** a list of indices corresponding to the selected frames.
5.  **Output Data**: The `ExplorerSampler` **produces** the final curated set of ASE `Atoms` objects. These are saved to the `AseDB`, overwriting or superseding the initial structure set, ready for **consumption** by `Module C`, the `LabelingEngine`.

This data flow acts as a funnel, starting with a few static structures, expanding to millions of candidates during exploration, and finally narrowing down to a few hundred or thousand highly valuable structures for the expensive DFT stage. This intelligent data reduction is the key value proposition of this cycle.

## 4. Implementation Approach

The implementation will focus on integrating the MACE model, building the sampling logic, and ensuring the process is performant.

1.  **Configuration Update**:
    *   Modify `data/models.py` to add the `ExplorerSamplerConfig` Pydantic model and include it in the `FullConfig`.
    *   Update the `ConfigExpander` service from Cycle 2 to populate this new section with sensible default values (e.g., a default MD temperature and number of steps).

2.  **Descriptor Utility Implementation**:
    *   Create `src/mlip_autopipec/utils/descriptors.py`.
    *   Implement a function to calculate descriptors for a set of `Atoms` objects. Start with a well-established descriptor like SOAP, using a library such as `dscribe`.
    *   Identify performance bottlenecks in the calculation loop and apply `Numba` JIT compilation to accelerate them.
    *   Write unit tests in `tests/unit/utils/test_descriptors.py` to verify that the descriptor calculation is correct and to benchmark its performance.

3.  **Explorer & Sampler Module Implementation**:
    *   Create the `ExplorerSampler` class in `src/mlip_autopipec/modules/b_explorer_sampler.py`.
    *   In the `execute` method, first load the MACE model using `mace-torch` APIs.
    *   Set up an ASE MD simulation (e.g., using `Langevin` dynamics) with the MACE model as the calculator. Run the simulation and save the trajectory.
    *   Load the trajectory and process it in batches. For each batch, call the utility function from `descriptors.py` to get the descriptor vectors.
    *   Implement the DIRECT sampling logic. This will involve:
        *   A clustering step (e.g., using a fast algorithm like `MiniBatchKMeans`).
        *   A selection step that picks representative structures from different clusters.
    *   Retrieve the selected `Atoms` objects from the trajectory and save them to the `AseDB`.
    *   Unit tests in `tests/unit/modules/test_b_explorer_sampler.py` will mock the MACE model and the MD run, providing a fixed trajectory file. The tests will assert that the sampling logic correctly selects a deterministic subset of structures.

4.  **Orchestration Integration**:
    *   Modify the `Orchestrator` class to include a call to the `ExplorerSampler.execute()` method at the correct point in the workflow pipeline (after `StructureGenerator` and before `LabelingEngine`).
    *   Update the end-to-end test, `test_cycle03_workflow.py`. This test will validate the entire `A -> B -> C -> D` workflow. It will need to use a mocked MACE model to avoid a real MD simulation, but it will verify that the data flows correctly from the generator, through a (mocked) exploration and sampling step, and into the (mocked) labeling and training engines.

## 5. Test Strategy

The testing for Cycle 3 must address the complexity of the MD simulation and the performance requirements of the sampling process.

**Unit Testing Approach (Min 300 words):**
*   **`descriptors.py`**: The unit tests for our descriptor utilities are crucial for ensuring both correctness and performance. In `test_descriptors.py`, we will create a small, fixed ASE `Atoms` object. We will calculate its SOAP descriptors using our implementation and compare the result against the known, correct output from the underlying library (`dscribe`) to ensure our logic is sound. Another important test will be performance-related. We will create a larger test case (e.g., a trajectory with a few hundred frames) and measure the execution time of the descriptor calculation function. We will have one version of the test for the pure Python implementation and another for the `Numba`-optimized version, asserting that the optimized version is significantly faster. This validates the effectiveness of our performance enhancements.
*   **`ExplorerSampler`**: Testing the `ExplorerSampler` requires mocking the computationally expensive parts. In `test_b_explorer_sampler.py`, we will not run a real MD simulation. Instead, our test setup will provide a small, pre-computed trajectory file as a fixture. The test will execute the module's logic on this fixed data. We will mock the MACE model loader. The primary assertions will be on the sampling algorithm's output. Given a fixed input trajectory, our DIRECT sampling implementation should produce the *exact same* list of selected frame indices every time. We will verify this deterministic behavior and also assert that the number of selected samples matches the configuration request. We can also test edge cases, such as what happens if the number of requested samples is larger than the number of frames in the trajectory.

**Integration Testing Approach (Min 300 words):**
*   **Module B and Database Integration**: We will test the interaction between the `ExplorerSampler` and the `AseDB`. The test will start with a temporary database populated by the `StructureGenerator`. The `ExplorerSampler` will then be run (with a mocked MD simulation). After it finishes, the test will connect to the database and verify its state. It should find that the initial structures have been replaced by or augmented with the new structures selected from the sampling process, and these new structures should be marked as ready for labeling. This ensures that the data handoff to the next stage of the pipeline via the database is working correctly.
*   **End-to-End Workflow Test**: The E2E test, `test_cycle03_workflow.py`, will validate the complete, updated pipeline. Using the `click.testing.CliRunner`, it will execute the main CLI command with an `input.yaml` configured for a Cycle 3 run. This test will require careful mocking.
    1.  The `StructureGenerator` (Module A) will run as normal.
    2.  The MD simulation within the `ExplorerSampler` (Module B) will be mocked to immediately return a small, fixed trajectory, avoiding a real simulation.
    3.  The DFT calculation within the `LabelingEngine` (Module C) will be mocked to return canned output.
    4.  The MLIP training within the `TrainingEngine` (Module D) will be mocked.
The assertions will verify the integrity of the entire data flow. We will check that the number of structures passed to the (mocked) `LabelingEngine` matches the number of samples selected by the `ExplorerSampler`, and that this number is different from the initial number of structures created by the `StructureGenerator`. This provides comprehensive proof that the new exploration and sampling stage is correctly integrated and is having the intended filtering effect on the pipeline.
