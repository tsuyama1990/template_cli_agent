# Specification: Cycle 3 - Efficient Exploration & Sampling

**Version:** 1.0.0
**Status:** Final

## 1. Summary

Cycle 3 of the MLIP-AutoPipe project introduces a critical layer of intelligence and efficiency to our automated workflow. While Cycles 1 and 2 established a pipeline that can autonomously generate initial structures and process them, it does so in a brute-force manner: every generated structure is sent for expensive DFT calculation. The primary objective of Cycle 3 is to rectify this inefficiency by implementing `Module B: Explorer & Sampler`. This module will act as a sophisticated filter, drastically reducing the number of DFT calculations required by intelligently selecting only the most informative atomic configurations for labeling. The core goal is to maximize the scientific value of each DFT calculation, thereby minimizing the computational cost and time-to-solution.

The scope of this cycle is to build a two-stage exploration and sampling system. The first stage, "Exploration," will leverage a pre-trained, universal "foundation model" MLIP (such as MACE-MP). This model, while not perfectly accurate for our specific system, is good enough to run large-scale, computationally cheap molecular dynamics simulations. By simulating the material at various temperatures and pressures, we can generate millions of diverse atomic configurations, effectively exploring a vast region of the potential energy surface without a single DFT calculation.

The second stage, "Sampling," will process the massive amount of data generated during exploration. We will implement the DIRECT sampling algorithm. This involves calculating a structural descriptor for each configuration, using dimensionality reduction techniques to visualize the coverage of the conformational space, and then applying a stratified sampling strategy to select a few hundred representative structures. This ensures that the selected structures are not redundant and provide a balanced coverage of the explored space, from stable equilibrium states to high-energy transition states. The output of this cycle will be a highly curated, information-rich dataset that will be passed to the labeling engine. The implementation of Module B will transform our pipeline from a simple automated script into an intelligent, resource-aware system. This cycle is pivotal in making the MLIP-AutoPipe a practical tool for real-world research, where computational resources are often a significant constraint. The successful completion of this cycle will result in a system that can generate high-quality training data with an order of magnitude fewer DFT calculations than a naive approach.

## 2. System Architecture

The architecture for Cycle 3 involves inserting `Module B: Explorer & Sampler` into the pipeline between the structure generation (`Module A`) and DFT labeling (`Module C`) stages. This placement is strategic: Module B will take the initial, diverse-but-unrefined structures from Module A and use them to seed a much broader, more efficient exploration of the conformational space.

**Architectural Placement and Data Flow:**
The new, enhanced data flow will be as follows:
1.  **Generation:** The pipeline begins with `Module A` generating a set of initial seed structures, as developed in Cycle 2.
2.  **Exploration & Sampling:** This set of structures is now passed to the new `Module B`.
    *   **Exploration Stage:** The module uses these seed structures to initialize a series of large-scale MD simulations using a fast, pre-trained universal potential. The simulations are run under various conditions (e.g., different temperatures) as specified in the configuration. The full trajectories from these simulations, containing potentially millions of frames, are temporarily stored.
    *   **Sampling Stage:** The module then processes these trajectories. It calculates descriptors for every frame, performs dimensionality reduction and clustering, and applies a stratified sampling algorithm to select a small (e.g., 200-500) subset of frames.
3.  **Handoff:** The resulting list of `ase.Atoms` objects—now a highly curated and information-dense set—is then passed to `Module C: Labeling Engine`.
4.  **Labeling and Training:** The rest of the pipeline (Modules C and D) proceeds as before, but now it operates on a much more intelligently selected dataset, ensuring that the expensive DFT and training resources are used to maximum effect.

**Internal Architecture of `Module B: Explorer & Sampler`:**
The module will be designed to be internally modular, separating the exploration and sampling logic.
*   **A Top-Level `ExplorerSampler` Class:** This will be the public interface. It will have a main `run_workflow()` method that takes the seed structures and returns the final sampled structures.
*   **An `Explorer` Component:** This component will be responsible for running the large-scale MD simulations. It will be configured with the universal potential model file and simulation parameters (time step, temperature, etc.). It will be designed to leverage GPU acceleration for the MLIP inference if a GPU is available, maximizing throughput.
*   **A `Sampler` Component:** This component will implement the DIRECT sampling algorithm. It will be internally divided into several sub-steps:
    *   **Descriptor Calculation:** A highly optimized routine to calculate a structural descriptor (like SOAP or ACE) for each frame in the MD trajectories. This part is performance-critical and will be accelerated with Numba.
    *   **Dimensionality Reduction:** A sub-component that uses standard libraries like scikit-learn to perform PCA or UMAP on the high-dimensional descriptor data.
    *   **Clustering and Sampling:** A sub-component that uses the low-dimensional representation to perform clustering (e.g., k-means) and then samples structures from each cluster to ensure diversity.

This architecture ensures that the computationally intensive parts (MD simulation, descriptor calculation) are well-defined and can be optimized independently. It also allows for the sampling strategy to be modified or replaced in the future without affecting the exploration part of the module.

## 3. Design Architecture

The design of `Module B` will be incorporated into the existing project structure. A new file will be added for the module, and the central configuration and orchestration logic will be updated to accommodate it.

**Updated Project Structure:**
```
src/mlip_autoflow/
├── __init__.py
├── main.py
├── config/
│   └── models.py
└── modules/
    ├── a_structure_generator.py
    ├── b_explorer_sampler.py      # New file
    ├── c_labeling_engine.py
    └── d_training_engine.py
    └── utils/
        └── descriptor_kernels.py # New file for Numba-optimized functions
```

**Class and Method Definitions:**

*   **`config.models.py`**: The Pydantic models will be updated.
    *   A new `ExplorerSamplerConfig` model will be added to `FullConfig`. This will contain fields like `universal_potential_path: str`, `md_temperatures: List[int]`, `num_samples_to_select: int`, and configuration for the DIRECT algorithm (e.g., descriptor type, clustering algorithm).

*   **`modules.utils.descriptor_kernels.py`**: This new file will contain performance-critical functions accelerated with Numba.
    *   Functions for calculating descriptors or parts of the descriptor algorithm that are slow in pure Python. This isolates the heavily optimized code.

*   **`modules.b_explorer_sampler.py`**: This file will contain the main `ExplorerSampler` class.
    *   `__init__(self, config: ExplorerSamplerConfig)`: Initializes with its configuration.
    *   `run_workflow(self, seed_structures: List[ase.Atoms]) -> List[ase.Atoms]`: The main public method. It orchestrates the internal call to the exploration and then the sampling stages.
    *   `_run_exploration_md(self, seed_structures: List[ase.Atoms]) -> List[ase.io.Trajectory]`: A private method that sets up and runs the MD simulations using the universal potential. It will use the ASE interface to a suitable MD library (like `ase.md.Langevin`). It will return a list of trajectory files or objects.
    *   `_run_direct_sampling(self, trajectories: List[ase.io.Trajectory]) -> List[ase.Atoms]`: A private method that implements the DIRECT algorithm.
        1.  It will read the frames from the trajectories.
        2.  It will call optimized functions (potentially from `descriptor_kernels.py`) to calculate descriptors for all frames.
        3.  It will use `scikit-learn` to perform PCA/UMAP and k-means clustering.
        4.  It will implement the stratified sampling logic to select a final list of `ase.Atoms` objects.

*   **`main.py`**: The main orchestrator script will be updated to insert this module into the pipeline. The new sequence of operations will be: `ConfigExpander` -> `StructureGenerator` -> `ExplorerSampler` -> `LabelingEngine` -> `TrainingEngine`. The output of each stage will become the input for the next.

This design ensures a clear separation between the simulation-heavy exploration phase and the data-processing-heavy sampling phase. The offloading of performance-critical code to a separate `utils` file keeps the main module logic clean and readable.

## 4. Implementation Approach

The implementation of Cycle 3 will focus on building the explorer and sampler components and integrating them seamlessly into the existing workflow.

**Step 1: Update Configuration**
We will start by adding the `ExplorerSamplerConfig` Pydantic model. This will involve defining the schema for all the new parameters needed for this cycle. We will also update the `ConfigExpander` to add a new section to `exec_config_dump.yaml` with sensible default values for these parameters (e.g., a default set of temperatures, a standard number of samples to select). We will also need to add logic to download or locate the specified universal potential model file.

**Step 2: Implement the Exploration Component**
The core of this step is the `_run_exploration_md` method. We will use the ASE library's molecular dynamics modules, as they provide a high-level interface that is compatible with many different MLIP evaluation libraries (e.g., those with an ASE calculator interface). The implementation will need to:
1.  Load the universal potential model and create an ASE calculator object from it.
2.  Iterate through the seed structures and the configured temperatures.
3.  For each combination, set up an `ase.md.Langevin` or similar dynamics object.
4.  Run the MD simulation for a specified number of steps, saving the trajectory to a file.
We will need to add a dependency on the universal potential's library (e.g., `mace-torch`) to `pyproject.toml`.

**Step 3: Implement the Descriptor Calculation**
This is a performance-critical step. We will implement the descriptor calculation logic inside the `_run_direct_sampling` method. We will choose a descriptor (e.g., SOAP) and use an existing library (like `dscribe`) to compute it. If performance is a bottleneck, we will identify the slow, numerical parts of the calculation and create a custom, simplified version in `modules/utils/descriptor_kernels.py`, using `@numba.jit` to compile it to machine code. This will involve working with NumPy arrays for maximum efficiency.

**Step 4: Implement the Sampling Algorithm**
With the descriptors calculated, we will implement the rest of the DIRECT algorithm. This will primarily involve using the `scikit-learn` library.
1.  **Dimensionality Reduction:** Use `sklearn.decomposition.PCA` to reduce the descriptor vectors to a lower-dimensional space (e.g., 2 or 3 dimensions for visualization, or ~10 for clustering).
2.  **Clustering:** Use `sklearn.cluster.KMeans` to group the data points in the reduced space into a predefined number of clusters.
3.  **Sampling:** Implement the stratified sampling logic. This involves iterating through each cluster and selecting one or more structures. A simple approach is to select the structure closest to the cluster centroid. A more advanced approach could be to select a random point or the point with the highest score on some metric.
The result of this process is the final list of `ase.Atoms` objects to be passed to the next stage.

**Step 5: Integration into the Main Workflow**
Finally, we will update `main.py` to instantiate and run the `ExplorerSampler` module. The orchestrator will be modified to pass the output of the `StructureGenerator` to the `ExplorerSampler`'s `run_workflow` method. The list of atoms returned by this method will then be fed into the `LabelingEngine`. The CLI and logging will be updated to reflect this new stage in the pipeline.

## 5. Test Strategy

Testing for Cycle 3 must address both the correctness of the complex sampling algorithm and the integration of the new module into the larger pipeline.

**Unit Testing Approach:**
*   **`ExplorerSampler`**: Due to the computational expense and stochastic nature of MD simulations, we will not run real MD in the unit tests. Instead, we will test the two components of the module separately.
    *   **Explorer Mock Test:** We will test the `_run_exploration_md` method by mocking the ASE MD classes. The test will verify that the method correctly sets up the calculator and dynamics objects with the parameters from the config (e.g., correct temperature, number of steps).
    *   **Sampler Test:** This will be the most critical unit test. We will create a fixed, pre-generated trajectory file with a small number of simple structures. We will write a test that runs the `_run_direct_sampling` method on this fixed input. The test will assert:
        1.  That the descriptor calculation runs without errors.
        2.  That the dimensionality reduction and clustering steps produce outputs of the expected shape and type.
        3.  That the final number of sampled structures is equal to the `num_samples_to_select` parameter.
        4.  We can also perform a sanity check on the diversity of the output, for example, by ensuring that the selected structures come from different clusters.

**Integration Testing Approach:**
The integration test will ensure that `Module B` correctly fits between `Module A` and `Module C`.
*   **Test Scenario:** A "dry run" of the new three-stage input pipeline.
    1.  **Setup:** Create a minimal `input.yaml` for a simple system. Configure the `ExplorerSamplerConfig` for a very short run: a tiny number of MD steps and a small number of samples to select (e.g., 5).
    2.  **Execution:**
        *   Run `Module A: StructureGenerator` to produce a few seed structures.
        *   Pass these seeds to `Module B: ExplorerSampler` and run its workflow. This will execute a real, but very short, MD simulation using the universal potential.
        *   Take the 5 sampled `ase.Atoms` objects returned by Module B and pass them to `Module C: LabelingEngine` (again, configured with a a dummy DFT script).
    3.  **Verification:**
        *   Assert that the entire process runs without errors.
        *   Verify that the number of structures passed from Module B to Module C is exactly 5.
        *   Check the database to confirm that the dummy DFT script was executed for these 5 specific structures.

This test will validate that the data handoffs between the modules are working correctly and that the new, more complex pipeline is correctly orchestrated from end to end. It will confirm that the intelligent sampling layer is successfully integrated and ready to provide significant efficiency gains to the overall workflow.
