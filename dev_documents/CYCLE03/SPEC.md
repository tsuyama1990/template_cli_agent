# Cycle 3 Specification: Efficient Exploration and Performance Optimisation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

This document provides the detailed technical specification for Cycle 3 of the MLIP-AutoPipe project. With the foundational data processing pipeline established in Cycle 1 and the user-friendly automation layer in Cycle 2, this cycle confronts one of the most fundamental challenges in computational materials science: the curse of dimensionality and the associated problem of **computational efficiency**. The brute-force or heuristic-based generation of structures in Cycle 2 is a significant improvement over manual creation, but it is fundamentally unguided. It operates without knowledge of the material's underlying potential energy surface (PES), and as such, is likely to generate many structures that are either thermodynamically irrelevant (prohibitively high in energy) or informationally redundant (very similar to other structures already in the dataset). Each of these suboptimal structures consumes a precious quantum of our most expensive resource: a high-fidelity DFT calculation. The primary objective of Cycle 3 is to replace this unguided approach with an intelligent, targeted exploration strategy that maximizes the scientific value obtained from every single DFT calculation.

To achieve this, this cycle focuses on the implementation of **Module B (Explorer & Sampler)**. This module introduces a hierarchical approach to exploring the vast, high-dimensional PES of a material. The core idea is to use a pre-trained, universal "foundation model" (such as MACE, M3GNet, or CHGNet) as a fast, low-cost surrogate for DFT. These models, trained on massive databases of existing materials calculations, provide a reasonable approximation of the true PES at a fraction of the computational cost (often with speedups of 10^5 to 10^6). This allows us to perform massive, long-timescale molecular dynamics simulations under a wide range of conditions, effectively exploring millions or even billions of potential atomic configurations—a scale utterly inaccessible to direct DFT simulation. From this enormous pool of candidate structures, the system will employ the sophisticated **DIRECT (Dimensionally-reduced, Information-Rich, and Equitable Trajectory) sampling** algorithm. This algorithm featurizes the entire trajectory, uses dimensionality reduction and clustering to map out the explored configurational landscape, and then performs a stratified sampling to select a small, diverse, and highly informative subset of structures to be passed to the expensive DFT labelling engine. A crucial secondary objective of this cycle is **performance optimisation**. The data processing involved in the DIRECT algorithm—calculating descriptors for millions of structures—can itself become a computational bottleneck. Therefore, this cycle includes a dedicated effort to profile the Python codebase, identify these bottlenecks, and re-implement them as high-performance, just-in-time (JIT) compiled kernels using **Numba**, ensuring that the entire pipeline remains fast and efficient.

## 2. System Architecture

Cycle 3 introduces a fundamental change to the pipeline's data generation strategy, shifting from a "generate-then-label" model to a more sophisticated "explore-and-select-then-label" architecture. Module B is inserted as a critical new stage between the initial seeding (a now much-reduced role for Module A) and the expensive labelling (Module C). This new architecture is designed to act as an intelligent information filter, ensuring that only the most valuable structures survive to be processed by the computationally intensive downstream modules. Module A is no longer responsible for creating the bulk of the initial dataset; instead, it is re-tasked to simply provide a few diverse "seed" structures (e.g., an equilibrium crystal, a defected structure, a liquid configuration) that act as starting points for Module B's large-scale exploration.

```mermaid
graph TD
    A[User Input: input.yaml] --> B{Config Expander};
    B --> C[Full Config];
    C --> D[Module A: Structure Generator];
    D -- Produces only a few --> E[Seed Structures];
    E --> F[Module B: Explorer & Sampler];
    C -- Configures Exploration --> F;
    F -- Uses pre-trained model --> G[(Foundation Model)];
    F -- 1. Runs massive surrogate MD --> H{Trajectory (Millions of frames)};
    H -- 2. Featurizes & Samples --> I[Information-Rich Candidate Structures];
    I -- A small, diverse set --> J[Module C: Labelling Engine];
    C -- Configures DFT --> J;
    J --> K[ASE Database];
    K --> L[Module D: Training Engine];
    C -- Configures Training --> L;
    L --> M[Trained MLIP Model];

    subgraph "New in Cycle 3"
        direction LR
        F;
        G;
    end
```

**Detailed Workflow Description:**

1.  **Seeding:** The workflow begins as in Cycle 2, with the `ConfigExpander` creating a `FullConfig`. However, this configuration now instructs Module A to generate only a small number of seed structures.
2.  **Exploration and Sampling (Module B):** This new, central module receives the seed structures and executes a multi-step process:
    a.  **Surrogate-Driven MD:** It loads the specified foundation model (e.g., MACE) and uses it as an ASE-compatible calculator. This calculator is then used to run extensive, high-temperature (and/or high-pressure) molecular dynamics simulations, starting from the seed structures. If a GPU is available, the model is automatically moved to it for maximum performance. This step generates a massive trajectory, a time-ordered series of atomic configurations, which represents a deep exploration of the material's accessible phase space.
    b.  **Featurization:** The module then processes this massive trajectory. For each and every frame, it computes a high-dimensional feature vector, or "descriptor" (such as SOAP, Smooth Overlap of Atomic Positions), that uniquely represents the local atomic environment around each atom. This step transforms the geometric information of the trajectory into a large numerical matrix. This is a primary target for Numba optimization.
    c.  **Dimensionality Reduction:** The resulting descriptor matrix is often too high-dimensional to be clustered effectively. The module applies a dimensionality reduction technique, like Principal Component Analysis (PCA), to project the data into a lower-dimensional space (e.g., 2-10 dimensions) while preserving the most important structural variances.
    d.  **Clustering:** In this low-dimensional space, the module uses a clustering algorithm, such as K-Means, to group the millions of points into a manageable number of distinct clusters. Each cluster represents a unique structural motif or "family" of configurations discovered during the MD run.
    e.  **Stratified Sampling:** This is the core of the DIRECT algorithm's selection logic. Instead of sampling randomly, the module performs a stratified or "equitable" sampling. It iterates through the identified clusters and selects a number of representative frames from each one. This ensures that even rare but structurally important configurations (which might form a small cluster) are represented in the final selection, alongside common configurations from large clusters.
3.  **Labelling and Training:** The final, curated list of information-rich candidate structures is then passed to Module C for DFT labelling and Module D for MLIP training, following the exact same robust procedure established in the previous cycles. The critical difference is the quality and efficiency of the input data: the DFT engine is now exclusively focused on a small set of structures that are guaranteed to be diverse, non-redundant, and highly informative, maximizing the return on investment for every CPU hour spent.

## 3. Design Architecture

This cycle introduces a computationally intensive new module, `ExplorerSampler`, and a new `performance` module to house optimized Numba kernels. The design must accommodate GPU usage and handle large in-memory datasets efficiently.

**Key Classes and Modules:**

-   **`src/mlip_autopipe/modules/explorer_sampler.py`:**
    -   **`ExplorerSampler` class:** This class orchestrates the entire exploration and sampling process.
        -   `__init__(self, config: dict)`: Its configuration will be extensive, specifying the surrogate model to use (e.g., `mace-mp-0`), MD parameters (temperature, pressure, timestep, number of steps), and the settings for the DIRECT sampler (descriptor type, number of clusters, final number of samples to select).
        -   `run(self, seed_structures: List[Atoms]) -> List[Atoms]`: The main public method that executes the full workflow.
        -   `_load_surrogate_model(self)`: A method responsible for dynamically loading the specified foundation model from a library like `mace-torch`. It will implement the necessary adapter logic to wrap the PyTorch model in an ASE `Calculator` API. Crucially, this method will contain the logic to check for `torch.cuda.is_available()` and automatically move the model to the GPU device if possible, providing transparent acceleration.
        -   `_run_surrogate_md(self, structure: Atoms) -> Trajectory`: This method will use the ASE library's MD integrators (e.g., `Langevin`) to execute the high-temperature simulation using the loaded surrogate calculator.
        -   `_calculate_descriptors_for_trajectory(self, trajectory: Trajectory) -> np.ndarray`: This method will be responsible for the featurization step. It will call a high-performance, Numba-optimised kernel from the `performance` module to do the actual computation, passing in the trajectory data and returning a large NumPy array of descriptor vectors.
        -   `_perform_direct_sampling(self, descriptors: np.ndarray) -> List[int]`: This method implements the core sampling logic. It will use robust, battle-tested implementations from the `scikit-learn` library for PCA and KMeans clustering. It will then contain the custom logic for the final stratified sampling step, which selects the indices of the frames to be passed downstream.

-   **`src/mlip_autopipe/common/performance.py` (New Module):** This new module is dedicated to housing performance-critical code.
    -   **`@jit(nopython=True, parallel=True)` `calculate_soap_descriptors_kernel(...)` function:** This will be a key function, a JIT-compiled implementation of a descriptor calculation. By decorating a pure Python function with Numba's decorators, we can achieve performance comparable to C or Fortran without the complexities of writing and compiling extension modules. The `nopython=True` flag ensures that the entire function is compiled to native machine code with no Python interpreter overhead, while `parallel=True` allows Numba to automatically parallelize loops over the atoms or frames, fully utilizing modern multi-core CPUs.
    -   Other potential JIT-compiled functions for tasks like fast distance matrix calculations or custom clustering algorithms may also be placed here.

## 4. Implementation Approach

The implementation for Cycle 3 will proceed in a logical sequence: first, integrate the external machine learning models; second, build the data processing and sampling pipeline; and third, profile and optimize the identified performance bottlenecks.

1.  **Surrogate Model Integration:**
    a.  The first step is to update the project dependencies in `pyproject.toml` to include `torch` and the chosen foundation model library, such as `mace-torch`.
    b.  The `_load_surrogate_model` method will be implemented. This is a critical integration point. The logic will need to handle finding the pre-trained model file, loading it with PyTorch, and then wrapping it in a class that inherits from `ase.calculators.calculator.Calculator`. This wrapper class will translate the ASE API calls (e.g., `get_potential_energy`, `get_forces`) into the corresponding PyTorch model inference calls. The implementation will include the `if torch.cuda.is_available(): model.to('cuda')` logic.
    c.  The `_run_surrogate_md` method will then be implemented to use this new calculator, and it will be tested by running a short MD simulation to ensure it produces a valid ASE `Trajectory` object.

2.  **DIRECT Sampling Implementation:**
    a.  New dependencies, `scikit-learn` and a descriptor library like `dscribe`, will be added to `pyproject.toml`.
    b.  An initial, unoptimized version of `_calculate_descriptors_for_trajectory` will be implemented using the `dscribe` library directly. This allows us to develop the correct logic before focusing on performance.
    c.  The `_perform_direct_sampling` method will be implemented using standard components from `scikit-learn`: `sklearn.decomposition.PCA` for dimensionality reduction and `sklearn.cluster.KMeans` for clustering. This will be followed by a custom Python implementation of the stratified sampling logic, which selects a proportional number of samples from each identified cluster.

3.  **Performance Optimisation with Numba:**
    a.  The `numba` library will be added as a project dependency.
    b.  The `_calculate_descriptors_for_trajectory` method, implemented in the previous step, will be profiled using a tool like `cProfile` to confirm that it is indeed a performance bottleneck when run on a large trajectory. This "measure, don't guess" approach is critical for effective optimization.
    c.  A custom, simplified descriptor algorithm will be implemented as a new function in `common/performance.py`. This function will be written in standard, loop-heavy Python but will be carefully decorated with `@jit(nopython=True, parallel=True)`.
    d.  The `_calculate_descriptors_for_trajectory` method in the `ExplorerSampler` will be refactored to call this new, high-performance JIT-compiled kernel instead of the external library.
    e.  A dedicated benchmark test will be created. This test will use `pytest-benchmark` to measure and compare the execution time of the Numba version against the original library-based version, and it will be configured to fail the CI build if the performance speedup is not significant, thus preventing performance regressions.

4.  **Integration into Main Workflow:**
    a.  The main workflow orchestrator in `orchestration/workflow.py` will be modified to insert a call to the new `ExplorerSampler` after the `StructureGenerator`.
    b.  The `ConfigExpander` will be updated with the logic to add and provide sensible defaults for a new `exploration` section in the `FullConfig`, allowing the user to control the key parameters of this new phase.
    c.  A full integration test will be written. It will start from a minimal `input.yaml`, run the entire new workflow (`ConfigExpander` -> `StructureGenerator` -> `ExplorerSampler`), and feed the selected structures to a mock `LabellingEngine`, ensuring the data flows correctly through the new architecture.

## 5. Test Strategy

The test strategy for Cycle 3 must rigorously validate both the scientific correctness of the complex sampling algorithm and the tangible performance gains achieved through JIT compilation.

**Unit Testing Approach:**
-   **Explorer & Sampler:** The `ExplorerSampler` is a complex module that requires thorough unit testing of its constituent parts. The `_load_surrogate_model` method will be tested to ensure it correctly identifies and prepares the model for both CPU and potential GPU environments (the latter will be tested via mocking). The `_run_surrogate_md` method will be tested with a mock ASE calculator to assert that it runs for the exact number of steps specified in the configuration and correctly attaches the calculator to the `Atoms` object. The core logic of the DIRECT sampler (`_perform_direct_sampling`) will be tested using a carefully constructed, synthetic dataset of descriptor vectors. This synthetic dataset will be generated using `numpy.random` to have a known, multi-modal structure (e.g., three distinct blobs of points in 2D space). The test will run the sampling method on this synthetic data and assert that the final list of selected indices includes representatives from each of the three known clusters, thereby confirming that the "equitable" and "diverse" sampling goal is being met.

-   **Numba Optimisation:** The JIT-compiled functions in `common/performance.py` will be subjected to two critical types of unit tests. The first is a **correctness test**. This test will run both the Numba-jitted function and a slow, simple, pure Python equivalent on the same small input and assert that their numerical outputs are identical to within a very tight tolerance. This is crucial to ensure that the optimization process has not introduced subtle bugs. The second is a **performance test**. A separate benchmark test will be created using the `pytest-benchmark` plugin. This test will run on a much larger, more realistic dataset and will be configured to fail the CI build if the execution time of the JIT-compiled function is not at least, for example, 10 times faster than the pure Python version. This creates an automated performance guardrail, ensuring that the benefits of the optimization are real and are protected against future code changes.

**Integration Testing Approach:**
-   **Module B End-to-End Test:** A full integration test will be created to validate the entire internal workflow of Module B. This test will start with a single, simple seed structure (e.g., a small cluster of Argon atoms). It will then run a short but real surrogate MD simulation using a real, loaded MACE model. It will proceed to calculate descriptors for the resulting trajectory and perform the final sampling step. The test will make several assertions to confirm success: that a valid ASE `Trajectory` object was created, that the number of selected structures exactly matches the number requested in the test's configuration, and that the selected structures are not all identical (e.g., by checking their RMSD against each other), which proves that the sampling is not just picking frames from the beginning of the trajectory. This test provides high confidence that all the internal components of Module B are working together correctly.

-   **GPU Usage Test:** To validate the GPU acceleration feature, a special integration test will be created. This test will be marked with a custom pytest marker (e.g., `@pytest.mark.gpu`) and will be configured to run only on CI runners that are equipped with a CUDA-enabled NVIDIA GPU. The test will initiate the `ExplorerSampler`, and while the surrogate MD phase is running, it will use Python's `subprocess` module to run the `nvidia-smi` command-line tool. It will then parse the output of `nvidia-smi` and assert that the Python process corresponding to the test runner is listed as having active GPU compute and memory usage. This provides direct, undeniable verification that the GPU acceleration feature is being correctly engaged on compatible hardware.
