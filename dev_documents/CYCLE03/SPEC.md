# CYCLE03 Specification: High-Throughput Exploration and Optimisation

## 1. Summary

CYCLE03 marks a critical optimisation phase for the MLIP-AutoPipe system, addressing the efficiency of training data selection. While CYCLE02 automated the generation of initial structures, the quality of these structures was still dependent on generalized, heuristic-driven methods. This cycle introduces the **Explorer & Sampler (Module B)**, a sophisticated component designed to replace brute-force randomness with intelligent, high-throughput exploration.

The core objective is to drastically reduce the number of expensive DFT calculations required to build a robust potential. This is achieved through the "DIRECT Sampling" methodology. Instead of sending all generated structures to the DFT engine, Module B will first use a fast, pre-trained, universal MLIP (a "surrogate model" like MACE-MP or M3GNet) to run very large-scale Molecular Dynamics (MD) simulations. This allows the system to explore a vast portion of the material's potential energy surface at a tiny fraction of the cost of DFT.

From the millions of frames generated during this surrogate MD run, the Sampler component will then select a small, maximally informative subset of structures for actual DFT labeling. This selection is not random; it is based on analyzing the structural diversity of the trajectory. To handle this data-intensive analysis efficiently, this cycle will introduce significant performance optimisations using **`Numba`**. Key algorithms, such as the calculation of structural descriptors (e.g., SOAP), will be JIT-compiled to achieve speeds comparable to native C/C++ code, which is essential for processing large trajectories in a reasonable timeframe. This cycle effectively ensures that our limited budget of expensive DFT calculations is spent only on the most valuable and unique atomic configurations.

## 2. System Architecture

The introduction of Module B fundamentally changes the data generation workflow, inserting a crucial exploration and filtering step between initial structure generation (Module A) and DFT labeling (Module C).

The refined data flow is as follows:
1.  **Configuration and Seeding (as before)**: The workflow starts with a minimal `input.yaml`, which is expanded by the `ConfigExpander`. The `StructureGenerator` (Module A) is then called, but its role is now to produce only a small number of *seed* structures (e.g., a perfect crystal, a slightly rattled version, a small amorphous cell).
2.  **Exploration (Module B - New)**: The new `ExplorerSampler` module takes the seed structures as starting points.
    *   It initializes large-scale MD simulations using a pre-trained universal potential (the surrogate model). This model is loaded and used as an ASE calculator, leveraging GPU acceleration if available.
    *   The simulation runs for hundreds of thousands or millions of steps, generating a massive trajectory file.
3.  **Sampling (Module B - New)**:
    *   The `ExplorerSampler` then analyzes this trajectory. It iterates through each frame and computes a high-dimensional descriptor vector (e.g., SOAP) that represents the local atomic environment. This computationally intensive loop is heavily accelerated with `Numba`.
    *   The set of all descriptor vectors is then processed. Dimensionality reduction (e.g., PCA) and clustering (e.g., k-means) are used to identify the main structural motifs present in the trajectory.
    *   Finally, a stratified sampling algorithm selects a few hundred representative structures from these clusters, ensuring a balanced representation of diverse configurations (e.g., solid-like, liquid-like, defects, transition states).
4.  **Labeling and Training (as before)**: This small, intelligently-selected set of candidate structures is then passed to the `LabelingEngine` (Module C) for high-precision DFT calculations. The rest of the workflow (database storage, training via Module D) proceeds as in the previous cycles.

This architectural change shifts the paradigm from "generate and check" to "explore, select, and verify," which is vastly more computationally efficient.

## 3. Design Architecture

This cycle's implementation will be centered on the new `b_explorer_sampler.py` module and the integration of performance-enhancing libraries.

**Key Files and Classes:**

*   `src/mlip_autoprope/modules/b_explorer_sampler.py`:
    *   `ExplorerSampler` class: The main orchestrator.
    *   `run(structures, config)`: The main public method.
    *   `_run_surrogate_md(start_structures, md_params)`: This method will set up and run the ASE MD simulation using the surrogate model as the calculator. It will leverage ASE's `Langevin` dynamics and save the trajectory. It will also handle the logic for automatically using a GPU if `torch.cuda.is_available()`.
    *   `_calculate_descriptors_numba(frames)`: A static, `Numba`-jitted function. It will contain the performance-critical loops for calculating descriptors for all atoms in all frames. The implementation will need to be `Numba`-compatible, relying primarily on NumPy arrays and avoiding unsupported Python features.
    *   `_select_diverse_subset(descriptors)`: This method will use libraries like `scikit-learn` to perform PCA, k-means clustering, and then sample from the resulting clusters.
*   `src/mlip_autoprope/config/core.py`: The `ExecConfig` Pydantic model will be updated with a new `explorer` section. This will allow the user to control parameters such as the choice of surrogate model, the MD simulation temperature and length, the type of descriptor to use, and the number of final structures to select.
*   `src/mlip_autoprope/cli.py`: The main `run` function will be updated to instantiate and execute the `ExplorerSampler` after the `StructureGenerator`.

**New Dependencies in `pyproject.toml`:**
*   `numba`: For JIT-compilation of performance-critical code.
*   `scikit-learn`: For dimensionality reduction and clustering algorithms.
*   `mace-torch` (or a similar library): To provide the pre-trained universal potential for the surrogate exploration.
*   `dscribe` (optional, for reference): Can be used to provide a reference implementation for SOAP descriptors to test the `Numba` version against.

## 4. Implementation Approach

1.  **Dependency Integration**: Add `numba`, `scikit-learn`, and `mace-torch` to `pyproject.toml` and ensure they can be installed correctly with `uv`.
2.  **Surrogate Model Wrapper**: Implement the logic to load the pre-trained MACE model and wrap it in an object that conforms to the ASE `Calculator` API. This will allow it to be used seamlessly with ASE's dynamics engines. Implement GPU detection logic at this stage.
3.  **MD Runner**: Implement the `_run_surrogate_md` method. Test this by running a short MD simulation on a simple system and ensuring a trajectory file is correctly produced.
4.  **Descriptor Calculation**:
    *   First, implement the descriptor calculation logic in pure Python/NumPy. This version will be slow but easy to debug.
    *   Write a unit test that compares the output of this pure Python version against a trusted library (like `dscribe`) to ensure correctness.
    *   Apply the `@jit(nopython=True)` decorator from `Numba`. This will likely require refactoring the code to use `Numba`-supported features (e.g., typed lists, explicit loops).
    *   Create a benchmark to measure the execution time before and after the `Numba` optimisation.
5.  **Sampling Logic**: Implement the `_select_diverse_subset` method using `scikit-learn`. This involves creating a `Pipeline` object that chains PCA and KMeans, and then implementing the logic to sample from the resulting clusters.
6.  **Module Integration**: Tie all the methods together in the `ExplorerSampler.run()` method.
7.  **CLI Integration**: Modify the main workflow in `cli.py` to insert the call to the `ExplorerSampler` in the correct position.

This approach isolates the most complex parts—the `Numba` optimisation and the ML model integration—allowing them to be developed and tested independently before being integrated into the larger system.

## 5. Test Strategy

Testing in CYCLE03 will have a strong focus on performance and numerical correctness, alongside the usual integration checks.

### Unit Testing Approach

*   **Surrogate Model**: Write a test to load the MACE model and perform a single energy/force calculation on a known structure. Assert that the results are numerical and have the correct shape.
*   **Descriptor Correctness**: The most critical unit test of this cycle. It will run three versions of the descriptor calculation on the same input structure: (1) the pure Python implementation, (2) the `Numba`-jitted implementation, and (3) a trusted reference implementation from a library like `dscribe`. The test will assert that the numerical outputs of all three are equal within a very tight tolerance.
*   **Sampling Logic**: Create a synthetic dataset of 2D points with clear clusters. Pass this data to the `_select_diverse_subset` method and assert that it selects at least one point from each of the predefined clusters. This verifies the sampling logic is not biased to a single dense region.

### Integration Testing Approach

*   **Performance Benchmark**: Create a dedicated test script that runs the `_calculate_descriptors_numba` function on a trajectory with several thousand frames. The test will run the function with and without the `@jit` decorator (this can be done by having two versions of the function). The test will assert that the `Numba`-enabled version is at least an order of magnitude (10x) faster than the pure Python version. This test directly validates the primary optimisation goal of the cycle.
*   **GPU Usage**: Write an integration test that can only be run on a machine with a GPU. This test will run the `_run_surrogate_md` method and will inspect the log output (or mock the `torch` library) to confirm that the MACE model was moved to the 'cuda' device.
*   **Full End-to-End Test**: Execute the entire pipeline with a minimal `input.yaml`. The key points to verify for this cycle are:
    1.  The system log should clearly show the `ExplorerSampler` module being executed.
    2.  The log should indicate the number of frames generated by the surrogate MD.
    3.  The log should report the number of structures selected for DFT (e.g., "Selected 250 out of 500,000 candidates").
    4.  The number of DFT calculations performed (which can be checked by counting entries in the database) should be equal to the number of selected structures, proving the filtering is working.
    5.  The workflow must complete successfully and produce a final MLIP model.
