# Cycle 3 Specification: Advanced Sampling & Optimisation

## 1. Summary

Cycle 3 is dedicated to enhancing the computational efficiency and intelligence of the MLIP-AutoPipe system. The core objective is to minimise the number of expensive DFT calculations required to build an accurate MLIP. This is achieved by implementing the **Explorer & Sampler (Module B)**, a sophisticated component designed to perform vast, low-cost explorations of the material's potential energy surface and then intelligently select only the most informative atomic configurations for high-precision DFT labeling.

This cycle introduces two critical innovations. First, it integrates a pre-trained, universal "foundation model" (such as MACE-MP or CHGNet) to act as a fast and cheap surrogate for DFT. This allows the system to run large-scale molecular dynamics simulations to explore a wide range of temperatures and pressures at a fraction of the cost of AIMD (Ab Initio Molecular Dynamics). Second, it implements an advanced sampling algorithm called DIRECT (Derivative-Informed Reduction of Configuration Space), which processes the massive amount of data from the surrogate simulations. It uses clustering and stratified sampling to pick a small, diverse subset of structures that optimally covers the explored conformational space. Furthermore, this cycle addresses performance bottlenecks by using Numba to apply Just-in-Time (JIT) compilation to computationally intensive parts of the code, such as the calculation of structural descriptors used in the sampling algorithm.

## 2. System Architecture

The architecture for Cycle 3 inserts Module B into the pipeline, positioning it between the initial structure generation (Module A) and the expensive DFT labeling (Module C).

**Module B: Explorer & Sampler**
This module acts as an intelligent filter, ensuring that the DFT engine's valuable resources are spent only on structures that provide new information.

1.  **Exploration Phase (The "Explorer"):**
    *   **Input:** The initial set of diverse structures from Module A.
    *   **Logic:** The Explorer uses a pre-trained foundation model as an MLIP. These models are trained on vast public datasets and can provide reasonably accurate energy and force predictions for a wide range of materials at very low computational cost. The system will run large-scale molecular dynamics (MD) simulations using this surrogate potential. This generates long trajectories, potentially containing millions of atomic configurations, that sample a wide range of thermodynamic conditions. If a GPU is available, this process will be accelerated by leveraging the PyTorch backend of the foundation model.
    *   **Output:** A massive trajectory file (or files) containing millions of atomic snapshots.

2.  **Sampling Phase (The "Sampler"):**
    *   **Input:** The large trajectory files from the Exploration phase.
    *   **Logic:** Processing millions of structures is computationally demanding. The Sampler's task is to select a small, representative subset (e.g., a few hundred structures) for the subsequent DFT calculations. This is done via the DIRECT sampling algorithm:
        1.  **Descriptor Calculation:** For each structure in the trajectory, it calculates a numerical "fingerprint" or descriptor (like SOAP or ACE) that represents the local atomic environment. This is a potential bottleneck and will be optimised with Numba.
        2.  **Dimensionality Reduction:** The high-dimensional descriptor space is reduced to a low-dimensional one (e.g., 2D or 3D) using techniques like PCA or UMAP.
        3.  **Clustering:** A clustering algorithm (e.g., k-means) is applied to the low-dimensional data to group similar structures together.
        4.  **Stratified Sampling:** Instead of picking structures randomly, the system selects a few representative samples from each cluster. This ensures that both common (low-energy) and rare (high-energy/transition-state) configurations are included in the training set, eliminating statistical bias.
    *   **Output:** A small, diverse list of `ase.Atoms` objects, ready to be passed to Module C for DFT labeling.

## 3. Design Architecture

The implementation of Module B will be encapsulated within its own class, with clear separation between the exploration and sampling logic.

**`modules/b_explorer_sampler.py`**:
-   **`ExplorerSampler` Class**:
    -   **`__init__(self, config: FullConfig)`**: The constructor will take the full configuration, which will specify the choice of foundation model, MD simulation parameters (temperature, pressure, duration), and sampling settings.
    -   **`run(self, initial_structures: list[ase.Atoms]) -> list[ase.Atoms]`**: The main public method. It orchestrates the entire process by calling the private methods for exploration and sampling.
    -   **`_run_exploration(self, structures: list[ase.Atoms]) -> Path`**: This method will set up and run the surrogate MD simulations. It will use a library like ASE, which can interface with various MD calculators. It will select the appropriate foundation model and run the simulation, saving the trajectory to a file. It will return the path to this trajectory file.
    -   **`_run_sampling(self, trajectory_path: Path) -> list[ase.Atoms]`**: This method implements the DIRECT algorithm. It will:
        1.  Load the trajectory.
        2.  Call a highly optimised descriptor calculation function on all frames.
        3.  Perform dimensionality reduction and clustering using libraries like Scikit-learn.
        4.  Execute the stratified sampling logic to select the final structures.
        5.  Return the final list of `ase.Atoms` objects.

**Performance Optimisation (`common/descriptors.py`)**:
A new file will be created to house the performance-critical descriptor calculation logic.
-   **`calculate_soap_descriptors_numba(atoms_list: list[ase.Atoms]) -> np.ndarray`**: This function will be decorated with `@numba.jit` or will use Numba's features for parallelisation (`@numba.njit(parallel=True)`). It will contain the raw loops for calculating the descriptors, avoiding slow Python overhead and achieving near-native performance.

## 4. Implementation Approach

1.  **Foundation Model Integration:** Research and select a suitable foundation model (e.g., MACE-MP). Write the code to load the pre-trained model and wrap it in an ASE-compatible calculator interface. This will allow it to be used seamlessly with ASE's MD simulation engines.
2.  **Exploration Logic:** Implement the `_run_exploration` method. This will involve setting up an ASE `MolecularDynamics` simulation object, attaching the foundation model calculator, and running the simulation for the specified number of steps.
3.  **Descriptor Calculation (Numba):** Implement the performance-critical descriptor calculation function. This is a key step. The initial implementation can be a pure Python version, which will then be systematically optimised with Numba. Profiling tools will be used to identify the exact bottlenecks.
4.  **Dimensionality Reduction & Clustering:** Implement the first parts of the `_run_sampling` method. Use the Scikit-learn library to perform PCA and k-means clustering. This is relatively straightforward as the library provides high-level APIs.
5.  **Stratified Sampling:** Implement the final part of the `_run_sampling` method. This involves writing the logic to iterate through the clusters and select a representative number of samples from each one.
6.  **Integration:** Update the main workflow orchestrator to call Module B after Module A and before Module C. Ensure the data flows correctly between the modules.

## 5. Test Strategy

Testing for Cycle 3 will focus on the correctness of the sampling algorithm and the performance gains from optimisation.

**Unit Testing Approach (Min 300 words):**
Unit tests will be crucial for the DIRECT sampling logic. We will create a synthetic, low-dimensional dataset where the clusters are known beforehand. For example, we can generate a 2D dataset with three distinct Gaussian blobs of points. We will write a test that runs our `_run_sampling` logic on this synthetic data and asserts that the final selection of points includes samples from all three blobs. This verifies the correctness of the clustering and stratified sampling implementation. We will also write a specific unit test for the Numba-optimised descriptor calculation. This test will run both the original pure Python version and the optimised Numba version of the function on a small set of `Atoms` objects. The test will assert two things: first, that the numerical output of both functions is identical (or very nearly identical, within a small tolerance), and second, that the Numba version is significantly faster. This can be done by timing the execution of both functions using `time.perf_counter`.

**Integration Testing Approach (Min 300 words):**
The primary integration test for this cycle will be a small end-to-end run of Module B. The test will start with a single `ase.Atoms` object (e.g., a small crystal cell). It will run a very short exploration MD simulation (e.g., just a few hundred steps) to generate a small trajectory file. It will then run this trajectory through the full sampling pipeline. The test will assert that the process completes without errors and returns a list of `ase.Atoms` objects that is smaller than the trajectory length but greater than zero. We will also check that the selected atoms have valid positions. This test ensures that all the sub-components of Module B (MD simulation, descriptor calculation, clustering, sampling) work together correctly. A second integration test will focus on the performance aspect. We will have a benchmark test that runs the sampling pipeline on a moderately sized trajectory file (e.g., 10,000 frames) and asserts that the total execution time is below a predefined threshold, proving that our Numba optimisations are effective in a real-world scenario.
