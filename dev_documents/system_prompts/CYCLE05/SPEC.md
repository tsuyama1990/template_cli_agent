# Specification: CYCLE05 - Advanced Sampling (DIRECT)

## 1. Summary

CYCLE05 addresses a critical challenge in the MLIP-AutoPipe workflow: data efficiency. The previous cycle, "Surrogate-Based Exploration," generated a massive number of atomic configurations, potentially millions of frames. Performing expensive DFT calculations on all of these frames would be computationally prohibitive and inefficient, as many frames are highly correlated and contain redundant information. The goal of this cycle is to implement the second half of `Module B: Explorer & Sampler`, the "Sampler," which uses an intelligent, descriptor-based algorithm called DIRECT to select a small, diverse, and maximally informative subset of structures for DFT labelling.

The scope of this cycle is to build the `Sampler` component. This component will take the raw trajectory data produced by the `Explorer`, calculate a structural descriptor (e.g., SOAP - Smooth Overlap of Atomic Positions) for each frame, and perform a clustering analysis in this high-dimensional descriptor space. Finally, it will use a stratified sampling technique to pick representative structures from different clusters, ensuring the final selection covers a wide range of atomic environments, from stable equilibrium states to rare transitional configurations. The successful completion of this cycle will mean the pipeline now has a powerful data reduction and selection mechanism, ensuring that every expensive DFT calculation provides the maximum possible value for training the MLIP. This is the core of the "intelligent sampling" philosophy of the project.

## 2. System Architecture

This cycle completes the `explorer_sampler.py` module by adding the `Sampler` class. It will require adding new dependencies for descriptor calculation and potentially for clustering algorithms.

**File Structure (CYCLE05 Focus):**

Files to be created or modified are in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── **workflow.py**         # Orchestrator now calls Sampler after Explorer
│       ├── **config.py**           # Add config options for Sampler (descriptor, num_clusters)
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── structure_generator.py
│           ├── **explorer_sampler.py** # Sampler implementation
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── simulation_engine.py
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_explorer.py
    │   └── **test_sampler.py**       # New test file
    └── integration/
        └── test_workflow.py
```

**Component Breakdown:**

*   **`pyproject.toml`**: New dependencies will be added. This will include a library for descriptor calculations, such as `dscribe` for SOAP, and potentially `scikit-learn` for clustering algorithms.
*   **`config.py`**: The `ExplorerConfig` will be renamed to `ExplorerSamplerConfig` or a new `SamplerConfig` model will be added.
    *   `SamplerConfig`: A new model will contain fields like `descriptor: str` (defaulting to `"soap"`), `num_clusters: int` (default 200), `sampling_strategy: str` (default "stratified").
*   **`modules/explorer_sampler.py`**: This file will now be completed with the `Sampler` class, or the logic will be added to the existing class.
    *   A main public method, `sample(trajectories: list[ase.Atoms])`, will be the entry point. It will take the large list of atoms from the `Explorer`.
    *   A private method, `_calculate_descriptors(atoms_list)`, will use the `dscribe` library to compute the SOAP descriptor for every atom in every structure. The average descriptor for each structure will be used.
    *   A private method, `_perform_clustering(descriptors)`, will use `scikit-learn`'s `KMeans` or a more efficient clustering algorithm to group the structures into `num_clusters`.
    *   A private method, `_stratified_sample(clusters)`, will implement the logic to pick a few representative structures from each cluster.
    *   The `sample` method will return a much smaller `list[ase.Atoms]` containing only the selected structures.
*   **`workflow.py`**: The `WorkflowOrchestrator`'s `run()` method will be updated significantly. The data flow will now be: `explored_frames = explorer.explore()` followed by `selected_frames = sampler.sample(explored_frames)`. Then, and only then, will the `selected_frames` be added to the database for labelling.
*   **`tests/unit/test_sampler.py`**: A new test file for the `Sampler` logic. These tests will use pre-computed descriptors and cluster assignments to avoid slow calculations.

## 3. Design Architecture

The design of the `Sampler` focuses on creating a modular and computationally efficient data processing pipeline.

*   **`Sampler` Class Design:**
    *   **Descriptor Calculation:** The `_calculate_descriptors` method will be optimized for performance. It will use `dscribe`'s parallelization capabilities (`n_jobs=-1`) to efficiently process large amounts of data. The SOAP parameters (cutoff, n_max, l_max) will be configurable via the Pydantic model, but with sensible defaults derived from the system's elements and cutoff radius.
    *   **Clustering:** While `KMeans` is a good starting point, for very large datasets, a more scalable algorithm like `MiniBatchKMeans` will be used to reduce computational overhead. The number of clusters will be a key parameter (`num_clusters`) that the user can tune to control the trade-off between dataset size and diversity.
    *   **Sampling Strategy:** The `_stratified_sample` method will ensure that the final dataset is well-balanced. It will sample from each cluster, potentially taking more samples from larger or more diverse clusters. This prevents the final dataset from being dominated by structures from a single, large, stable basin on the potential energy surface.
    *   **Data Flow:** The `Sampler` acts as a pure function. Its input is a large list of `ase.Atoms` objects, and its output is a much smaller list of `ase.Atoms` objects. It has no side effects and does not interact with the database directly, which makes it very easy to test.

*   **Workflow Orchestration:**
    *   The `WorkflowOrchestrator` is responsible for managing the data flow between the modules. The new sequence will be:
        1.  `StructureGenerator` -> populates DB.
        2.  `Explorer` -> reads from DB, produces a large list of frames *in memory*.
        3.  `Sampler` -> takes the in-memory list, produces a smaller list of frames *in memory*.
        4.  The `Orchestrator` then takes the final small list from the sampler and calls `db_wrapper.add_atoms()` for each structure, marking them for labelling.
    *   This design cleanly separates the exploration/sampling phase from the labelling/training phase.

## 4. Implementation Approach

The implementation will focus on building the sampling pipeline step-by-step, with thorough testing at each stage.

1.  **Dependencies and Config:** Add `dscribe` and `scikit-learn` to `pyproject.toml`. Update `config.py` with the new `SamplerConfig` Pydantic model.
2.  **TDD for Descriptor Calculation:** Create `test_sampler.py`. Write the first test, `test_descriptor_calculation`.
    *   This test will create a few simple `ase.Atoms` objects.
    *   It will call the private method `_calculate_descriptors`.
    *   It will assert that the output is a NumPy array of the correct shape (number of atoms x descriptor size).
    *   It will use a known, pre-computed SOAP vector for a simple configuration to check for numerical correctness.
3.  **Implement Descriptor Calculation:** Implement the `_calculate_descriptors` method using `dscribe` to make the test pass.
4.  **TDD for Clustering and Sampling:** Write a test `test_clustering_and_sampling`.
    *   This test will create a dummy dataset of pre-computed descriptor vectors where the cluster assignments are obvious (e.g., three distinct groups of points).
    *   It will call `_perform_clustering` and `_stratified_sample`.
    *   It will assert that the clustering algorithm correctly assigns the points to the three clusters.
    *   It will assert that the sampling method returns a list of indices that includes at least one index from each of the three groups, proving the stratified approach is working.
5.  **Implement Clustering and Sampling:** Implement the remaining private methods using `scikit-learn` to make the tests pass.
6.  **Integrate into `sample` Method:** Implement the public `sample` method, which calls the private methods in the correct order.
7.  **Update Workflow and Integration Test:**
    *   Modify the `WorkflowOrchestrator` to include the call to the `Sampler`. The data should now flow from explorer to sampler before going to the database.
    *   Update the main integration test in `test_workflow.py`. The mock for the `Explorer` will still return a list of frames. A mock for the `Sampler` will be added. The test will assert that the `Sampler.sample` method is called with the list of frames from the explorer. The mock sampler will return a smaller, fixed list of atoms. The test will then assert that only this smaller list of atoms is added to the database and passed to the `LabelingEngine`.

## 5. Test Strategy

Testing for this cycle is focused on the correctness of the data reduction pipeline and its algorithms.

**Unit Testing Approach (Min 300 words):**
The unit tests in `test_sampler.py` will be critical for verifying the numerical and logical correctness of the sampling algorithm without the overhead of real calculations.

*   **Descriptor Test:** The `test_descriptor_calculation` will use a very simple, fixed `ase.Atoms` object (e.g., a single Argon atom with a neighbor at a fixed distance). The expected SOAP vector for this configuration can be computed once manually and stored in the test. The test will then call the `_calculate_descriptors` method and assert that the computed vector is numerically close (`np.allclose`) to the stored, correct vector. This provides a strong guarantee that the interface with the `dscribe` library is being used correctly.
*   **Clustering Logic Test:** The `test_clustering_and_sampling` will not use real descriptors. It will construct a simple 2D NumPy array of points, for example: 10 points centered around `(0,0)`, 10 points centered around `(10,10)`, and only 2 points centered around `(20,20)`. This mimics a scenario with two large, common clusters and one small, rare cluster. The test will then call the clustering and sampling logic. The assertions will be:
    1.  The number of returned samples should be much smaller than the input.
    2.  The returned sample indices **must** contain at least one point from the `(20,20)` cluster. This is the key validation of the stratified approach—it ensures that rare configurations are not lost, which is something that random sampling would likely miss.

**Integration Testing Approach (Min 300 words):**
The integration test will confirm that the `Sampler` is correctly "plugged in" to the main data workflow, effectively reducing the amount of data sent to the expensive DFT stage.

*   **`test_workflow_with_sampling_step`**: This test in `test_workflow.py` will be an evolution of the previous cycle's test.
    1.  **Setup**: Minimal `input.yaml`, empty database.
    2.  **Mocking**:
        *   `Explorer.explore`: Mocked to return a list of 100 identical `ase.Atoms` objects. This simulates a long, boring MD trajectory where nothing happens.
        *   `Sampler.sample`: Mocked to inspect its input and return a list containing just *one* of those `ase.Atoms` objects.
        *   `subprocess.run` (QE): Mocked as always.
    3.  **Execution**: Run the workflow via the CLI runner.
    4.  **Assertions**:
        *   Assert that `Explorer.explore` was called and returned its list of 100 atoms.
        *   Assert that `Sampler.sample` was called, and crucially, that its input was the list of 100 atoms from the explorer.
        *   The most important assertion: Assert that the mock for `subprocess.run` (representing the `LabelingEngine`) was called **exactly once**. This proves that the `Sampler` successfully filtered the 100 redundant structures down to just one, and the orchestrator respected this decision. This confirms the entire data flow and the sampler's role as a critical data reduction filter.
