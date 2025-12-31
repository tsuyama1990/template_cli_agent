# SPEC.md: Cycle 03 - Efficient Exploration & Optimisation

## 1. Summary

Cycle 03 marks a crucial turning point in the development of the MLIP-AutoPipe, addressing the single greatest challenge in automated potential generation: the immense, often prohibitive, computational cost associated with thoroughly exploring a material's vast and complex potential energy surface. The primary, non-negotiable goal of this cycle is to imbue the pipeline with a high degree of intelligence and efficiency, transforming it into a "DFT-frugal" system. This will be achieved through the design, implementation, and optimisation of the **Explorer & Sampler (Module B)**, a sophisticated component engineered to rapidly scan a wide gamut of atomic configurations and, most importantly, to select only the most unique and informative structures for the subsequent, and extremely expensive, DFT labelling stage. This cycle is explicitly designated as performance-critical. It will involve handling massive datasets, potentially comprising millions of atomic configurations, and will therefore necessitate a significant and focused effort on software optimisation to prevent this exploration phase from becoming an intractable bottleneck.

The core technical strategy of Module B revolves around the use of a pre-trained, general-purpose surrogate model—specifically, the state-of-the-art MACE-MP-0 model—which will serve as a fast, reliable, and approximate calculator. This powerful tool allows the system to run large-scale molecular dynamics simulations, generating long, information-rich trajectories that effectively explore diverse temperature and pressure regimes, all without incurring any DFT cost. However, generating this vast amount of data is only half the battle; the key intellectual challenge is to avoid the naive and inefficient trap of labelling redundant configurations. To solve this, the system will implement the **DIRECT (Descriptor-Informed Representative Environmental Clustering for Trajectories)** sampling algorithm. This advanced, multi-step process involves: first, calculating a high-dimensional local atomic environment descriptor (such as SOAP) for each frame in the trajectory; second, using an efficient clustering algorithm to group these descriptors and thereby identify unique structural motifs; and third, applying a stratified sampling technique to select a small, representative, and maximally diverse set of structures for labelling.

Recognising that these descriptor calculations and clustering operations can themselves become a significant computational bottleneck when processing millions of frames, a major focus of this cycle will be on aggressive performance optimisation. Computationally intensive "hot spots" in the Python code will be systematically identified through profiling. These bottlenecks will then be targeted for acceleration using Numba's powerful Just-In-Time (JIT) compilation capabilities, translating critical code paths into high-performance machine code. Upon the successful completion of Cycle 03, the MLIP-AutoPipe will have evolved from a simple, linear data generation tool into a smart, efficient, and resourceful system that intelligently allocates its most precious computational resources, bringing it a significant and measurable step closer to its ultimate goal of fully autonomous operation.

## 2. System Architecture

In Cycle 03, the newly developed Explorer & Sampler module (Module B) is strategically inserted into the main workflow, positioning it directly between the initial structure generation phase (developed in Cycle 02) and the DFT labelling phase (developed in Cycle 01). This placement is fundamental to its purpose as an intelligent data filter, ensuring that the computationally expensive labelling engine only ever receives a small, highly-curated set of structures.

The updated architectural workflow for this phase is as follows:
1.  **Initiation**: The overall workflow begins as defined at the end of Cycle 02. The system expands the user's `input.yaml`, and the `Structure Generator (Module A)` produces an initial set of seed structures, which are then saved to the database.
2.  **Exploration Phase Trigger**: The `Orchestrator`, having completed the initialisation, now transitions the system state and initiates the exploration phase. It instantiates the `Explorer & Sampler (Module B)`, providing it with the full system configuration.
3.  **Surrogate-Model MD Simulation**: Module B's first task is to run one or more molecular dynamics simulations using the fast, pre-trained MACE model as its ASE-compatible calculator. These simulations are typically seeded from the diverse structures generated in Cycle 02 and run at various temperatures to broadly sample the conformational space. The result of this step is a very long trajectory file (or multiple files) containing a massive number of atomic configurations (frames).
4.  **DIRECT Sampling Execution**: This is the core data processing step of the module, which proceeds in a sub-sequence:
    *   **Descriptor Calculation**: The module iterates through the entire, potentially massive, trajectory. For each frame, it calculates a local atomic environment descriptor (e.g., SOAP) for a representative set of atoms. This part of the code is computationally intensive and is where the Numba optimisations will be critical.
    *   **Clustering**: The resulting high-dimensional descriptor vectors, each representing a snapshot of an atomic environment, are then clustered using an efficient, scalable algorithm (e.g., a variant of k-means or DBSCAN). The goal is to group structurally similar environments together, effectively creating a catalogue of the unique motifs discovered during the MD run.
    *   **Stratified Sampling**: Finally, a small number of representative frames are intelligently selected from each cluster. This stratified approach is crucial because it ensures that the final selection is balanced and diverse, capturing not only the common, low-energy equilibrium states but also the rare, high-energy, or transition-state configurations that are essential for a robust MLIP.
5.  **Data Persistence of Selected Candidates**: The final, carefully selected ASE Atoms objects (frames) are then passed to the `AseDB`, where they are saved with a distinct state, such as 'selected_for_labelling'. This clearly distinguishes them from the initial seed structures.
6.  **Hand-off to the Labelling Engine**: The Orchestrator, having been notified of the completion of the exploration phase, now proceeds as before. It queries the database for all structures that are ready to be labelled (which now includes this newly selected, high-value set) and begins passing them to the `Labelling Engine (Module C)` for the expensive ground-truth DFT calculations.

This refined architecture is a powerful implementation of a data "funnelling" strategy. It strategically uses a cheap-but-approximate model to perform the computationally intensive "heavy lifting" of exploration. This allows the system to reserve the precious, expensive-but-accurate DFT engine for a much smaller, more valuable set of structures. This principle of filtering and information enrichment is the key architectural pattern of this cycle and is absolutely fundamental to the overall efficiency and feasibility of the entire MLIP-AutoPipe project.

## 3. Design Architecture

The design for Cycle 03 is centered on the implementation of the self-contained `ExplorerSampler` class (Module B) and the creation of a new set of high-performance utility functions for descriptor calculations and clustering. This continues the project's adherence to modular design principles.

**New and Updated Classes/APIs:**

1.  **`mlip_autopipec.modules.explorer.ExplorerSampler`** (New Class, located in `modules/explorer.py`)
    *   **Purpose**: To completely encapsulate and manage the surrogate-model-driven exploration and the subsequent DIRECT sampling logic.
    *   **Public API**:
        *   `__init__(config: FullConfig)`: The constructor will take the full system configuration object, from which it will extract all necessary parameters, such as the path to the surrogate model, MD simulation parameters (temperature, timestep, duration), and the settings for the DIRECT sampling algorithm.
        *   `run_exploration() -> List[ase.Atoms]`: This is the main public method that the `Orchestrator` will call. It is responsible for orchestrating the entire internal workflow of the module and returning a list of the final, selected `ase.Atoms` objects.
    *   **Key Internal Methods**:
        *   `_load_surrogate_model() -> ase.calculators.Calculator`: A private method responsible for loading the pre-trained MACE model from the filesystem and wrapping it in an ASE-compatible calculator interface, ready to be attached to an `Atoms` object.
        *   `_run_md_simulation(atoms: ase.Atoms) -> ase.io.Trajectory`: This method will take a starting `Atoms` object and use standard ASE MD libraries to run a molecular dynamics simulation using the loaded surrogate model. It will return a trajectory object.
        *   `_perform_direct_sampling(trajectory: ase.io.Trajectory) -> List[ase.Atoms]`: This private method will contain the core logic of the DIRECT algorithm, coordinating the calls to the descriptor calculation and clustering utilities.

2.  **`mlip_autopipec.utils.descriptors`** (New Module)
    *   **Purpose**: To provide highly optimised, standalone functions for calculating atomic environment descriptors.
    *   **Public API**:
        *   `@numba.jit(nopython=True, parallel=True)`
        *   `calculate_soap_descriptors(positions: np.ndarray, numbers: np.ndarray, cell: np.ndarray, pbc: np.ndarray, ...) -> np.ndarray`: A high-performance function designed for speed. It will take raw NumPy arrays (positions, atomic numbers, cell vectors, etc.) as input to be compatible with Numba's `nopython` mode. This will be the main target for Numba optimisation. The implementation will likely need to be written from scratch or by carefully adapting existing pure-Python libraries to be Numba-compatible.

3.  **`mlip_autopipec.utils.clustering`** (New Module)
    *   **Purpose**: To provide efficient, scalable clustering algorithms for the high-dimensional descriptor data.
    *   **Public API**:
        *   `cluster_and_sample(descriptors: np.ndarray, n_clusters: int, samples_per_cluster: int) -> np.ndarray`: A function that takes the large N x D array of descriptors (where N is the number of frames and D is the descriptor dimension) and returns a NumPy array of integer indices, representing the frames selected by the stratified sampling process. This may also be a target for Numba optimisation, particularly for the distance calculation components of the clustering algorithm.

4.  **`mlip_autopipec.orchestrator.Orchestrator`** (Updated)
    *   **New Method**: `run_exploration_phase()`: This new public method will be added to the `Orchestrator`'s API to manage the exploration stage of the workflow. It will be called in sequence after the `run_initialization_phase`.
    *   **Updated Workflow Logic**: The main `run` method of the orchestrator will be updated to chain the phases together in the correct order: `run_initialization_phase()` -> `run_exploration_phase()` -> `run_labelling_phase()`.

**Key Dependencies**:
*   `mace-torch`: This library will be added as a core dependency in `pyproject.toml`, as it is required for loading and running the surrogate model.
*   `numba`: This will be added as a core dependency to provide the JIT compilation capabilities essential for the performance goals of this cycle.
*   A library for SOAP descriptors, such as `dscribe`, may be added as a reference or for testing, but the production implementation will rely on our custom, Numba-compatible version for maximum performance.

The design of the `ExplorerSampler` class as a self-contained unit is a key architectural choice. It encapsulates all the complexity of this new phase, meaning the `Orchestrator`'s logic remains simple and high-level: it only needs to know to trigger the exploration phase and what to do with the resulting list of selected structures. The low-level, computationally intensive code is further isolated in the `utils` directory, making it much easier to profile, optimise, and test independently of the main application logic.

## 4. Implementation Approach

The implementation for Cycle 03 will follow a deliberate, two-stage approach for the core algorithms: first, build a functionally correct reference implementation, and second, aggressively optimise it. This ensures correctness before tackling the complexities of high-performance computing.

1.  **Surrogate Model Integration**:
    *   The first step is to add `mace-torch` to the project's dependencies in `pyproject.toml` and install it.
    *   The `_load_surrogate_model` method in the `ExplorerSampler` will be implemented. This will involve using the documented MACE API to load the pre-trained MACE-MP-0 model from its standard distribution location.
    *   The `_run_md_simulation` method will then be implemented. This will leverage standard ASE libraries for molecular dynamics, such as `ase.md.velocitydistribution` for initialising temperatures and `ase.md.langevin` for thermostatting the simulation, all while using the loaded MACE calculator.

2.  **Reference Descriptor Implementation (Pure Python)**:
    *   Before any optimisation is attempted, a simple, clear, and correct pure-Python implementation for the SOAP descriptor calculation will be created in `utils/descriptors.py`. This implementation can be based on textbook algorithms or by referencing the logic in existing libraries like `dscribe`. This version will be intentionally slow but will serve as the crucial "golden standard" for verifying the numerical correctness of the optimised version later on.

3.  **Functional DIRECT Sampling Logic**:
    *   With the reference descriptor function in place, the `_perform_direct_sampling` method in the `ExplorerSampler` will be implemented.
    *   The clustering of descriptors can initially be performed using a well-tested, standard library like `scikit-learn`'s `KMeans` implementation, which is easy to use and sufficient for the functional prototype.
    *   The stratified sampling logic will then be implemented to select a specified number of items from each of the clusters returned by `KMeans`.

4.  **Early Integration into Orchestrator**:
    *   At this stage, the new `ExplorerSampler` will be wired into the `Orchestrator`. The system will now have a fully functional, end-to-end workflow from initialization to labelling, passing through the new exploration phase. It will be functionally correct, albeit likely very slow in the exploration part. This early integration is crucial as it allows for testing of the data flow and component interactions before investing heavily in optimisation.

5.  **Performance Optimisation with Numba**:
    *   A separate test and benchmarking script will be created specifically for profiling the descriptor calculation and clustering functions to get hard data on the bottlenecks.
    *   The new, Numba-accelerated function `calculate_soap_descriptors` will be created. This is the most technically demanding part of the cycle. It involves carefully translating the pure-Python logic into a form that Numba can compile effectively in `nopython` mode. This often means working directly with NumPy arrays for all data structures and avoiding unsupported Python features like dynamic-typed lists or dictionaries in the core loops. The `@numba.jit(nopython=True, parallel=True)` decorator will be applied to leverage multi-core processing.
    *   A rigorous suite of unit tests will be written to compare the output of the Numba version against the pure-Python reference version for a variety of atomic structures, ensuring perfect numerical consistency down to a tight floating-point tolerance.
    *   Finally, the call to the slow descriptor function in the `ExplorerSampler` will be replaced with a call to the new, highly optimised Numba version.
    *   If profiling indicates it is necessary, this process of profiling and targeted optimisation will be repeated for the clustering utility functions.

This phased approach—build correctly first, then make it fast—is a standard, risk-reducing software engineering practice. It ensures that the complex scientific logic of the DIRECT sampler is validated and correct before the equally complex task of performance tuning begins.

## 5. Test Strategy

The test strategy for Cycle 03 must be multi-faceted, covering not only the functional correctness of the sophisticated sampling algorithm but also the numerical correctness and, critically, the performance gains of the optimised code. The testing will involve unit tests, integration tests, and a dedicated performance benchmarking suite.

**Unit Testing Approach (Min 600 words):**

Unit tests in this cycle will be meticulously designed to verify the logic of the `ExplorerSampler` and the numerical accuracy of the utility functions in complete isolation.

*   **`ExplorerSampler`**: The surrogate MD simulation, which is a slow external process, will be completely mocked. Instead of running a real MD simulation during the test, the test setup will provide a short, pre-computed `ase.io.Trajectory` object as a direct input to the sampling method.
    *   The test for `_perform_direct_sampling` will be the most important functional test. We will create a carefully designed toy trajectory file where, for example, the first half of the frames represent a solid crystal and the second half represent a disordered liquid. The test will then assert that the clustering and sampling logic correctly identifies these two distinct structural groups and selects a balanced number of representative samples from each. This will be done by mocking the `descriptors` and `clustering` utility functions, allowing the test to focus solely on the correctness of the sampling and selection logic within the `ExplorerSampler`.

*   **`utils.descriptors`**: This is the most performance-critical and numerically sensitive part of the cycle, and it requires the most rigorous testing.
    *   A comprehensive suite of tests will be established for the Numba-optimised `calculate_soap_descriptors` function. The core of this suite will be a comparison against a pure-Python reference implementation. The tests will take several reference `ase.Atoms` objects (representing diverse environments like a single atom, a dimer, a simple bulk crystal structure, and a molecule) and calculate the SOAP descriptors for them using the slow, simple, pure-Python version.
    *   The tests will then call the Numba-optimised function with the exact same `Atoms` objects.
    *   The powerful `np.testing.assert_allclose` function will be used to assert that the output arrays of descriptors from both functions are numerically identical to within a very tight tolerance (e.g., 1e-9). This provides a very high degree of confidence that the complex optimisation process did not introduce any subtle numerical bugs.
    *   To further ensure robustness, property-based testing using the `hypothesis` library will be employed. This will generate thousands of random but physically valid atomic configurations and feed them into the Numba-optimised code. This is excellent for finding edge cases that might cause crashes or incorrect calculations that were missed by the hand-picked test cases.

**Integration Testing and Performance Benchmarking (Min 300 words):**

The integration test for Cycle 03 will confirm that Module B works correctly as part of the larger pipeline, while a separate, dedicated benchmarking process will be used to validate the crucial performance gains.

*   **Integration Test**:
    *   The test will target the `Orchestrator.run_exploration_phase` method, treating the entire phase as a black box.
    *   It will be an end-to-end test for this specific phase, starting from a temporary database that has been pre-populated with a few initial seed structures (simulating the output of Cycle 02).
    *   The test will use the *real*, unmocked MACE model to ensure that the model loading and calculator interface are working correctly. However, to keep the test duration manageable (i.e., seconds, not minutes), it will be configured to run a very short MD simulation (e.g., only 10-20 steps). It will use the real, optimised descriptor and clustering functions.
    *   The primary assertion will be on the final state of the `AseDB`. After the `run_exploration_phase` method completes, the test will query the database and verify that a new set of atoms has been successfully added with the correct state, 'selected_for_labelling'. It will also check that the number of newly added atoms is reasonable and consistent with the sampling parameters configured for the test.

*   **Performance Benchmark**:
    *   A separate benchmarking script, which is not part of the regular, fast-running test suite, will be created in a `benchmarks/` directory.
    *   This script will create or load a large, realistic `ase.io.Trajectory` (e.g., several thousand frames of a 64-atom system).
    *   It will then use a library like `timeit` to accurately measure the execution time of the pure-Python reference descriptor function and the Numba-optimised function on this large dataset. Each function will be run multiple times to get a stable average.
    *   The script will print a clear report comparing the two timings and will include an assertion that the speed-up factor is significant (e.g., `assert speedup > 50`). This provides the concrete, quantitative evidence that the primary non-functional requirement of this cycle—making the exploration phase computationally feasible—has been successfully achieved.
