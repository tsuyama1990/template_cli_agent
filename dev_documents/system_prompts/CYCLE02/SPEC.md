# SPEC: CYCLE 02 - Advanced Features, Web UI, and Robustness

## 1. Summary

This document provides the detailed technical specification for the second development cycle of MLIP-AutoPipe. Building on the core engine from Cycle 01, this cycle introduces advanced scientific features to enhance the diversity and physical realism of the generated data. Key deliverables include the implementation of a hybrid MD/MC exploration engine, the integration of ZBL potential mixing for simulation stability, a more sophisticated `FarthestPointSampler` (FPS), and support for ionic materials via an `IonicGenerator`. This cycle also includes a dedicated task for performance profiling and optimization of the simulation engine and the development of a proof-of-concept Web UI for improved usability.

## 2. System Architecture

This cycle extends the existing architecture with new components and enhanced logic within the explorer.

**File Structure (Cycle 02 Focus):**

Files/directories to be created or modified are in **bold**.

```
.
└── src
    └── mlip_autopipec
        ├── **main_gui.py**         # Web UI application
        ├── generators
        │   └── **ionic.py**        # Ionic-specific generator
        ├── explorers
        │   └── md_engine.py      # Existing file, will be HEAVILY MODIFIED
        └── sampling
            └── **fps.py**          # Farthest Point Sampler
```

## 3. Implementation Tasks with Scientific Validation

1.  **Implement `generators/ionic.py`**:
    *   **Task 1a (Implementation):** Create the `IonicGenerator`.
    *   **Task 1b (Scientific Validation):** Inside the `IonicGenerator`, implement a **charge balance check**. The method must verify that the total charge of the generated structure is neutral. It will raise a `PhysicsViolationError` if the composition and oxidation states do not result in a neutral cell.

2.  **Enhance `explorers/md_engine.py`**:
    *   **Task 2a (Hybrid MD/MC):** Implement logic to perform Monte Carlo (MC) atom swaps at pre-defined intervals during an MD run. The implementation must ensure that the overall composition of the system remains constant after a swap.
    *   **Task 2b (ZBL Potential Mixing):** Create a custom ASE calculator that wraps the MLIP calculator and a ZBL calculator. This wrapper will delegate energy and force calculations to the MLIP for long-range interactions and switch to the ZBL potential for short-range interactions (when atoms are closer than a defined cutoff). This is a critical stability feature.

3.  **Implement `sampling/fps.py`**:
    *   **Task 3a (Implementation):** Create the `FarthestPointSampler`. This will require adding a dependency for calculating SOAP descriptors (e.g., `dscribe`).
    *   **Task 3b (Scientific Validation):** The FPS logic must correctly identify and select a structurally diverse subset of configurations from a larger trajectory.

4.  **Performance Profiling and Optimization**:
    *   **Task 4a (Profiling):** Profile the `md_engine.py` module under a realistic workload (e.g., a medium-sized system for several thousand steps). Identify the primary performance bottlenecks (CPU, memory, I/O).
    *   **Task 4b (Optimization):** Based on the profiling results, implement optimizations. This could include refining the parallelization strategy, improving data handling to reduce memory usage, or optimizing critical loops.

5.  **Implement `main_gui.py`**:
    *   **Task 5a (Implementation):** Develop a simple proof-of-concept Web UI using Streamlit that allows a user to upload a configuration file, trigger a pipeline run as a subprocess, and view the final generated structures.

## 4. Detailed Test Strategy

Testing for Cycle 02 focuses on the new, complex scientific features and performance improvements.

**Unit Testing Approach:**
*   **`generators/ionic.py`**:
    *   **Test 1 (Scientific Validation):** Write a test that attempts to generate an ionic structure with a non-neutral composition. Assert that this correctly raises a `PhysicsViolationError`.
*   **`explorers/md_engine.py`**:
    *   **Test 2 (Hybrid MD/MC):** Write a test for the MC swap logic. It will provide a simple `ase.Atoms` object, run the swap function, and then assert that the `collections.Counter` of atomic numbers in the structure is identical before and after the swap, proving that composition is conserved.
    *   **Test 3 (ZBL Potential):** Write a test that creates a structure with two atoms placed very close together. Calculate the energy using only the MLIP (which might be low or cause an error) and then calculate it again using the ZBL-mixed potential. Assert that the energy from the mixed potential is very high and positive, correctly reflecting the physical repulsion.
*   **`sampling/fps.py`**:
    *   **Test 4 (Scientific Validation):** Create a small, dummy trajectory where some structures are identical (zero distance in descriptor space) and others are very different. Run the FPS algorithm and assert that it correctly selects the unique, diverse structures and discards the redundant ones.

**Integration Testing Approach:**
*   **End-to-End ZBL Test**:
    *   A new integration test will be created that runs the full pipeline with a configuration designed to make a standard MD simulation fail (e.g., very high temperature).
    *   The test will be run twice: once with ZBL potential mixing disabled, and once with it enabled.
    *   The test will assert that the first run fails (or produces structures with absurdly high energies), while the second run with ZBL enabled completes successfully and produces a database with physically plausible energies. This validates the effectiveness of the ZBL integration as a stability feature.
