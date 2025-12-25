# Cycle 3 User Acceptance Tests: Efficient Exploration and Performance Optimisation

**Version:** 1.0.0
**Status:** Final

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 3. The primary focus of these tests is to validate the significant new capabilities introduced in **Module B (Explorer & Sampler)** from an end-user's perspective. The core value propositions to the user in this cycle are a massive increase in computational efficiency and a more intelligent, automated approach to dataset creation. These UAT scenarios are designed to confirm that these benefits are real, measurable, and reliable. We will verify that the system can correctly use pre-trained foundation models to explore a material's phase space, that the DIRECT sampling algorithm intelligently selects a diverse and valuable subset of structures, and that the performance optimizations using Numba provide a tangible speedup for the user. Successful completion of these tests will instill confidence that the pipeline is not just automated, but is making smart, efficient decisions that lead to better results with less computational waste.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C3-001** | **Successful Exploration Run:** This is the fundamental end-to-end test for the new module. It ensures that a user can configure and successfully execute a run that includes the surrogate-driven exploration and sampling phase. It validates that the new, complex data pipeline (MD -> Descriptors -> PCA -> Clustering -> Sampling) is correctly integrated and functions without crashing, producing a valid set of candidate structures for the downstream modules. This is the primary indicator of the cycle's successful implementation. | **High** |
| **UAT-C3-002** | **Correct Number of Selected Samples:** A key user control is the ability to define the trade-off between the cost of the DFT calculations and the size of the training set. This test confirms that the system precisely respects the user's choice. It validates that if the user requests a specific number of structures to be selected, the DIRECT sampling algorithm delivers exactly that number, demonstrating that the user has reliable control over the scope and cost of the subsequent, expensive labelling phase. | **High** |
| **UAT-C3-003** | **Measurable Performance Improvement:** The introduction of Numba for JIT compilation is a technical detail, but its benefit to the user—speed—is a critical feature. This test is designed to make that benefit tangible. By comparing the runtime of a key data processing step (descriptor calculation) with and without the Numba optimization, this test provides concrete, measurable proof that the optimization work translates into a significantly faster user experience, especially when dealing with the large datasets this cycle is designed to handle. | **High** |
| **UAT-C3-004** | **GPU Acceleration:** For many users, particularly those in research institutions or with modern workstations, GPU access is common. This test verifies that the pipeline can automatically detect and leverage this powerful hardware. It confirms that the most computationally intensive part of this cycle, the massive surrogate MD simulation, is offloaded to the GPU when available. This is a crucial feature for making the exploration phase feasible in a reasonable timeframe. | **Medium** |
| **UAT-C3-005** | **Diversity of Selected Samples:** This test validates the "intelligence" of the DIRECT sampling algorithm. It's not enough to select the right number of samples; they must be the *right* samples. This test ensures that the selection is scientifically valuable by verifying its diversity. It confirms that the algorithm does not simply pick redundant, highly similar structures, but instead provides a broad and representative sample of all the different structural motifs and phases discovered during the exploration. This is key to building trust that the automated selection is superior to random sampling. | **Medium** |

---

## 2. Behavior Definitions

### **UAT-C3-001: Successful Exploration Run**

**Scenario:** A user wants to generate a potential for a material that undergoes a solid-liquid phase transition. They want to use the new exploration capability to automatically find important structures in both the solid and liquid phases to include in the training set.
> **GIVEN** an `input.yaml` file that is configured to enable the exploration phase.
> **AND** the configuration specifies a valid, publicly available, pre-trained surrogate model (e.g., `mace-mp-0`).
> **AND** the MD parameters are set to heat the material from a solid state to a molten state.
>
> **WHEN** the user executes the main pipeline command.
>
> **THEN** the system's log output must clearly and explicitly indicate that it is starting the "Surrogate Exploration Phase".
> **AND** the log should show real-time progress for the surrogate MD run, for example by printing the current simulation step and temperature every few seconds.
> **AND** after the MD run is complete, the log must indicate that it is starting the featurization and DIRECT sampling process.
> **AND** the subsequent "DFT Labelling Phase" section of the log must show that it is processing the exact number of structures selected by the sampler.
> **AND** the entire pipeline must complete successfully from start to finish and produce a final, trained MLIP model file.

### **UAT-C3-002: Correct Number of Selected Samples**

**Scenario:** A user is working on a tight computational budget and wants to limit the number of expensive DFT calculations. They decide that for their initial run, they can only afford to label exactly 50 new structures.
> **GIVEN** an `input.yaml` where the `exploration` section is configured with the parameter `num_samples_to_select: 50`.
>
> **WHEN** the user runs the pipeline, and the exploration and sampling phase completes.
>
> **THEN** the log file must contain a clear, explicit, and easily searchable statement confirming the outcome, such as: "INFO: DIRECT sampling complete. Selected 50 information-rich structures for DFT labelling."
> **AND** when the `LabellingEngine` (Module C) begins its work, the log must show that it is processing a total of 50 structures (e.g., showing progress like "Labelling structure 1/50", "Labelling structure 2/50", etc.). This provides direct confirmation that the user's constraint was respected.

### **UAT-C3-003: Measurable Performance Improvement**

**Scenario:** A user is working with a large system, and the surrogate MD run has produced a very large trajectory file containing 100,000 frames. They want to be sure that the system can process this large amount of data efficiently.
> **GIVEN** a large trajectory file (`large.traj`) containing 100,000 frames.
> **AND** the user has access to two versions of the `explorer_sampler.py` module: one that uses a standard library for descriptor calculation, and one that uses the new Numba-optimised kernel.
>
> **WHEN** the user first runs the descriptor calculation and sampling step on `large.traj` using the unoptimized, standard library version, and measures the wall-clock time for this step.
> **AND** the user then runs the exact same step on the same file using the Numba-enabled version and measures the new wall-clock time.
>
> **THEN** the measured wall-clock time for the Numba-enabled version must be at least 10 times shorter than the time for the unoptimized version.
> **AND** the final list of selected frame indices produced by both versions must be identical, proving that the significant performance increase was achieved without sacrificing the correctness of the result.

### **UCT-C3-004: GPU Acceleration**

**Scenario:** A user is running the pipeline on their lab's workstation, which is equipped with a powerful NVIDIA RTX 4090 GPU. They expect the software to automatically use this hardware to accelerate the most demanding computations.
> **GIVEN** the user is working on a machine with a compatible NVIDIA GPU, and the drivers and CUDA-enabled version of PyTorch are correctly installed.
> **AND** the pipeline is configured to run the surrogate exploration phase, which involves a long MD simulation.
>
> **WHEN** the user starts the pipeline from their terminal.
>
> **THEN** the initial log messages, before the MD run begins, must contain a clear statement confirming that the hardware was detected, for example: "INFO: CUDA device 'NVIDIA RTX 4090' detected. Moving surrogate model to GPU for accelerated performance."
> **AND** during the surrogate MD run, if the user opens a new terminal and runs the `nvidia-smi` command, the output table must list the Python process running the pipeline and show that it has significant GPU memory allocated and a non-zero GPU-Util percentage.
> **AND** the total time taken for the surrogate MD run should be dramatically shorter than if the user were to force the same run on a CPU-only machine.

### **UAT-C3-05: Diversity of Selected Samples**

**Scenario:** A user is studying the melting of aluminum. The surrogate MD simulation has explored the crystalline solid phase at low temperature, a "slushy" mixed phase near the melting point, and a fully liquid phase at high temperature. The user needs to be confident that the DIRECT sampler is selecting structures representative of all three important regimes.
> **GIVEN** the exploration phase has completed on the aluminum melting simulation and has selected a set of 100 structures.
>
> **WHEN** the user runs a post-processing analysis script that takes the selected structures and the original MD trajectory as input.
> **AND** the script computes a structural fingerprint (like the radial distribution function, RDF) for each of the 100 selected structures.
>
> **THEN** the set of 100 RDFs must not all be identical.
> **AND** the analysis must show that some of the RDFs exhibit the sharp peaks characteristic of a crystalline solid.
> **AND** the analysis must also show that some of the RDFs exhibit the broad, smooth features characteristic of a liquid.
> **AND** this qualitative difference proves that the sampling algorithm did not get stuck in one region of the phase space, but successfully identified and extracted diverse structures representing the different physical states encountered during the simulation, thus creating a well-balanced and scientifically valuable training set.
