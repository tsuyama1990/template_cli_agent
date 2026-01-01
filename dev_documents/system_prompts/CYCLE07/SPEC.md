# Specification: CYCLE07 - Advanced Active Learning & Boundary Treatment

## 1. Summary

CYCLE07 focuses on refining and improving the robustness of the active learning loop developed in the previous cycle. While the foundational OTF mechanism is in place, this cycle addresses two sophisticated challenges that can impact the quality of the training data and the efficiency of the learning process. The first goal is to implement a *dynamic uncertainty threshold*. A fixed threshold is often too restrictive in the early stages of learning and too lenient in the later stages. A dynamic threshold that adapts to the model's current maturity will make the active learning process more efficient.

The second, and more critical, goal is to implement advanced boundary treatment for the structures trapped during OTF simulations. Simply cutting a cluster of atoms out of a large periodic simulation can introduce serious physical artefacts, such as dangling bonds in covalent materials or charge imbalances in ionic systems. These artefacts can corrupt the DFT training data and poison the MLIP. This cycle will implement a suite of techniques—including buffer regions, force masking, and atomic passivation—to ensure that the structures sent to the `LabelingEngine` are physically meaningful and free from unphysical boundary effects. The successful completion of this cycle will elevate the pipeline's data quality, leading to more accurate and robust final models.

## 2. System Architecture

This cycle's changes are primarily concentrated within the `SimulationEngine` and the `WorkflowOrchestrator`, enhancing the logic of the active learning loop.

**File Structure (CYCLE07 Focus):**

Files to be created or modified are in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── cli.py
│       ├── **workflow.py**         # Orchestrator manages dynamic threshold
│       ├── **config.py**           # Config for dynamic threshold and boundary treatment
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── structure_generator.py
│           ├── explorer_sampler.py
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── **simulation_engine.py** # Implement boundary treatments
└── tests/
    ├── conftest.py
    ├── unit/
    │   └── **test_simulation_engine_advanced.py** # New tests for boundary treatment
    └── integration/
        └── test_active_learning_loop.py
```

**Component Breakdown:**

*   **`config.py`**: The `ActiveLearningConfig` Pydantic model will be updated.
    *   The `uncertainty_threshold` field will be changed to `Union[float, str]`. This allows the user to provide a fixed float or a string like `"dynamic_95percentile"`.
    *   A new `BoundaryTreatmentConfig` model will be added, with boolean flags like `use_buffer_region: bool`, `use_force_masking: bool`, and `use_passivation: bool`.
*   **`workflow.py`**: The `WorkflowOrchestrator`'s active learning loop will be updated to manage the dynamic threshold.
    *   Before starting the simulation in each generation, it will calculate the current uncertainty distribution on the *existing* training data.
    *   It will compute the 95th percentile (or other configured value) of this distribution.
    *   It will then pass this *newly calculated float value* as the threshold to the `SimulationEngine` for that specific generation's run.
*   **`modules/simulation_engine.py`**: The `SimulationEngine` will be significantly enhanced.
    *   The `run_otf_md` method will now accept the `threshold` as an argument from the orchestrator.
    *   When a high-uncertainty structure is trapped, it is not returned immediately. Instead, it is passed to a new private method, `_apply_boundary_treatments(atoms)`.
    *   `_apply_boundary_treatments` will be a dispatcher that, based on the config and material type, applies a series of cleaning operations:
        1.  **Buffer Region:** It will select a larger sphere of atoms around the high-uncertainty atom, creating a core/shell structure.
        2.  **Passivation:** If the material is covalent, it will detect any bonds that were cut at the edge of the sphere and attach hydrogen atoms to passivate them.
        3.  **Force Masking:** It will store information about which atoms are in the buffer region (the shell) in the `atoms.info` dictionary. The `LabelingEngine` will later use this information to ignore forces on these atoms.
    *   Only this "cleaned" `ase.Atoms` object is returned to the orchestrator.
*   **`tests/unit/test_simulation_engine_advanced.py`**: A new test file to unit-test the boundary treatment logic in isolation.

## 3. Design Architecture

The design of this cycle focuses on creating a flexible and configurable data cleaning pipeline for trapped structures and a responsive, adaptive control loop.

*   **Dynamic Threshold Logic:**
    *   The responsibility for calculating the threshold is placed in the `WorkflowOrchestrator`. This is because the orchestrator has access to the full training dataset and can decide *when* to re-calculate the threshold (once per generation).
    *   The orchestrator will query the database for all labelled structures, run the current MLIP model on each of them to get their uncertainty scores, and then use `numpy.percentile` to find the desired threshold.
    *   This design keeps the `SimulationEngine` simple: its job is only to respect the threshold it is given for a particular run.

*   **Boundary Treatment Pipeline:**
    *   The `_apply_boundary_treatments` method in the `SimulationEngine` will act as a mini-pipeline.
    *   **Input:** An `ase.Atoms` object representing the large simulation cell and the index of the atom with the highest uncertainty.
    *   **Step 1: Identify Core and Buffer:** Find all atoms within a "core" radius of the central atom and all atoms within a larger "buffer" radius.
    *   **Step 2: Create Cluster:** Create a new, non-periodic `ase.Atoms` object containing only the core and buffer atoms.
    *   **Step 3: Passivate (if applicable):** Analyze the connectivity of the cluster. If any core atom has a broken bond at the buffer boundary, add a passivating atom (e.g., H).
    *   **Step 4: Tag Atoms:** Create a `mask` array. Set mask value to `1` for core atoms and `0` for buffer/passivation atoms. Store this mask in `cluster.info['force_mask']`.
    *   **Output:** The new, smaller, cleaned `cluster` object.
    *   This modular pipeline design allows for each step to be tested independently and for new treatment steps to be added in the future. The `LabelingEngine` will also need to be slightly modified to look for the `force_mask` in `atoms.info` and use it if present.

## 4. Implementation Approach

The two main features of this cycle can be developed and tested in parallel before being integrated.

1.  **Configuration:** Update `config.py` with the new Pydantic models and fields for dynamic thresholds and boundary treatments.
2.  **TDD for Boundary Treatment:** Create the new `test_simulation_engine_advanced.py` file.
    *   **`test_passivation`**: Create a test case with a simple covalent structure (e.g., a chunk of silicon). Manually cut a bond to create a dangling bond. Call the `_apply_boundary_treatments` method. Assert that a new hydrogen atom has been added at the correct position and bonded to the silicon atom.
    *   **`test_force_masking`**: Create a test case with a simple cluster. Call `_apply_boundary_treatments`. Assert that the returned `atoms` object has an array in its `.info['force_mask']` dictionary. Assert that the mask is `1` for atoms near the center and `0` for atoms near the edge.
3.  **Implement Boundary Treatment:** Implement the logic in `_apply_boundary_treatments` within the `SimulationEngine` to make these unit tests pass.
4.  **TDD for Dynamic Threshold:** In the integration test file (`test_active_learning_loop.py`), modify the main loop test.
    *   Configure the workflow to use a dynamic threshold.
    *   The test will need to mock the MLIP model's uncertainty predictions. In Generation 1, make the uncertainties low. In Generation 2, make them higher.
    *   Assert that the threshold value passed to the mock `SimulationEngine.run_otf_md` is different in Generation 1 and Generation 2, and that it corresponds to the percentile of the mock uncertainties on the training set.
5.  **Implement Dynamic Threshold:** Implement the threshold calculation logic in the `WorkflowOrchestrator`'s `run_active_learning` loop.
6.  **Final Integration:** Ensure the `SimulationEngine` correctly receives the dynamic threshold and that the `LabelingEngine` is prepared to handle the `force_mask`. A slight modification to the `LabelingEngine`'s spec is needed: when preparing the training data, if an `atoms` object has a `force_mask` in its `.info` dict, the loss function's weight for the forces on the masked atoms should be set to zero.

## 5. Test Strategy

Testing in this cycle is about verifying the correctness of the new, sophisticated data processing and control logic.

**Unit Testing Approach (Min 300 words):**
The unit tests will be highly targeted at the new algorithms for boundary treatment.

*   **`test_simulation_engine_advanced.py`**:
    *   A key test will focus on the creation of the core/buffer cluster. It will start with a large (e.g., 10x10x10) periodic cell of a crystal. It will designate a central atom as the point of uncertainty. After calling `_apply_boundary_treatments`, the test will assert several things about the returned `cluster` object:
        1.  The object is no longer periodic (`all(cluster.pbc) == False`).
        2.  The number of atoms in the cluster is much smaller than the original cell but larger than 1.
        3.  The atom farthest from the center of the cluster is no more than the specified buffer radius away.
        4.  The `force_mask` array exists and correctly identifies the core atoms that were within the inner radius.
    *   Another test will focus specifically on covalent passivation. It will construct a diatomic molecule like Si2, place it in a box, and then create a cluster where one Si atom is "core" and the other is "buffer". This artificially breaks the bond. The test will assert that a new 'H' atom is added and its position is along the vector of the broken Si-Si bond.

**Integration Testing Approach (Min 300 words):**
The integration test will verify the dynamic threshold mechanism over a multi-generation run.

*   **`test_dynamic_threshold_adaptation`**: This test in `test_active_learning_loop.py` will simulate a scenario where the model becomes more confident over time.
    1.  **Setup**: Configure a 3-generation run with a `"dynamic_95percentile"` threshold.
    2.  **Mocking**:
        *   Mock the MLIP model's uncertainty prediction method.
        *   Mock the `TrainingEngine`, `LabelingEngine`, and `SimulationEngine` as before.
    3.  **Execution & Assertions (Generation 1)**:
        *   When the orchestrator calculates the threshold, the mock MLIP will return a set of low uncertainties on the initial data, resulting in a low threshold (e.g., `threshold_gen1 = 1.5`).
        *   Assert that the mock `SimulationEngine` is called with `threshold=1.5`.
    4.  **Execution & Assertions (Generation 2)**:
        *   When the orchestrator calculates the threshold for the next generation, configure the mock MLIP to return a wider distribution of uncertainties on the now-larger dataset. This will result in a higher threshold (e.g., `threshold_gen2 = 2.5`).
        *   Assert that the mock `SimulationEngine` is called the second time with `threshold=2.5`.
    *   This test proves that the orchestrator is correctly re-evaluating the model's confidence at each generation and feeding this adaptive threshold back into the simulation, fulfilling the core requirement of the dynamic learning loop.
