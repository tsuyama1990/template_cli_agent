# Specification: CYCLE02 - Configuration & Heuristics

## 1. Summary

CYCLE02 builds directly upon the foundation laid in the previous cycle by introducing a sophisticated "Two-Tier Configuration" strategy. The primary goal of this cycle is to significantly enhance user experience and automation by abstracting away the complexities of setting up DFT and simulation parameters. This will be achieved by creating a `ConfigExpander` service, a heuristic engine that takes a very simple, user-provided `input.yaml` and intelligently expands it into a comprehensive `exec_config_dump.yaml` file. This expanded file will contain all the detailed, explicit parameters required for the pipeline to run, many of which a non-expert user would not know how to set.

The core of this cycle is the development of the physics- and materials-science-based heuristics. This includes logic to analyze the chemical composition of the target system to automatically determine key properties. For example, the system will learn to identify the material type (e.g., "alloy"), predict a reasonable melting point, and use this information to set appropriate temperature ranges for simulations. It will automatically select high-quality pseudopotential protocols (SSSP) and determine the necessary DFT plane-wave cutoffs based on the elements present. It will also make intelligent decisions about magnetic properties, for instance, by enabling ferromagnetic calculations if elements like Iron are detected. The successful completion of this cycle means the user is no longer required to provide a long, complex configuration file. Instead, they can initiate a sophisticated workflow by specifying only the most essential information: the elements and composition of their material. This is a critical step towards the project's core philosophy of "removing the human expert from the loop."

## 2. System Architecture

This cycle primarily impacts the configuration and workflow components of the architecture. The `ConfigExpander` will be a new, central piece of logic that sits between the user input and the main workflow execution.

**File Structure (CYCLE02 Focus):**

Files to be modified in this cycle are marked in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── cli.py
│       ├── **workflow.py**         # Updated to call ConfigExpander
│       ├── **config.py**           # Major changes: Add Expander and full Pydantic models
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── structure_generator.py
│           ├── explorer_sampler.py
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── simulation_engine.py
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── **test_config_expander.py** # New test file
    │   └── test_database.py
    └── integration/
        └── test_workflow.py
```

**Component Breakdown:**

*   **`config.py`**: This file will see the most significant changes.
    *   **Minimal Models:** New Pydantic models will be created to represent the user's `input.yaml` (e.g., `UserInputSystem`, `UserInputSimulation`). These models will have very few required fields.
    *   **Full Models:** The existing Pydantic models (`DFTComputeConfig`, `MLIPTrainingConfig`, etc.) will be expanded to become the "full" or "execution" configuration, containing every single parameter needed for the workflow, with no optional fields.
    *   **`ConfigExpander` Class**: A new class, `ConfigExpander`, will be created. Its main method, `expand(user_input: UserInputConfig) -> FullConfig`, will contain all the heuristic logic.
*   **`workflow.py`**: The `WorkflowOrchestrator` will be modified. Its `__init__` method will no longer take a simple config path. Instead, it will be initialized with a fully-formed `FullConfig` object.
*   **`cli.py`**: The CLI entry point will be updated. It will now be responsible for:
    1.  Loading the user's `input.yaml` into the minimal Pydantic models.
    2.  Instantiating the `ConfigExpander`.
    3.  Calling the `expand()` method to generate the `FullConfig` object.
    4.  Saving the `FullConfig` to `exec_config_dump.yaml` for user inspection and traceability.
    5.  Instantiating the `WorkflowOrchestrator` with the newly created `FullConfig` object.
*   **`tests/unit/test_config_expander.py`**: A new test file dedicated to unit-testing the heuristic logic of the `ConfigExpander`.

## 3. Design Architecture

The design of this cycle revolves around creating a robust and testable heuristic engine. The data contracts defined by Pydantic are essential for ensuring that the transformation from a minimal user input to a full execution configuration is reliable and predictable.

*   **Pydantic Schema Design:**
    *   **`UserInputConfig`**: This will be the top-level model for `input.yaml`. It will have two fields: `system: UserInputSystem` and `simulation: UserInputSimulation`.
        *   `UserInputSystem`: Will contain `elements: list[str]` and `composition: str`.
        *   `UserInputSimulation`: Might contain an optional `temperature: list[int]`.
    *   **`FullConfig`**: This will be the top-level model for the expanded configuration. It will contain fully-populated instances of `SystemConfig`, `SimulationConfig`, `DFTComputeConfig`, `MLIPTrainingConfig`, etc.
        *   `SystemConfig`: Will contain expanded fields like `structure_type: str` (e.g., "alloy") and `melting_point_guess: float`.
        *   `DFTComputeConfig`: Will now have all fields as required, such as `ecutwfc: float`, `ecutrho: float`, `kpoints_density: float`, and `magnetism: str`.
    *   **Producers/Consumers:** The user is the producer of the `UserInputConfig`. The `ConfigExpander` is the consumer of `UserInputConfig` and the producer of `FullConfig`. All other modules in the pipeline are consumers of `FullConfig` or its sub-models.

*   **`ConfigExpander` Heuristic Logic:** The `expand` method will be a pure function that encapsulates the following logic steps:
    1.  **System Analysis:**
        *   Read the `elements` and `composition`.
        *   Use a materials science library (like `pymatgen` or internal logic) to classify the `structure_type` based on electronegativity differences.
        *   If magnetic elements (`Fe`, `Co`, `Ni`) are present, set `DFTComputeConfig.magnetism` to `"ferromagnetic"`. Otherwise, set it to `None`.
    2.  **DFT Parameter Heuristics:**
        *   Implement a data structure (e.g., a dictionary) that maps element symbols to their recommended `ecutwfc` and `ecutrho` from the SSSP Precision protocol.
        *   The `expand` method will look at the `elements` in the input, find the maximum required `ecutwfc` and `ecutrho` among them, and set these values in the `DFTComputeConfig`.
    3.  **Simulation Parameter Heuristics:**
        *   If the user provides a temperature range, the expander will interpolate a reasonable number of steps (e.g., 3-5 steps) within that range.
        *   If no temperature is given, it will use the `melting_point_guess` to establish a default range (e.g., from 300K up to 0.7 * T_melt).
    4.  **Defaulting:** For all other parameters not covered by a specific heuristic (e.g., `loss_weights`, `r_cut`), the `ConfigExpander` will apply sensible, hard-coded default values.

## 4. Implementation Approach

The implementation will be focused and sequential, starting with the data models and then building the logic that transforms them.

1.  **Refactor Pydantic Models:** The first step is to modify `config.py`. The existing models will be renamed and adjusted to represent the `FullConfig`. The new, simpler `UserInputConfig` models will be created.
2.  **Create `ConfigExpander` Skeleton:** In `config.py`, create the `ConfigExpander` class with an `expand` method. Initially, this method can just return a hard-coded `FullConfig` object.
3.  **Write Unit Tests First (TDD):** Create the new test file `tests/unit/test_config_expander.py`. Write the first test case, e.g., `test_expand_fept_alloy`. In this test, create a `UserInputConfig` for "FePt". Define the exact `FullConfig` you expect it to produce, including the correct `structure_type`, magnetism settings, and SSSP-derived cutoffs. The test will initially fail.
4.  **Implement Heuristic Logic:** Now, implement the logic inside the `ConfigExpander.expand` method to make the test pass. This involves adding the element data for SSSP cutoffs and writing the `if 'Fe' in elements:` logic.
5.  **Iterate with More Tests:** Add more unit tests for different scenarios: a non-magnetic system (e.g., Si), a system requiring a different SSSP cutoff (e.g., Al), a system where a temperature range is provided, and a system where it is not. Implement the corresponding heuristic logic for each test until all tests pass.
6.  **Update CLI and Workflow:** Once the `ConfigExpander` is fully tested and functional, modify `cli.py`. Update the main function to perform the load->expand->save->execute sequence. Change the `WorkflowOrchestrator` to accept the `FullConfig` object.
7.  **Integration Test:** Finally, adapt the existing integration test in `test_workflow.py`. The test setup should now involve creating a minimal `input.yaml` in a temporary directory, and the test will verify that the `exec_config_dump.yaml` is created correctly before proceeding with the (mocked) labelling and training workflow.

## 5. Test Strategy

Testing in CYCLE02 is paramount to ensure the heuristics, which replace expert knowledge, are reliable. The focus shifts heavily towards unit testing the decision-making logic of the `ConfigExpander`.

**Unit Testing Approach (Min 300 words):**
The unit tests for the `ConfigExpander` will be comprehensive and data-driven. Each test will represent a specific user input and will assert the entire structure of the resulting `FullConfig` to prevent unintended regressions.

*   **`test_magnetic_alloy`**: This test will provide an input for a system like FePt. The assertions will be highly specific:
    *   `assert full_config.system.structure_type == "alloy"`
    *   `assert full_config.dft_compute.magnetism == "ferromagnetic"`
    *   `assert full_config.dft_compute.ecutwfc == 90.0` (Assuming 90 Ry is the max required by Fe or Pt in the SSSP data).
    *   `assert full_config.dft_compute.ecutrho == 720.0`
*   **`test_nonmagnetic_covalent`**: This test will use an input like Si.
    *   `assert full_config.system.structure_type == "covalent"`
    *   `assert full_config.dft_compute.magnetism is None`
    *   `assert full_config.dft_compute.ecutwfc == 40.0` (Using the SSSP value for Si).
*   **`test_temperature_interpolation`**: This test will provide a user input with `temperature = [300, 1000]`.
    *   `assert full_config.simulation.temperature_steps == [300, 650, 1000]` (or a similar interpolated list).
*   **`test_default_temperature`**: This test will provide no temperature input.
    *   The test will assert that `full_config.simulation.temperature_steps` is a reasonable default based on the system's estimated melting point.

These tests ensure that every key piece of heuristic logic is validated independently.

**Integration Testing Approach (Min 300 words):**
The integration test will verify the new end-to-end workflow, from the minimal user input to the final (mocked) execution.

*   **`test_full_workflow_from_minimal_input`**: This test, located in `test_workflow.py`, will replace the previous integration test.
    1.  **Setup**: The test will create a temporary directory. Inside it, it will create a minimal `input.yaml` file with content like `system: {elements: [Fe, Pt], composition: FePt}`.
    2.  **Execution**: It will use `click.testing.CliRunner` to invoke the main `cdd` command, pointing to the temporary `input.yaml`. As before, `subprocess.run` and `pacemaker.fit` will be mocked.
    3.  **Assertions**:
        *   First, the test will assert that the file `exec_config_dump.yaml` now exists in the temporary directory.
        *   It will load this YAML file and parse it into the `FullConfig` Pydantic model to ensure it is valid and contains the expected expanded values (e.g., magnetism is on, cutoffs are correct). This validates the file-based part of the new workflow.
        *   Finally, it will assert that the mocked `subprocess.run` (for QE) and `pacemaker.fit` (for training) were called. This confirms that the generated `FullConfig` object was correctly passed to the orchestrator and used by the downstream modules.
This single integration test covers the entire orchestration chain, confirming that the new configuration strategy is correctly wired into the application's execution flow.
