# CYCLE 02: SPEC.md - Structure Generation and Configuration Expansion

## 1. Summary

Cycle 02 builds directly upon the core engine established in Cycle 01, shifting focus to the crucial initial stages of the workflow: configuration simplification and automated data generation. The primary goal of this cycle is to enhance user experience and further the project's core philosophy of "removing the human expert from the loop." This is achieved by implementing two key features: the `ConfigExpander` and the `StructureGenerator`. The `ConfigExpander` will introduce a "Two-Tier Configuration" strategy, allowing users to provide a minimal `input.yaml` file with only the essential information (e.g., elements). The system will then use a heuristic engine to expand this into the full, detailed `exec_config_dump.yaml` required by the pipeline, automatically determining optimal DFT parameters and other settings.

The second major feature, the `StructureGenerator` (Module A), automates the creation of the initial training dataset. Instead of requiring the user to provide a set of atomic structures, this module will programmatically generate a diverse and physically relevant set of configurations. For this cycle, the implementation will focus on two important classes of materials: alloys and molecules. It will use the Special Quasirandom Structures (SQS) method for alloys to model compositional disorder, and Normal Mode Sampling (NMS) for molecules to explore vibrational degrees of freedom. By the end of this cycle, the user will be able to initiate a complete workflow—from structure creation to model training—by providing just a few lines of configuration, significantly lowering the barrier to entry and streamlining the entire process.

## 2. System Architecture

This cycle introduces a new functional module, `a_structure_generator.py`, and requires significant additions to the configuration models and the workflow orchestrator to handle the new "expansion" step.

**File Structure (Cycle 02 Focus):**

```
.
├── src/
│   └── mlip_autopipec/
│       ├── cli.py              # Modified to accept minimal config
│       ├── configs/
│       │   └── **models.py**       # Add ConfigExpander and minimal input models
│       ├── modules/
│       │   ├── **a_structure_generator.py** # NEW module implementation
│       │   ├── c_labeling_engine.py
│       │   └── d_training_engine.py
│       └── **workflow.py**         # Modified to include expansion & generation steps
└── pyproject.toml
```

The files marked in **bold** are the primary deliverables for this cycle. The `cli.py` will be updated to handle the new, simpler user input. The `configs/models.py` file will be expanded to include the Pydantic model for the minimal `input.yaml` and the logic for the `ConfigExpander`. The `a_structure_generator.py` file is a new addition. Finally, the `WorkflowOrchestrator` will be modified to incorporate these new components at the beginning of its execution flow.

## 3. Design Architecture

The design of this cycle focuses on creating intelligent, automated front-end components that seamlessly integrate with the existing core engine.

**Pydantic Schema (`configs/models.py`):**

We will introduce a new set of Pydantic models for the minimal user input (`input.yaml`).

*   `UserInputSystem`:
    *   `elements: List[str]`
    *   `composition: Optional[Dict[str, float]] = None`
    *   `structure_file: Optional[str] = None`

*   `UserInputSimulation`:
    *   `temperature: List[int]`

*   `MinimalInputConfig`:
    *   `system: UserInputSystem`
    *   `simulation: UserInputSimulation`

The main addition will be the `ConfigExpander` class.

*   **`ConfigExpander` (`configs/models.py`):**
    *   The constructor will take a `MinimalInputConfig` object.
    *   A public method, `expand()`, will return a fully populated `MainConfig` object (from Cycle 01).
    *   **Heuristic Logic:**
        1.  **Material Type Classification:** It will analyze the `elements` and `composition` to classify the material as 'alloy' or 'molecule' using heuristics (e.g., presence of metallic elements suggests alloy; a small, finite number of atoms suggests a molecule). This determination is crucial for selecting the correct structure generation algorithm.
        2.  **DFT Parameter Automation:** It will use the elements list to look up recommended pseudopotentials and cutoff energies from a built-in table representing a protocol like SSSP.
        3.  **Default Parameter Injection:** It will set reasonable defaults for all other parameters not specified by the user (e.g., `smearing`, `loss_weights`, `r_cut`).

**Class and Module Design:**

*   **`StructureGenerator` (`modules/a_structure_generator.py`):**
    *   The constructor will accept a `MainConfig` object (the *expanded* config) and an `AseDBWrapper` instance.
    *   Its main public method, `run()`, will execute the generation logic.
    *   **Internal Logic:**
        1.  It will read the `structure_type` ('alloy' or 'molecule') determined by the `ConfigExpander`.
        2.  **If 'alloy'**: It will implement the SQS generation logic. This will likely involve wrapping an external tool or library (like `icet` or `ase.build.sqs`) to create a supercell of a target size (e.g., 32-64 atoms) with the specified composition. It will also generate distorted versions of the SQS cell (e.g., applying volumetric and shear strain) to create a diverse dataset.
        3.  **If 'molecule'**: It will implement NMS. This involves performing a preliminary geometry optimization using a cheap calculator (like ASE's built-in EMT), calculating the Hessian matrix to get the vibrational modes, and generating new structures by displacing the atoms along these modes.
    *   Finally, it will use the `AseDBWrapper` to save all generated structures to the database with the status 'unlabeled'.

## 4. Implementation Approach

The implementation will add the new features to the front of the existing workflow.

1.  **Update Configuration Models:** Add the new `MinimalInputConfig` Pydantic models to `src/mlip_autopipec/configs/models.py`.

2.  **Implement ConfigExpander:** Create the `ConfigExpander` class. Start by implementing the data structures for the heuristic rules (e.g., a dictionary mapping elements to SSSP cutoff recommendations). Implement the logic for material type classification and parameter injection.

3.  **Implement StructureGenerator (Alloys):** Create the `StructureGenerator` class. First, implement the SQS generation path. Identify a suitable library for SQS generation and write a wrapper to invoke it based on the configuration. Add the logic for applying strains.

4.  **Implement StructureGenerator (Molecules):** Add the NMS logic to the `StructureGenerator`. Use ASE's `Vibrations` module to calculate vibrational modes and then write the function to generate displaced structures.

5.  **Update WorkflowOrchestrator:** Modify the `run` method in `WorkflowOrchestrator`. The new sequence will be:
    1.  The orchestrator is now initialized with the `MinimalInputConfig`.
    2.  It first instantiates and runs the `ConfigExpander` to get the full `MainConfig`.
    3.  It saves the `exec_config_dump.yaml` for user inspection and reproducibility.
    4.  It instantiates the `StructureGenerator` with the full config and runs it.
    5.  It then proceeds with the existing execution flow from Cycle 01: `LabelingEngine` -> `TrainingEngine`.

6.  **Update CLI:** Modify `cli.py` so that the main command now accepts the minimal `input.yaml`. The logic for loading the full config will be replaced by the logic to load the minimal config and start the updated orchestrator. The option to provide a pre-existing structure file will be preserved for cases where the user wants to skip automated generation.

## 5. Test Strategy

Testing will focus on the new heuristic logic and the correctness of the structure generation modules.

**Unit Testing Approach:**
(Located in `tests/unit/`)
*   **`configs/models.py` (`ConfigExpander`):**
    *   Write a test for each major heuristic. For example, provide a minimal config with `elements: ["Fe", "Pt"]` and assert that the expanded config correctly identifies the `structure_type` as 'alloy'.
    *   Provide `elements: ["H", "O"]` with composition `{"H": 2, "O": 1}` and assert the type is 'molecule'.
    *   Assert that the DFT parameters (e.g., `ecutwfc`) are correctly populated based on the elemental inputs.
    *   Test edge cases, such as when a user provides an element not supported by the heuristic table.
*   **`modules/a_structure_generator.py` (`StructureGenerator`):**
    *   Test the SQS and NMS methods independently.
    *   **SQS Test:** Mock the external SQS library call. Provide a config for an alloy and assert that the `StructureGenerator` calls the mock with the correct parameters (composition, cell size). The test should also verify that the returned `Atoms` objects are correctly processed and added to the mocked database wrapper.
    *   **NMS Test:** Use a simple, known molecule (like H2O). Mock the database dependency. The test will run the NMS generation logic and assert that the correct number of structures are generated and that their geometries are physically reasonable (i.e., they are distorted versions of the original molecule).

**Integration Testing Approach:**
(Located in `tests/e2e/`)
*   **Minimal Config Workflow Test:** This test will verify the integration of the new components with the existing pipeline.
    *   **Setup:** Create a temporary directory. Write a minimal `input.yaml` for a simple alloy (e.g., CuAu).
    *   **Execution:** The test will use the `click.testing.CliRunner` to invoke the CLI. It will mock the `subprocess.run` call for the `LabelingEngine` and the training call for the `TrainingEngine` as in Cycle 01.
    *   **Verification:** After the command completes, the test will:
        1.  Assert that the `exec_config_dump.yaml` file was created and that its contents are valid and reflect the expected expansion (e.g., `structure_type: "alloy"`).
        2.  Connect to the temporary database and assert that it has been populated with multiple unlabeled structures generated by the `StructureGenerator`.
        3.  Although the rest of the pipeline is mocked, the test verifies that the `WorkflowOrchestrator` correctly executes the new steps in the right order before proceeding, confirming the end-to-end data flow from minimal input to the start of the labeling phase.
