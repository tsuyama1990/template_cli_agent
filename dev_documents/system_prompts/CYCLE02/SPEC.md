# CYCLE02/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 2 of the MLIP-AutoPipe project. Building upon the core DFT and training engine established in Cycle 1, the primary objective of this cycle is to eliminate the need for users to provide their own initial atomic structures. This represents a major step towards the system's core philosophy of "removing the human expert from the loop." We will achieve this by implementing **Module A: Structure Generator**, a sophisticated component capable of generating diverse, physically plausible initial structures based on minimal user input. This cycle also introduces the crucial **Config Expander (Heuristic Engine)**, which will automatically translate a user's high-level request (e.g., "I want a potential for FePt") into a detailed, explicit configuration file for the entire pipeline.

The scope of this cycle is to create a fully-fledged "zero-to-first-model" pipeline. A user will provide a minimal `input.yaml` file containing only the chemical elements and composition. The system will then autonomously perform a series of steps: first, the Config Expander will analyze the input and generate a complete `exec_config_dump.yaml`, inferring optimal DFT parameters and workflow settings. Second, Module A will use this configuration to generate a set of initial structures using physics-based algorithms like Special Quasirandom Structures (SQS) for alloys. Finally, these newly created structures will be passed directly to the existing labeling and training engines from Cycle 1 to produce a baseline MLIP. By the end of this cycle, MLIP-AutoPipe will have evolved from a data processor into a true generative tool, significantly lowering the barrier to entry for creating new machine learning potentials and marking a significant milestone in the project's journey toward full automation.

## 2. System Architecture

The architecture for Cycle 2 expands on the foundation of Cycle 1 by introducing the new structure generation module and the heuristic configuration engine. The file structure will be augmented to accommodate these new components, and the orchestrator's logic will be updated to incorporate this new initial step.

**File Structure:**

The following ASCII tree highlights the new or significantly modified files for this cycle. Bolded files are the focus of implementation for Cycle 2.

```
.
├── src/
│   └── mlip_autopipec/
│       ├── cli.py          # Modified to accept minimal input.yaml
│       ├── **configs/        # Directory for storing default configs/templates**
│       │   └── **sssp_defaults.json** # Example file for DFT heuristics
│       ├── data/
│       │   └── models.py     # Modified to include new config sections
│       ├── modules/
│       │   ├── **a_structure_generator.py**
│       │   ├── c_labeling_engine.py
│       │   └── d_training_engine.py
│       ├── orchestrator.py # Modified to include generation step
│       └── **services/**
│           └── **__init__.py**
│           └── **config_expander.py** # The new Heuristic Engine
├── tests/
│   ├── unit/
│   │   ├── modules/
│   │   │   └── **test_a_structure_generator.py**
│   │   └── services/
│   │       └── **test_config_expander.py**
│   └── e2e/
│       └── **test_cycle02_workflow.py**
├── **input.yaml**              # User's minimal input (primary input now)
├── exec_config_dump.yaml   # Auto-generated full config (artifact now)
└── pyproject.toml
```

**Component Blueprint:**

*   **`cli.py` (Modified)**: The main CLI command will now be simplified. Instead of requiring a full configuration file, it will primarily expect the minimal `input.yaml`. The logic will be updated to first invoke the `ConfigExpander` service before proceeding with the `Orchestrator`.
*   **`services/config_expander.py`**: This new file will house the `ConfigExpander` class. This service is the "heuristic engine." Its main function, `expand_config`, will take the path to the minimal `input.yaml`, parse it, analyze the material composition, and then programmatically build the full, explicit Pydantic configuration object. It will infer parameters like DFT cutoffs (by consulting a data file like `sssp_defaults.json`), k-point densities, and reasonable default settings for training. The result will be written to `exec_config_dump.yaml`.
*   **`modules/a_structure_generator.py`**: This new file contains the `StructureGenerator` class. Its `execute` method will be responsible for creating the initial database of atomic structures. It will read the configuration to determine the material type and then call the appropriate generation algorithm. For Cycle 2, the focus will be on implementing SQS generation for alloys, potentially by wrapping an external library or implementing a simplified version. It will also generate simple strained and perturbed variations of the structures to ensure data diversity. The generated structures will be saved to the database via the `AseDBWrapper`.
*   **`orchestrator.py` (Modified)**: The `Orchestrator` will be updated with a new primary workflow method, `run_full_pipeline()`. This method will now execute the modules in the sequence: `Generate -> Label -> Train`. It will be responsible for initializing and calling the new `StructureGenerator` before proceeding with the existing engines from Cycle 1.
*   **`data/models.py` (Modified)**: The Pydantic models will be restructured to support the two-tier configuration system. A `MinimalConfig` model will be created to validate the `input.yaml`. The existing `Cycle01Config` will be renamed to `FullConfig` or similar, representing the fully expanded configuration.

## 3. Design Architecture

The design of Cycle 2 is focused on creating a clear and robust separation between user intent (the minimal config) and system execution (the full config). Pydantic schemas will be used to enforce this separation and ensure data integrity at every step.

**Pydantic Schema Design:**

*   **`MinimalSystem` (Model)**: Represents the user's input in `input.yaml`.
    *   `elements`: A list of strings (e.g., `["Fe", "Pt"]`).
    *   `composition`: A string or dictionary representing the stoichiometry (e.g., `"FePt"` or `{"Fe": 1, "Pt": 1}`).
    *   `simulation_temperature`: An optional list of floats representing the target temperature range.

*   **`StructureGeneration` (Model, part of FullConfig)**: Contains settings for `Module A`.
    *   `generation_strategy`: A string literal (e.g., `"sqs"`).
    *   `supercell_size`: An integer or list of integers specifying the target number of atoms for the generated structures.
    *   `strains`: A list of floats representing the lattice strains to be applied.

*   **`FullConfig` (Top-level Model, modified from Cycle 1)**: This is the comprehensive model that governs the entire execution.
    *   `system`: A validated system model containing atom types, composition, and inferred properties like `structure_type` (e.g., "alloy").
    *   `generation`: An instance of the `StructureGeneration` model.
    *   `dft_compute`: The `DFTCompute` model from Cycle 1, now with all fields fully populated by the `ConfigExpander`.
    *   `mlip_training`: The `MLIPTraining` model from Cycle 1.

**Data Flow and Consumers/Producers:**

1.  **Producer (User)**: The user creates the `input.yaml` file.
2.  **Consumer/Producer (`ConfigExpander`)**: This service is the central component of the new workflow.
    *   It **consumes** the `input.yaml`, parsing it into the `MinimalSystem` Pydantic model.
    *   It performs heuristic logic:
        *   Analyzes `elements` and `composition` to infer `structure_type` as "alloy".
        *   Looks up the recommended `ecutwfc` for "Fe" and "Pt" from an internal data source (e.g., a JSON file containing SSSP protocol recommendations) and selects the maximum.
        *   Determines other DFT parameters based on best practices.
        *   Selects "sqs" as the default `generation_strategy` for an alloy.
    *   It **produces** the `exec_config_dump.yaml` file and the corresponding `FullConfig` Pydantic object.
3.  **Consumer (`Orchestrator`)**: The orchestrator consumes the `FullConfig` object. It passes the `generation` section to the `StructureGenerator`.
4.  **Consumer (`StructureGenerator`)**: This module consumes the `StructureGeneration` config to guide its logic.
5.  **Producer (`StructureGenerator`)**: It produces a set of ASE `Atoms` objects, which are then passed to the `AseDBWrapper` and stored in the database, ready for consumption by the `LabelingEngine` as in Cycle 1.

This design cleanly separates the "what" (user's goal) from the "how" (the detailed execution plan), which is the essence of the heuristic engine.

## 4. Implementation Approach

The implementation will focus on first building the heuristic engine and then the structure generation module, finally integrating them into the existing workflow.

1.  **Refactor Configuration**:
    *   Modify `data/models.py` to define the new `MinimalSystem` and `StructureGeneration` Pydantic models. Rename the top-level configuration model to `FullConfig` and incorporate the new sections.
    *   Create the `configs` directory and add a data file (e.g., `sssp_defaults.json`) with recommended pseudopotential cutoffs for a few elements to test the heuristic logic.

2.  **Implement ConfigExpander Service**:
    *   Create the `ConfigExpander` class in `src/mlip_autopipec/services/config_expander.py`.
    *   Implement the `expand_config` method. It should load the `input.yaml`, validate it with the `MinimalSystem` model, perform the heuristic calculations to determine all parameters for the `FullConfig`, and instantiate the `FullConfig` object.
    *   Write thorough unit tests in `tests/unit/services/test_config_expander.py` to verify that for a given minimal input, the generated full configuration contains the correct and expected values.

3.  **Implement StructureGenerator Module**:
    *   Create the `StructureGenerator` class in `src/mlip_autopipec/modules/a_structure_generator.py`.
    *   The `execute` method will be the main entry point. It should take the `StructureGeneration` config and an `AseDBWrapper` instance.
    *   Implement the SQS generation logic. This may involve using a library like `ase.build.sqs` or `icet`. The goal is to generate one or more supercells with the target composition.
    *   Implement logic to apply a series of strains (e.g., -2%, -1%, 0%, 1%, 2%) to the generated SQS cell, creating a set of new structures.
    *   Save all generated structures to the database using `AseDBWrapper`.
    *   Unit tests in `tests/unit/modules/test_a_structure_generator.py` will verify that the module produces the correct number of structures with the expected properties (composition, strain).

4.  **Integrate into Orchestrator and CLI**:
    *   Modify `orchestrator.py` to add the `StructureGenerator` to the main workflow. The `run_full_pipeline` method will now call the modules in the order A -> C -> D.
    *   Modify `cli.py` to change the user interaction flow. The command should now primarily use the `ConfigExpander` to generate the full config before passing it to the `Orchestrator`.
    *   Update the end-to-end test. `tests/e2e/test_cycle02_workflow.py` will now start from a minimal `input.yaml`. It will invoke the CLI, and then assert that the `exec_config_dump.yaml` is created correctly, that the database is populated with generated structures, and that the rest of the (mocked) pipeline runs successfully.

## 5. Test Strategy

The testing for Cycle 2 will validate the new generative and heuristic capabilities while ensuring they integrate smoothly with the existing engine.

**Unit Testing Approach (Min 300 words):**
*   **`ConfigExpander`**: This is the most critical component for unit testing in this cycle. In `test_config_expander.py`, we will create several example `input.yaml` files as test fixtures (e.g., a simple binary alloy, a mono-elemental system). For each input, we will run the `expand_config` method and perform detailed assertions on the output. We will check that the `structure_type` is correctly inferred, that DFT parameters like `ecutwfc` are chosen correctly based on our test data file, that default values are populated for omitted settings, and that the final `FullConfig` Pydantic object is valid. This ensures our heuristic logic is sound and robust.
*   **`StructureGenerator`**: The tests in `test_a_structure_generator.py` will validate the structure creation algorithms. We will mock the `AseDBWrapper` to act as a sink for the generated structures. The tests will call the `execute` method with a specific `StructureGeneration` config and then inspect the structures that were "saved" to the mock database. We will assert that the correct number of structures were created (e.g., number of SQS cells * number of strains), that each structure has the correct chemical composition, and that the applied strains are correctly reflected in the simulation cell parameters. This verifies the physical correctness of the generated data.

**Integration Testing Approach (Min 300 words):**
*   **ConfigExpander and Orchestrator Integration**: We will test the flow of configuration from the new service to the orchestrator. A test will run the `ConfigExpander` to produce a `FullConfig` object, then pass this object to the `Orchestrator`. We will use mock modules to verify that the `Orchestrator` correctly initializes the `StructureGenerator` with the `generation` section of the config and the `LabelingEngine` with the `dft_compute` section. This confirms that the configuration data is correctly routed to its consuming module.
*   **Generator and Database Integration**: An integration test will run the `StructureGenerator` against a real, temporary `AseDB`. After the generator's `execute` method is called, the test will connect to the database and query it to ensure that the atomic structures have been written correctly and can be read back without errors. This validates the serialization of `Atoms` objects into the database.
*   **End-to-End Workflow Test**: The `test_cycle02_workflow.py` E2E test will provide the ultimate validation. It will be set up with a minimal `input.yaml` for a system like "Cu" (mono-elemental to keep the SQS part simple). The `click.testing.CliRunner` will execute the main CLI command. The test will perform a series of assertions:
    1.  Check that `exec_config_dump.yaml` has been created.
    2.  Check that this dump file contains the expected heuristically-determined parameters.
    3.  Connect to the temporary database and verify that it has been populated with a set of generated structures (e.g., strained versions of a bulk Cu cell).
    4.  Verify that the (mocked) `LabelingEngine` and `TrainingEngine` were subsequently called, demonstrating that the data generated by Module A flows correctly into the rest of the pipeline.
This test comprehensively validates the entire "zero-to-first-model" user story of Cycle 2.
