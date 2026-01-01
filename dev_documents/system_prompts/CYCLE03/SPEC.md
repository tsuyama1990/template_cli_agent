# Specification: CYCLE03 - Initial Structure Generation

## 1. Summary

CYCLE03 marks a significant leap forward in the pipeline's automation capabilities by implementing `Module A: Structure Generator`. Until now, the workflow has assumed the existence of a pre-defined atomic structure. This cycle removes that assumption. The primary goal is to empower the system to generate a diverse and physically plausible set of initial atomic structures based solely on the chemical composition provided in the user's `input.yaml`. This is the first critical step in seeding the entire active learning workflow.

This cycle will implement several heuristic-based structure generation algorithms, each tailored to a specific type of material. The system will use the `structure_type` ("alloy", "molecule", "ionic", "covalent") determined by the `ConfigExpander` in CYCLE02 to dispatch the appropriate generation method. For alloys, we will implement a Special Quasirandom Structures (SQS) generator. For molecular materials, Normal Mode Sampling (NMS) will be used to create distorted geometries. For ionic and covalent systems, we will implement foundational versions of Ab Initio Random Structure Searching (AIRSS) and Deep Rattling, respectively.

The output of this module will be a set of `ase.Atoms` objects. These objects will be passed to the `AseDBWrapper` to be saved in the database with the status `'unlabeled'`. The successful completion of this cycle means that the MLIP-AutoPipe pipeline can now be initiated from a truly minimal input—just the chemical formula—and autonomously generate the first batch of candidate structures for the subsequent exploration and labelling stages.

## 2. System Architecture

This cycle introduces the first major piece of scientific logic into the `modules` directory. The `StructureGenerator` will become a key component called by the `WorkflowOrchestrator` at the very beginning of a new project.

**File Structure (CYCLE03 Focus):**

Files to be created or modified are in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── cli.py
│       ├── **workflow.py**         # Orchestrator now calls the StructureGenerator first
│       ├── config.py
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── **structure_generator.py** # Module A - Main implementation
│           ├── explorer_sampler.py    # Stub
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── simulation_engine.py
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_config_expander.py
    │   └── **test_structure_generator.py** # New test file
    └── integration/
        └── test_workflow.py
```

**Component Breakdown:**

*   **`modules/structure_generator.py`**: This is the core of the cycle.
    *   A `StructureGenerator` class will be created. It will be initialized with the `FullConfig` object and a reference to the `AseDBWrapper`.
    *   A main public method, `generate()`, will be the entry point. This method will check if the database is empty. If it is, it will call a private dispatch method.
    *   A private method, `_dispatcher()`, will read `config.system.structure_type` and call the appropriate generation method (e.g., `_generate_sqs_for_alloy()`).
    *   Private methods for each generation type (`_generate_sqs_for_alloy`, `_generate_nms_for_molecule`, etc.) will contain the specific logic for that material class. Each of these methods will return a list of `ase.Atoms` objects.
*   **`workflow.py`**: The `WorkflowOrchestrator`'s `run()` method will be updated. The very first action it takes will be to instantiate the `StructureGenerator` and call its `generate()` method. This ensures that a baseline set of structures exists before any other step is attempted.
*   **`pyproject.toml`**: Dependencies required for the generation algorithms will be added. This will likely include libraries like `pymatgen` for crystallographic analysis and potentially a dedicated SQS generation library if a suitable one is found.
*   **`tests/unit/test_structure_generator.py`**: A new test file to unit-test the `StructureGenerator` class. These tests will be crucial for verifying the correctness of the generated structures.

## 3. Design Architecture

The design of the `StructureGenerator` focuses on modularity and testability. The dispatcher pattern allows for the easy addition of new generation algorithms in the future, and each algorithm is implemented as a separate, self-contained method.

*   **`StructureGenerator` Class Design:**
    *   The `generate()` method will act as a gatekeeper, preventing accidental regeneration of structures if the database is already populated. This is important for resuming a workflow.
    *   The internal methods, like `_generate_sqs_for_alloy`, will encapsulate all the logic for that specific task.
    *   **SQS Logic:** This method will need to determine a suitable supercell size (e.g., 32-96 atoms) and use a library (or a built-in algorithm) to generate the SQS configuration that best mimics a random alloy for the given composition. It should also generate distorted versions of the SQS cell by applying volumetric and shear strains to enrich the initial dataset.
    *   **NMS Logic:** This method will first need to find a reasonable equilibrium geometry for the molecule (e.g., from a public database or a quick geometry optimization with a cheap classical potential). Then, it will calculate the Hessian and normal modes, and generate new structures by displacing the atoms along these mode vectors with varying amplitudes.
    *   **AIRSS/Rattling Logic:** These methods will be simpler initially. They will generate random atomic positions within a cell (respecting ionic radii for AIRSS) or take a known structure and apply significant random displacements ("rattling").
    *   **Data Contracts:** All private generation methods have a clear contract: they must return a `list[ase.Atoms]`. The `generate()` method's contract is to take this list and call `database.add_atoms(atoms, state='unlabeled')` for each one.

*   **Pydantic Schema Interaction:**
    *   The `StructureGenerator` is a primary consumer of the `SystemConfig` model produced by the `ConfigExpander`. It directly uses the `structure_type`, `elements`, and `composition` fields to make its decisions.
    *   It might also consume parts of the `SimulationConfig`, for instance, to determine how many initial structures to generate. A new field, `initial_structures_to_generate: int`, could be added to the `SimulationConfig` model with a sensible default.

## 4. Implementation Approach

The implementation will focus on one generation method at a time, ensuring each is well-tested before moving to the next.

1.  **Skeleton and Workflow Integration:** First, create the `structure_generator.py` file with the `StructureGenerator` class skeleton. Implement the public `generate()` method and the `_dispatcher()` method. Update the `WorkflowOrchestrator` to call `structure_generator.generate()` at the beginning of its `run()` method. At this point, no structures are actually generated, but the architectural integration is complete.
2.  **TDD for SQS:** Create the `test_structure_generator.py` file. Write the first test, `test_generate_sqs_for_fept`. This test will:
    *   Create a mock `FullConfig` object with `structure_type = "alloy"` for FePt.
    *   Instantiate the `StructureGenerator` with this config and a mock database wrapper.
    *   Call `generate()`.
    *   Assert that the mock database's `add_atoms` method was called a certain number of times.
    *   Inspect the `ase.Atoms` objects that were passed to the mock. Assert that they have the correct composition (e.g., 50% Fe, 50% Pt), are within the expected size range, and that some of them have strained lattice vectors.
3.  **Implement SQS:** Implement the `_generate_sqs_for_alloy` method. This may involve integrating a third-party library like `icet` or `ase.build.supercells`. Write the logic to apply strains to the generated SQS cell. Run the test until it passes.
4.  **TDD for NMS:** Write a new test, `test_generate_nms_for_h2o`. This test will mock the initial geometry lookup and the Hessian calculation. It will assert that the generated structures are displaced from the original geometry and that the number of generated structures is correct.
5.  **Implement NMS:** Implement the `_generate_nms_for_molecule` method. This will involve using `scipy.linalg.eigh` to diagonalize a Hessian matrix and generate displacement vectors.
6.  **Implement Simpler Methods:** Implement the AIRSS and Deep Rattling methods. These are algorithmically simpler, and their tests will focus on ensuring the correct number of atoms are placed randomly within a cell, respecting basic constraints.
7.  **Final Integration Test:** Update the main integration test in `test_workflow.py`. The test should now start from an even more minimal `input.yaml` (without any initial structure file). It should assert that after the `StructureGenerator` runs, the database contains a set of 'unlabeled' structures, which are then correctly picked up by the (mocked) `LabelingEngine`.

## 5. Test Strategy

Testing for this cycle is crucial as it validates the scientific correctness of the initial data. Errors in structure generation can lead to fundamentally flawed potentials.

**Unit Testing Approach (Min 300 words):**
The unit tests in `test_structure_generator.py` will verify each algorithm in isolation, with a strong focus on the physical and chemical properties of the generated `ase.Atoms` objects.

*   **Testing SQS:** The `test_generate_sqs_for_fept` will not just check if atoms were added to the mock database. It will perform detailed assertions on the generated structures. It will get the list of `Atoms` objects captured by the mock `add_atoms` method. For each `Atoms` object, it will call `atoms.get_chemical_formula()` and assert it is `'FePt'` (or equivalent, like `'Fe16P16'`). It will check that the total number of atoms is within a reasonable range (e.g., `32 <= len(atoms) <= 96`). It will also test the strain application by comparing the cell vectors of the base SQS structure to the strained ones, asserting that the volume has changed as expected.
*   **Testing NMS:** The `test_generate_nms_for_h2o` will use a fixed, known geometry for a water molecule as its starting point. It will mock the calculation of normal modes, providing a fixed set of displacement vectors. The test will then assert that the generated structures have positions that are linear combinations of the original positions and these displacement vectors. It will also verify that bond lengths have changed, proving the geometry was actually distorted.
*   **Testing the Dispatcher:** A separate test will verify the `_dispatcher` logic. It will involve creating four different mock configs (one for each `structure_type`), and asserting that the correct internal generation method is called for each one, using `mocker.spy`.

**Integration Testing Approach (Min 300 words):**
The main integration test will confirm that the `StructureGenerator` correctly seeds the database for the rest of the pipeline.

*   **`test_end_to_end_from_composition`**: This test in `test_workflow.py` will simulate the full, updated workflow.
    1.  **Setup**: It will start with a minimal `input.yaml` (e.g., for FePt) and a completely empty database file in a temporary directory.
    2.  **Execution**: It will invoke the CLI runner. The mocks for `subprocess.run` (QE) and `pacemaker.fit` (training) will remain active.
    3.  **Assertions**:
        *   The test will first let the `StructureGenerator` run completely.
        *   After the generation step, it will use the real `AseDBWrapper` to connect to the database file. It will perform a query, `get_all_unlabeled_ids()`, and assert that this list is not empty, proving the generator populated the database. It will retrieve one of the generated `Atoms` objects and verify its composition, as in the unit test.
        *   The test will then allow the (mocked) `LabelingEngine` to run. It will assert that the `subprocess.run` mock was called a number of times equal to the number of structures the generator created.
        *   Finally, it will verify that the (mocked) `TrainingEngine` was called and that the data passed to it corresponds to the structures that were initially generated. This confirms the seamless flow of data from initial, automated generation all the way to the training step.
