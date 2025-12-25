# Specification: Cycle 2 - Initial Structure Generation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

Cycle 2 of the MLIP-AutoPipe project transitions our focus to the very beginning of the pipeline: the autonomous generation of initial atomic structures. The successful completion of Cycle 1 provided the core automated workflow for labelling and training. However, it still required the user to provide a set of input structures. Cycle 2 aims to remove this dependency, bringing us one major step closer to a truly "zero-touch" system. The primary goal of this cycle is to implement `Module A: Structure Generator`. This module will be responsible for taking the minimal user input (chemical composition) and intelligently creating a diverse, physically plausible, and computationally relevant set of initial configurations to seed the entire workflow.

The scope of this cycle is to design and build a robust system for programmatic structure creation. A key feature will be the implementation of a material type classifier. This component will analyze the provided chemical formula and, using heuristics based on fundamental chemical principles (like electronegativity), automatically classify the material as an alloy, a molecule, or an ionic/covalent crystal. This classification will then determine which specialized generation algorithm is most appropriate to use. For alloys, we will implement a Special Quasirandom Structures (SQS) generator. For molecular systems, we will use Normal Mode Sampling (NMS) to explore vibrational degrees of freedom. For ionic and covalent materials, we will implement methods like random structure searching and deep rattling to create a variety of plausible polymorphs and defect structures.

By the end of this cycle, the MLIP-AutoPipe system will be capable of initiating its workflow from nothing more than a chemical formula. The user will no longer be required to have any pre-existing structural data. This enhancement not only improves usability but also contributes to the scientific rigour of the process by ensuring that the initial training data is generated in a systematic and unbiased manner. This cycle will deliver a fully integrated `StructureGenerator` that seamlessly feeds its output into the labelling and training engines developed in Cycle 1, creating a more complete and powerful automated pipeline. The successful completion of this cycle will mean that a user can start with a single line of input, the chemical composition of their material, and the system will be able to generate a rich and diverse dataset ready for the next stages of the pipeline. This is a critical step towards realizing the project's vision of a fully autonomous materials simulation platform.

## 2. System Architecture

The architecture for Cycle 2 involves integrating `Module A: Structure Generator` at the very beginning of the existing pipeline. The workflow will now start with the user's `input.yaml`, proceed to structure generation, and then flow directly into the labelling and training modules from Cycle 1.

**Architectural Placement and Data Flow:**
`Module A` will be the first active module in the pipeline, executed immediately after the `ConfigExpander`. The overall data flow will be as follows:
1.  **Configuration:** The `ConfigExpander` runs as before, creating the `FullConfig` object. This object will now be enriched with a new `StructureGenerationConfig` section, which includes the auto-detected material type and parameters for the relevant generation algorithm.
2.  **Instantiation:** The main orchestrator will use the `StructureGenerationConfig` to instantiate the `StructureGenerator` module.
3.  **Execution:** The orchestrator will invoke the main method of the `StructureGenerator`. The generator will perform its internal logic (classification and generation) and will return a list of `ase.Atoms` objects. This list represents the initial, diverse set of structures for the workflow.
4.  **Handoff:** This list of `ase.Atoms` objects is then passed directly to `Module C: Labeling Engine` as its primary input.

From this point onwards, the pipeline proceeds exactly as it did in Cycle 1. The Labeling Engine will calculate DFT properties for the newly generated structures, and the Training Engine will use this data to train the MLIP.

**Internal Architecture of `Module A: Structure Generator`:**
The module itself will be designed with a clear, hierarchical structure:
*   **A Top-Level `StructureGenerator` Class:** This class will serve as the public interface to the module. It will have a single primary method, `generate()`, which encapsulates the entire logic.
*   **A `MaterialClassifier` Component:** Internally, the `generate()` method will first call upon a dedicated `MaterialClassifier`. This component will take the elemental composition and use a rules-based engine to determine the material type. This decision-making logic will be isolated to make it easy to refine the heuristics in the future.
*   **Specialized Generation Sub-modules:** Based on the classifier's output, the main class will delegate the actual generation task to one of several specialized sub-modules or functions, each responsible for a single generation algorithm (e.g., `_generate_sqs_for_alloy`, `_generate_nms_for_molecule`).

This internal architecture makes the module highly extensible. Adding a new material type or a new generation algorithm in the future would simply involve adding a new classification rule and a corresponding private generation method, without altering the public interface or the overall system workflow. The module will be stateless; its output will depend only on its configuration, ensuring that the entire process remains reproducible. The clear separation of concerns also allows for easier testing and maintenance of the individual components. The `StructureGenerator` will be designed to be robust to a wide range of inputs, and will include error handling for cases where the user provides an invalid or unsupported composition.

## 3. Design Architecture

The design of `Module A` will be integrated into the existing `src/mlip_autoflow` package. We will introduce a new file for this module and a supporting file for the classification logic.

**Updated Project Structure:**
```
src/mlip_autoflow/
├── __init__.py
├── main.py
├── config/
│   ├── expander.py
│   └── models.py
├── data/
│   └── ase_db_manager.py
└── modules/
    ├── a_structure_generator.py
    ├── c_labeling_engine.py
    └── d_training_engine.py
    └── utils/
        └── material_classifier.py # New file
```

**Class and Method Definitions:**

*   **`config.models.py`**: The Pydantic models will be updated.
    *   A new `StructureGenerationConfig` model will be added to `FullConfig`. This will include fields like `material_type: str` (read-only, set by the expander), `sqs_config: SQSConfig`, `nms_config: NMSConfig`, etc.
    *   The `UserConfig` will be updated to accept structure information, like a CIF file path, as an optional alternative to just a composition string.

*   **`modules.utils.material_classifier.py`**: This new file will contain the classification logic.
    *   A `MaterialClassifier` class.
    *   `__init__(self, elements: List[str])`: Takes the list of elements.
    *   `classify(self) -> str`: This method will contain the heuristic logic. It will use libraries like `pymatgen` or `mendeleev` to access element properties (e.g., electronegativity, atomic radii) and apply a set of `if/elif/else` rules to return a string: `"alloy"`, `"molecule"`, `"ionic"`, or `"covalent"`.

*   **`modules.a_structure_generator.py`**: This file will house the main `StructureGenerator` class.
    *   `__init__(self, config: StructureGenerationConfig, composition: dict)`: Initializes with its configuration section and the system composition.
    *   `generate(self) -> List[ase.Atoms]`: The main public method. It will first call the `MaterialClassifier`, then call the appropriate internal generation method based on the result.
    *   `_generate_sqs_for_alloy(self) -> List[ase.Atoms]`: A private method that will interface with a library like `icet` to generate Special Quasirandom Structures. It will also apply various strains (volumetric and shear) to the generated structures to create a more diverse dataset.
    *   `_generate_nms_for_molecule(self) -> List[ase.Atoms]`: A private method that will perform Normal Mode Sampling. This involves first performing a geometry optimization with a cheap potential (or using a provided initial structure), calculating the Hessian, and then displacing the atoms along the normal modes.
    *   `_generate_for_ionic_covalent(self) -> List[ase.Atoms]`: A private method that will implement techniques like Ab Initio Random Structure Searching (AIRSS) and "deep rattling" to generate a variety of crystal structures and defect configurations.

*   **`main.py`**: The orchestration logic in the main CLI script will be updated to include this new step. The sequence will now be: `ConfigExpander` -> `StructureGenerator` -> `LabelingEngine` -> `TrainingEngine`.

This design cleanly separates the logic for "what type of material is it?" (`MaterialClassifier`) from "how do I make structures for it?" (`StructureGenerator`), promoting modularity and testability.

## 4. Implementation Approach

The implementation of Cycle 2 will be a step-by-step process, starting with the foundational classification logic and then building out each of the specialized generation methods.

**Step 1: Update Configuration Models**
The first task is to update the Pydantic models in `config/models.py`. We will add the new `StructureGenerationConfig` and its sub-models. This defines the data contract for the new module and ensures that all necessary parameters are passed through the system correctly. The `ConfigExpander` will be updated to populate these new configuration fields, including the logic to call the `MaterialClassifier`.

**Step 2: Implement the Material Classifier**
We will create the `MaterialClassifier` class in `modules/utils/material_classifier.py`. The initial implementation will use simple, robust heuristics. For example:
*   If all elements are metals, classify as `"alloy"`.
*   If electronegativity differences are large, classify as `"ionic"`.
*   If the input is a recognized molecular formula (e.g., H2O), classify as `"molecule"`.
We will add dependencies like `pymatgen` or `mendeleev` to `pyproject.toml` to provide the necessary elemental data. We will write unit tests to verify that the classifier correctly identifies a range of example compositions.

**Step 3: Implement SQS Generation for Alloys**
Next, we will tackle the first generation method in `StructureGenerator`. We will implement the `_generate_sqs_for_alloy` method. This will involve adding a dependency on a suitable SQS generation library, such as `icet`. The implementation will involve:
1.  Setting up the primitive lattice for the alloy system.
2.  Calling the SQS generation function from the library to get the best structure.
3.  Creating copies of the structure and applying a series of deterministic volumetric and shear strains to it.
4.  Returning the collection of original and strained structures as a list of `ase.Atoms` objects.
Unit tests will be written to ensure that this method returns a list of valid, strained structures with the correct composition.

**Step 4: Implement Normal Mode Sampling for Molecules**
After the alloy generator, we will implement `_generate_nms_for_molecule`. This is a more complex multi-step process:
1.  Ensure a single, optimized molecular structure is available. If the user provides a structure file, use that; otherwise, perform a quick geometry optimization using a universal force field.
2.  Calculate the Hessian matrix (second derivatives of the energy). This can be done numerically using finite differences.
3.  Diagonalize the Hessian to obtain the normal modes (eigenvectors) and their frequencies (eigenvalues).
4.  Generate new structures by displacing the atoms along each normal mode by various amplitudes.
This method will be tested to ensure it produces a set of structures that correctly represent vibrations around the equilibrium geometry.

**Step 5: Implement Generators for Ionic/Covalent Systems**
The final generation method will be for ionic and covalent systems. We will implement `_generate_for_ionic_covalent`. This will initially focus on a "deep rattling" approach, where atoms in a known crystal structure are displaced by a large random amount, followed by a local geometry optimization with a cheap potential. This is effective at finding nearby local minima. We can also implement a simplified random structure search, where atoms are placed randomly in a box (respecting ionic radii) and then relaxed.

**Step 6: Integrate into the Main Workflow**
With all the generation methods implemented, the final step is to update the main orchestrator in `main.py`. The script will be modified to instantiate and run the `StructureGenerator` after the `ConfigExpander`. The list of `ase.Atoms` objects returned by the generator will then be passed directly to the `LabelingEngine`. The CLI will be updated so the user can now run the pipeline by providing a composition string instead of a structure file.

## 5. Test Strategy

The testing for Cycle 2 will focus heavily on `Module A`, ensuring that each component of the structure generation process is reliable and produces scientifically valid outputs.

**Unit Testing Approach:**
*   **`MaterialClassifier`**: This will be tested extensively. We will create a test suite with a variety of chemical compositions and assert that the `classify()` method returns the correct material type. For example, `{"Fe": 0.5, "Pt": 0.5}` should return `"alloy"`, `{"H": 2, "O": 1}` should return `"molecule"`, and `{"Na": 1, "Cl": 1}` should return `"ionic"`.
*   **`StructureGenerator`**: Each of the private generation methods (`_generate_sqs_for_alloy`, etc.) will be unit tested.
    *   The SQS test will call the method and assert that the output is a list of `ase.Atoms` objects, that all objects have the correct composition, and that the cell structures correspond to the applied strains.
    *   The NMS test will provide a simple molecule (e.g., H2O) and verify that the generated structures represent displacements along the known vibrational modes. We can check this by projecting the displacements onto the true normal modes.
    *   The ionic/covalent generator test will verify that the output structures are chemically plausible (e.g., atoms are not unrealistically close to each other) and have the correct composition.
*   We will also test the main `generate()` method by mocking the classifier. For example, we will force the classifier to return `"alloy"` and verify that the `_generate_sqs_for_alloy` method is called.

**Integration Testing Approach:**
The primary integration test for Cycle 2 will verify the seamless flow of data from the new `StructureGenerator` into the existing `LabelingEngine`.
*   **Test Scenario:** The test will cover the initial part of the full pipeline.
    1.  **Setup:** Create a minimal `input.yaml` file containing only a simple composition, for example, `composition: "Si2"`.
    2.  **Execution:**
        *   Run the `ConfigExpander` on this input file.
        *   Instantiate and run the `StructureGenerator` using the expanded configuration.
        *   Take the list of `ase.Atoms` objects produced by the generator and pass it to the `LabelingEngine` (which will again be configured to use a "dummy" DFT wrapper script).
    3.  **Verification:**
        *   Assert that the `StructureGenerator` produces a non-empty list of `ase.Atoms` objects.
        *   Verify that all generated structures have the correct composition ("Si2").
        *   Check that the `LabelingEngine` runs without errors and that the dummy DFT script was called for each of the generated structures.
        *   Inspect the temporary database to ensure that it now contains entries corresponding to the structures created by `Module A`.

This integration test will provide high confidence that the new structure generation module is correctly "wired" into the pipeline and that the data it produces is in the correct format for consumption by the downstream modules. This lays the groundwork for a fully autonomous workflow.
