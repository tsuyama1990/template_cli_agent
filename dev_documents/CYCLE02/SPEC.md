# Cycle 2 Specification: Structure Generation & Configuration

## 1. Summary

Cycle 2 represents a significant step forward in the MLIP-AutoPipe system's journey towards full automation. The primary goal of this cycle is to eliminate the need for manual preparation of inputs. This will be achieved by developing two key components: the **Structure Generator (Module A)** and the **Config Expander (Heuristic Engine)**. The Structure Generator is responsible for creating a diverse and physically plausible set of initial atomic structures from minimal user input. The Config Expander introduces the "Two-Tier Configuration Strategy," a cornerstone of the system's user-friendly design. It takes a very simple, user-provided `input.yaml` file and automatically expands it into a comprehensive `exec_config_dump.yaml` file, which contains all the detailed parameters required for the entire pipeline.

By the end of this cycle, the user experience will be vastly improved. Instead of needing to provide a set of atomic structures and a detailed configuration file, a user will only need to specify the chemical composition (e.g., "FePt") and a few high-level goals. The system will then autonomously generate the initial dataset and all the necessary computational settings. This cycle effectively builds the "front door" to the pipeline, making it accessible to non-experts. It embodies the core philosophy of "removing the human expert from the loop" by encapsulating expert knowledge about structure generation and DFT parameterisation into a robust, automated engine.

## 2. System Architecture

The architecture for Cycle 2 introduces the first two major components of the automated workflow, which precede the modules developed in Cycle 1.

**Config Expander (Heuristic Engine):**
This component is the brain of the initial setup phase. It is not a module in the same sense as A-E but rather a critical orchestrator that translates user intent into a precise execution plan. It will be designed as a standalone engine that is invoked at the very beginning of the pipeline.
1.  **Input:** It takes a path to a minimal `input.yaml` file.
2.  **Processing:**
    *   It parses the `input.yaml` to understand the basic system (elements, composition).
    *   It applies a set of physics-based heuristics. For example, it will analyze the elements involved to classify the material as an alloy, ionic crystal, molecule, or covalent network.
    *   Based on this classification, it determines the appropriate strategy for initial structure generation.
    *   It determines optimal DFT parameters (e.g., cutoff energies from the SSSP protocol, k-point density) and simulation settings (e.g., temperature ranges based on estimated melting points).
3.  **Output:** It produces a `exec_config_dump.yaml` file. This file is a complete, unambiguous specification of every parameter needed for the rest of the pipeline, ensuring reproducibility.

**Module A: Structure Generator (Initial Seeding):**
This module is responsible for creating the initial population of atomic structures that will form the seed for the training dataset. It is designed to be a fast, DFT-free process. The module will be driven by the configuration produced by the Config Expander.
-   **Input:** It will receive the `FullConfig` data, which specifies the structure type and other parameters.
-   **Logic:** It will contain several sub-modules, each specialised for a different type of material:
    *   **Alloys:** It will use a Special Quasirandom Structures (SQS) generation algorithm to model random alloys effectively. It will also apply various strains (volumetric and shear) to the structures to ensure the resulting potential can model elastic properties.
    *   **Molecules:** It will use Normal Mode Sampling (NMS), where the molecule is distorted along its vibrational modes to generate a diverse set of geometries.
    *   **Ionic/Covalent:** For these, it will employ techniques like random structure searching (AIRSS-like) or "deep rattling" followed by a quick relaxation with a cheap surrogate model to create a wide variety of plausible atomic configurations.
-   **Output:** It will produce a list of `ase.Atoms` objects, representing the initial, unlabeled dataset.

## 3. Design Architecture

The design for this cycle introduces a new configuration handling system and the first of the main pipeline modules.

**Configuration System:**
-   **`config/models.py`**: This file will contain two Pydantic models:
    -   `MinimalConfig`: Defines the very simple structure of `input.yaml` (e.g., `system: {elements: ["Fe", "Pt"], composition: "FePt"}`). It will have very few required fields.
    -   `FullConfig`: Defines the comprehensive structure of `exec_config_dump.yaml`. This will be a large, nested model containing detailed sections for DFT parameters, MLIP training settings, simulation parameters, etc.
-   **`config/expander.py`**: A `ConfigExpander` class will be implemented.
    -   **`__init__(self)`**: The constructor might load any necessary heuristic data (e.g., a table of elemental properties).
    -   **`expand(self, minimal_config: MinimalConfig) -> FullConfig`**: This method will contain the core logic. It will take a `MinimalConfig` object and programmatically build and return a `FullConfig` object, filling in all the required details.

**Module A Implementation:**
-   **`modules/a_structure_generator.py`**: This file will contain the `StructureGenerator` class.
    -   **`__init__(self, config: FullConfig)`**: The constructor will take the full configuration object.
    -   **`generate(self) -> list[ase.Atoms]`**: The main public method. It will read the material type from the config and call the appropriate private helper method (e.g., `_generate_sqs()`, `_generate_nms()`).
    -   **Private Methods (`_generate_sqs`, `_generate_nms`, etc.)**: Each of these methods will encapsulate the logic for a specific structure generation algorithm. They will use external libraries (e.g., `pymatgen`, `ase`, `icet`) where appropriate.

**Integration:**
The overall workflow will be that the main `cli.py` script will first instantiate and run the `ConfigExpander`. The resulting `FullConfig` will then be used to instantiate the `StructureGenerator`, which is then run to produce the initial structures. These structures would then be ready to be passed to the modules from Cycle 1.

## 4. Implementation Approach

1.  **Pydantic Models:** Define the `MinimalConfig` and `FullConfig` Pydantic models. Start with a minimal set of fields and expand them as the heuristic logic is developed.
2.  **Config Expander Logic:** Implement the `ConfigExpander` class. Begin with the material type classification logic. For a given set of elements, it should correctly identify the material type. Then, implement the logic for populating the DFT parameters based on the SSSP protocol rules.
3.  **File I/O:** Implement the code to read the `input.yaml` and write the `exec_config_dump.yaml`. Use a standard YAML library like `ruamel.yaml` to preserve comments and formatting if needed.
4.  **StructureGenerator Shell:** Create the `StructureGenerator` class and its main `generate` method.
5.  **Implement SQS Generation:** Integrate a library like `icet` or `ase.build.sqs` to implement the `_generate_sqs` method. Add the logic for applying volumetric and shear strains to the generated structures.
6.  **Implement NMS Generation:** Implement the `_generate_nms` method. This will involve using `ase` to perform a vibrational analysis on a reference molecular structure to get the normal modes, and then displacing the atoms along these modes.
7.  **Implement Other Generators:** Implement the basic logic for the ionic and covalent structure generators.
8.  **Update Integration Script:** Modify the main script to chain the new components together. It should now be able to run the entire process from an `input.yaml` file to a set of generated `ase.Atoms` objects.

## 5. Test Strategy

Testing for Cycle 2 focuses on the correctness of the heuristic logic and the output of the structure generation algorithms.

**Unit Testing Approach (Min 300 words):**
Unit tests for the `ConfigExpander` are critical. We will create a suite of tests that cover its heuristic decision-making. For example, we will create several `MinimalConfig` objects for different chemical systems (`FePt` for an alloy, `H2O` for a molecule, `NaCl` for an ionic crystal) and assert that the `expand` method produces a `FullConfig` object with the correct `structure_type` field set. Another set of tests will verify the DFT parameter generation. We can test that for an element like Silicon, the expander correctly looks up the recommended SSSP cutoff energy and places it in the `FullConfig`.

For the `StructureGenerator`, we will write unit tests for each generation method. For `_generate_sqs`, we can create a simple binary system and assert that the returned list of `ase.Atoms` objects is not empty, that each structure has the correct chemical composition, and that the cell dimensions vary if strain has been applied. For `_generate_nms`, we can test with a water molecule and verify that the generated structures show atomic displacements that are consistent with known vibrational modes (e.g., symmetric stretch, bending). These tests will verify the correctness of each algorithm in isolation.

**Integration Testing Approach (Min 300 words):**
The primary integration test for this cycle will verify the entire "front-end" workflow. The test will involve creating a temporary `input.yaml` file on disk. The test will then run a function that orchestrates the `ConfigExpander` and the `StructureGenerator`. The main assertion will be to check that this end-to-end process runs without errors and produces a non-empty list of `ase.Atoms` objects. We will inspect the properties of these generated objects to ensure they are plausible (e.g., atoms are not unreasonably close to each other). We will also check the intermediate `exec_config_dump.yaml` file to ensure it is syntactically correct and contains the expected sections. This test confirms that the `FullConfig` object produced by the expander is a valid input for the `StructureGenerator`. A second integration test will connect the output of this cycle to the input of Cycle 1. We will take the list of `Atoms` objects generated by `StructureGenerator` and pass it to a mocked `LabelingEngine` to ensure the data structures are compatible.
