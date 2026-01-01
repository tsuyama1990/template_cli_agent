# User Acceptance Testing (UAT): CYCLE03 - Initial Structure Generation

## 1. Test Scenarios

This UAT focuses on the system's ability to autonomously generate its own starting data. The user, a materials scientist, will test the `StructureGenerator` module by providing only a chemical formula. They will then verify that the generated structures are physically reasonable and correctly stored in the database, ready for the next stages of the pipeline.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C03-01  | Successful generation of SQS structures for a binary alloy        | High     |
| UAT-C03-02  | Verification of generated structures in the database              | High     |
| UAT-C03-03  | Successful generation of distorted molecules via NMS              | Medium   |

---

### **Scenario UAT-C03-01: Successful generation of SQS structures for a binary alloy**

**(Min 300 words)**
This scenario tests the most common use-case for the structure generator: creating initial configurations for a metallic alloy. The user wants to start a workflow for the FePt system but has no initial structure files. They will rely entirely on the system's SQS (Special Quasirandom Structures) generation capability.

The user will create a minimal `input.yaml` containing just `system: {elements: [Fe, Pt], composition: FePt}`. They will then run the pipeline. The primary acceptance criterion is the creation of a set of valid and diverse atomic structure files or database entries. The user is not expected to have the `LabelingEngine` run in this UAT. They can stop the process after the generation step. They will then connect to the `asedb.db` file using the `ase db` command-line tool or a custom script. They will query the database to see how many structures were generated. They should expect a reasonable number, for example, 10 to 20 initial structures.

The user will then extract one or two of these structures and save them to a file format like CIF or POSCAR for visualization. By opening the file in a visualization software (like VESTA or OVITO), they will visually inspect the structure to confirm that it represents a plausible, disordered FePt alloy. They will check the chemical composition to ensure a 50/50 ratio of Fe and Pt atoms. Furthermore, they will check the cell parameters of several of the generated structures and verify that some of them are slightly different, confirming that volumetric and shear strains have been applied to create a richer dataset.

---

### **Scenario UAT-C03-02: Verification of generated structures in the database**

**(Min 300 words)**
This scenario focuses on the data integrity aspect of the structure generation process. It's not enough to just create the structures; they must be correctly registered in the project's database to be picked up by the rest of the pipeline. The user will act as a data steward, carefully checking the database state after the generator has run.

This test follows directly from UAT-C03-01. After running the pipeline on the FePt input, the user will use a database browser (e.g., `sqlite3`) to directly inspect the contents of `asedb.db`. The first acceptance criterion is to count the number of rows in the main table and confirm it matches the number of structures the system reported generating. For each row, the user will check the key-value pairs. Crucially, they must verify that every single new row has a key named `state` with the value `'unlabeled'`. This is critical, as it flags these structures for processing by the `LabelingEngine` in the next step of the full workflow.

The user will also perform a spot-check on the data itself. They will retrieve the `atoms` object for a few rows and verify that the data stored in the database (atomic numbers, positions, cell) is correct. They will also check for the absence of other state-related keys. For example, there should be no `dft_result` key associated with these new rows. This confirms that the `StructureGenerator` is correctly populating the database with pristine, ready-to-be-labeled structures, ensuring a clean handoff to the next module in the pipeline.

---

### **Scenario UAT-C03-03: Successful generation of distorted molecules via NMS**

**(Min 300 words)**
This scenario tests the generator's versatility, ensuring it can handle non-crystalline, molecular systems in addition to alloys. The user will test the Normal Mode Sampling (NMS) functionality using a simple molecule like water (H2O).

The user will create an `input.yaml` file with the content: `system: {elements: [H, O], composition: H2O}`. As the heuristic engine should identify this as a 'molecule', this will trigger the NMS generation path. After running the generation step, the user will again inspect the database. They will expect to find a series of new, `'unlabeled'` entries. The core of the UAT is to verify the geometries of these generated structures.

The user will extract all the generated H2O structures from the database. One of these should represent the approximate equilibrium geometry. The user will measure its bond lengths and angle to confirm they are close to the known experimental values for water (~0.96 Å and ~104.5°). Then, they will inspect the other structures. The acceptance criterion is that these other structures are clearly distorted versions of the equilibrium geometry. The user will measure their bond lengths and angles and confirm that they deviate from the equilibrium values. Some structures should show stretched O-H bonds, others compressed bonds, and others a change in the H-O-H angle. This variation is the expected output of NMS and confirms the system is correctly exploring the potential energy surface around the minimum.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C03-01 & UAT-C03-02 - Successful generation and database verification for a binary alloy**

*   **GIVEN** a clean project environment with a fully configured pipeline.
*   **AND** the user has created a minimal `input.yaml` for FePt: `system: {elements: [Fe, Pt], composition: FePt}`.
*   **AND** the project database `asedb.db` is completely empty.
*   **WHEN** the user executes the command `cdd run input.yaml`.
*   **AND** the `ConfigExpander` correctly identifies the system as an `"alloy"`.
*   **THEN** the `WorkflowOrchestrator` should call the `StructureGenerator` module.
*   **AND** the `StructureGenerator` should execute its SQS algorithm.
*   **AND** the system should log a message indicating that it is generating an initial set of structures.
*   **AND** after the process finishes, the user should be able to connect to the `asedb.db` file.
*   **AND** the database should now contain a non-zero number of entries (e.g., at least 10).
*   **AND** every entry in the database created by the generator must have a key-value pair of `'state': 'unlabeled'`.
*   **AND** when an entry is retrieved as an `ase.Atoms` object, its chemical formula must be equivalent to 'FePt'.
*   **AND** at least one of the generated structures should have a cell matrix that is different from the others, indicating that strain was applied.

---

**Scenario: UAT-C03-03 - Successful generation of distorted molecules via NMS**

*   **GIVEN** a clean project environment.
*   **AND** the user has created a minimal `input.yaml` for H2O: `system: {elements: [H, O], composition: H2O}`.
*   **AND** the project database `asedb.db` is empty.
*   **WHEN** the user executes the command `cdd run input.yaml`.
*   **AND** the `ConfigExpander` correctly identifies the system as a `"molecule"`.
*   **THEN** the `WorkflowOrchestrator` should call the `StructureGenerator` module.
*   **AND** the `StructureGenerator` should execute its NMS algorithm.
*   **AND** the system should log that it is generating structures using Normal Mode Sampling.
*   **AND** the database `asedb.db` should be populated with a set of new, `'unlabeled'` structures.
*   **AND** when the geometries of these structures are analyzed, one of them should have bond lengths and angles close to the equilibrium values for water.
*   **AND** the other generated structures should show a distribution of bond lengths and angles that deviate from the equilibrium values, confirming that distortions have been successfully generated.
