# User Acceptance Test (UAT): Cycle 2

**Version:** 1.0.0
**Status:** Ready for Testing

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 2. The focus is on verifying that `Module A: Structure Generator` can correctly interpret a user's minimal input, classify the material type, and generate a scientifically valid and diverse set of initial atomic structures.

| Scenario ID | Test Scenario Description                                                                | Priority |
|-------------|------------------------------------------------------------------------------------------|----------|
| UAT-C2-001  | **Verify SQS Generation for an Alloy:** The user provides a composition for a binary alloy (e.g., FePt). The system must identify it as an alloy and generate a set of Special Quasirandom Structures (SQS), including strained variations. | High     |
| UAT-C2-002  | **Verify NMS Generation for a Molecule:** The user provides a composition for a simple molecule (e.g., H2O). The system must identify it as a molecule and generate a set of structures corresponding to displacements along its normal vibrational modes. | High     |
| UAT-C2-003  | **Verify Random Structure Generation for an Ionic Solid:** The user provides a composition for a simple ionic solid (e.g., NaCl). The system must identify the material type and generate a set of plausible, varied crystal structures. | Medium   |
| UAT-C2-004  | **Verify Correct Material Classification:** The user provides several different compositions (e.g., "SiGe", "CH4", "LiF"). The user inspects the logs or configuration output to confirm that the system correctly classifies each one as "alloy", "molecule", and "ionic" respectively. | High     |
| UAT-C2-005  | **Verify Pipeline Integration:** The user runs the full pipeline (Cycles 1 & 2) with a simple composition string. The system must run end-to-end, from generating structures to training a final MLIP, without crashing. | High     |


## 2. Behavior Definitions

This section provides detailed Gherkin-style behavior definitions for the test scenarios, outlining the exact steps and expected outcomes.

---

### **Scenario: UAT-C2-001 - Verify SQS Generation for an Alloy**

*   **GIVEN** a minimal `input.yaml` file with `composition: "FePt"` and `elements: ["Fe", "Pt"]`.
*   **WHEN** the user runs the `StructureGenerator` module.
*   **THEN** the system should produce a list or file containing multiple `ase.Atoms` objects.
*   **AND** all generated structures must have a 50/50 composition of Iron and Platinum atoms.
*   **AND** the generated structures should include one unstrained SQS structure and several others with various volumetric and shear strains applied.
*   **AND** no two atoms in any generated structure should be closer than a physically realistic distance (e.g., ~1.5 Angstroms).
*   **AND** the log output should clearly indicate that the material was classified as an alloy and that the SQS generation method was used.

---

### **Scenario: UAT-C2-002 - Verify NMS Generation for a Molecule**

*   **GIVEN** a minimal `input.yaml` file with `composition: "H2O"` and `elements: ["H", "O"]`. An initial structure of a water molecule is also provided.
*   **WHEN** the user runs the `StructureGenerator` module.
*   **THEN** the system should produce a set of `ase.Atoms` objects.
*   **AND** one of the structures should be the original, optimized geometry of the water molecule.
*   **AND** the other structures should show atomic displacements that visually correspond to the known symmetric stretch, asymmetric stretch, and bending modes of a water molecule.
*   **AND** all structures must contain exactly two Hydrogen atoms and one Oxygen atom.
*   **AND** the log output should clearly indicate that the material was classified as a molecule and that the NMS generation method was used.

---

### **Scenario: UAT-C2-003 - Verify Random Structure Generation for an Ionic Solid**

*   **GIVEN** a minimal `input.yaml` file with `composition: "NaCl"` and `elements: ["Na", "Cl"]`.
*   **WHEN** the user runs the `StructureGenerator` module.
*   **THEN** the system should produce a set of varied crystalline structures.
*   **AND** each structure must contain an equal number of Sodium and Chlorine atoms.
*   **AND** the structures should represent different plausible crystal packings (polymorphs).
*   **AND** a visual inspection of the structures should show ordered, crystal-like arrangements, not amorphous or random clusters.
*   **AND** the log output should clearly indicate that the material was classified as ionic and that the appropriate random structure generation method was used.

---

### **Scenario: UAT-C2-004 - Verify Correct Material Classification**

*   **GIVEN** an `input.yaml` file with `composition: "SiGe"`.
*   **WHEN** the user runs the `ConfigExpander` and `StructureGenerator` initialization.
*   **THEN** the `exec_config_dump.yaml` or log output must explicitly state the classified material type is `"alloy"`.
*   **GIVEN** an `input.yaml` file with `composition: "CH4"`.
*   **WHEN** the user runs the `ConfigExpander` and `StructureGenerator` initialization.
*   **THEN** the `exec_config_dump.yaml` or log output must explicitly state the classified material type is `"molecule"`.
*   **GIVEN** an `input.yaml` file with `composition: "LiF"`.
*   **WHEN** the user runs the `ConfigExpander` and `StructureGenerator` initialization.
*   **THEN** the `exec_config_dump.yaml` or log output must explicitly state the classified material type is `"ionic"`.
*   **AND** in each case, the subsequent structure generation method used should be the one appropriate for the classified material type.

---

### **Scenario: UAT-C2-005 - Verify Pipeline Integration**

*   **GIVEN** a directory containing only a minimal `input.yaml` with `composition: "Si2"`.
*   **WHEN** the user executes the full pipeline command (e.g., `mlip-pipe run input.yaml`).
*   **THEN** the system should execute the entire workflow (Structure Generation, Labeling, Training) without any crashes or fatal errors.
*   **AND** the log files should show that `Module A` generated a set of structures, which were then passed to `Module C`.
*   **AND** the final output directory must contain a valid, trained MLIP file (e.g., `final_model.ace`).
*   **AND** the final database must contain entries for the structures that were programmatically generated by the system.
*   **AND** the number of entries in the database should match the number of structures generated by `Module A`.
*   **AND** the entire process should be completed in a reasonable amount of time.
