# Cycle 02 User Acceptance Test (UAT): Structure Generation and Configuration

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 02. The focus is on verifying the new user-facing capabilities of the pipeline from an end-user's perspective. These tests are designed to confirm the successful implementation of two major features: the ability to start a workflow from a minimal, human-friendly configuration file, and the system's capacity to automatically generate a relevant and diverse set of initial atomic structures without manual intervention. Passing these tests will demonstrate that the system has become significantly more autonomous and user-friendly, fulfilling the core goals of this development cycle.

| Scenario ID | Description                                                              | Priority |
|-------------|--------------------------------------------------------------------------|----------|
| UAT-02-001  | **Successful Configuration Expansion:** Verify that the system can take a minimal `input.yaml` and correctly expand it into a complete, well-formed, and physically sensible `exec_config_dump.yaml`. This test validates the core logic of the heuristic engine and its ability to correctly set up a complex calculation from simple inputs. | High     |
| UAT-02-002  | **Automatic Structure Generation for an Alloy:** Verify that, given an input for a binary alloy (e.g., FePt), the `StructureGenerator` correctly identifies it as an alloy and generates a set of Special Quasirandom Structures (SQS), including strained variations. This test confirms the context-specific behaviour of the generation module for a major class of materials. | High     |
| UAT-02-003  | **Automatic Structure Generation for a Molecule:** Verify that, given an input for a molecule (e.g., H2O), the `StructureGenerator` correctly identifies it as a molecule and generates a set of physically plausible distorted structures by sampling its normal modes. This test confirms the module's correct behaviour for a different, chemically distinct class of materials. | High     |
| UAT-02-004  | **End-to-End Run from Minimal Input:** Verify that the entire pipeline (generation, labelling, training) can be successfully and automatically executed starting *only* from a minimal `input.yaml` for a simple, well-behaved system. This is the capstone test for this cycle, proving that all new components are correctly integrated into the main workflow. | High     |
| UAT-02-005  | **Error on Invalid `input.yaml`:** Verify that the system provides a clear, informative, and user-friendly error message if the user provides a malformed or incomplete `input.yaml` file. This test ensures the system is robust against common user errors and provides helpful feedback. | Medium   |

## 2. Behaviour Definitions

The following Gherkin-style definitions describe the expected behaviour for each test scenario in detail.

---

### **Scenario: UAT-02-001 - Successful Configuration Expansion**

This scenario tests the core functionality of the `ConfigExpander` heuristic engine. The user experience of the entire application depends on this component's ability to correctly translate a simple, high-level request into a detailed, technically correct set of instructions for the pipeline's modules. It must not only create a complete file but also make intelligent, physically-based decisions on behalf of the user.

**GIVEN** a user has created a minimal `input.yaml` file containing only the following lines:
```yaml
system:
  elements: ["Fe", "Pt"]
  composition: "FePt"
```
**WHEN** the user executes the main orchestrator script, pointing to this input file.
**THEN** the system should, as its very first action, create a new file named `exec_config_dump.yaml` in the output directory.
**AND** this file should be a syntactically valid YAML file.
**AND** the file should contain all the necessary top-level keys for a complete run: `system`, `simulation`, `dft_compute`, and `mlip_training`.
**AND** under the `system` key, the `structure_type` parameter should be automatically and correctly inferred and set to the string `"alloy"`, based on the metallic nature of Fe and Pt.
**AND** under the `dft_compute` key, crucial parameters like `pseudopotentials`, `ecutwfc`, and `kpoints_density` should be populated with sensible default values derived from the system's internal SSSP protocol data for Iron and Platinum.
**AND** the `magnetism` setting should be appropriately set, for example to `"ferromagnetic"`, as the heuristic engine should recognise that Iron is a magnetic element.

---

### **Scenario: UAT-02-002 - Automatic Structure Generation for an Alloy**

This scenario verifies that the `StructureGenerator` correctly implements the appropriate physical model for generating disordered alloy structures. For an MLIP to be effective for an alloy, it must be trained on structures that correctly capture the random distribution of atomic species, as well as the effect of lattice strain on the local atomic environments. SQS is the standard method for achieving this in a computationally efficient manner.

**GIVEN** an `input.yaml` file has been created for the composition "SiGe".
**AND** a clean and empty ASE database has been initialised.
**WHEN** the user runs the pipeline, and the orchestrator executes the `StructureGenerator` module.
**THEN** the system should log a clear, human-readable message indicating that it has classified the material as an alloy and is now generating SQS structures.
**AND** upon completion of the module, the ASE database should contain multiple new entries (e.g., more than 10 structures).
**AND** each of these new entries in the database should have been assigned the state "unlabelled".
**AND** a user, by inspecting the generated structures (e.g., by exporting them to a visualisable format like XYZ), should be able to confirm that they represent varied atomic arrangements of Silicon and Germanium atoms on a diamond-like lattice, consistent with the SQS methodology.
**AND** the user should also confirm that some of the generated structures have had different strains applied, resulting in lattice vectors that are slightly different from a standard equilibrium cell.

---

### **Scenario: UAT-02-003 - Automatic Structure Generation for a Molecule**

This scenario tests the generator's ability to handle a completely different type of material: a molecule. The physics of molecular interactions is dominated by the stretching and bending of covalent bonds. Therefore, the training data must include geometries that sample these vibrational modes. Normal Mode Sampling is a physically motivated way to achieve this.

**GIVEN** an `input.yaml` file has been created for the composition "H2O".
**AND** a clean and empty ASE database has been initialised.
**WHEN** the user runs the pipeline, and the orchestrator executes the `StructureGenerator` module.
**THEN** the system should log a message indicating that it has correctly classified the system as a molecule and is now generating structures via Normal Mode Sampling (NMS).
**AND** upon completion of the module, the ASE database should contain multiple new entries.
**AND** a user, by inspecting the generated structures, should observe distorted versions of the water molecule. The O-H bond lengths and H-O-H bond angles should be slightly perturbed from their known equilibrium values, in patterns consistent with the symmetric stretch, asymmetric stretch, and bending modes of a water molecule.

---

### **Scenario: UAT-02-004 - End-to-End Run from Minimal Input**

This is the capstone test for Cycle 02, demonstrating that all the new components work together and are correctly integrated with the core engine from Cycle 01. It proves that the system can now perform a complete, scientifically meaningful task starting from nothing more than a high-level user request.

**GIVEN** a minimal `input.yaml` for a simple elemental system like bulk Silicon (`elements: ["Si"]`).
**AND** an empty output directory and a clean ASE database.
**WHEN** the user executes the main orchestrator script with this single input file.
**THEN** the system should execute the full, four-stage workflow without crashing or requiring any further user interaction.
**AND** the console log should show the sequential execution of the main stages in the correct order: 1. Configuration Expansion, 2. Structure Generation, 3. DFT Labelling, 4. MLIP Training.
**AND** at the end of the process, a final MLIP model file (e.g., `model.ace`) should exist in the output directory and be non-empty.
**AND** the ASE database should contain a set of structures that are all in the "labelled" state, indicating that the entire initial dataset was successfully processed.

---

### **Scenario: UAT-02-005 - Error on Invalid `input.yaml`**

A robust tool must be able to handle user error gracefully. This scenario tests the validation layer of the configuration system. The use of Pydantic models for configuration validation should allow the system to catch errors early and provide specific, helpful feedback, rather than failing with an obscure internal traceback.

**GIVEN** a user has created an `input.yaml` file where a required key is mis-spelled (e.g., `elementss` instead of `elements`).
**WHEN** the user executes the main orchestrator script with this invalid file.
**THEN** the system should terminate gracefully within the first few seconds, before starting any long-running computational tasks like structure generation or DFT.
**AND** the system should print a user-friendly and specific error message to the console, for instance, "ERROR: Invalid configuration in `input.yaml`. Field `elementss` is not a recognised field."
**AND** no `exec_config_dump.yaml` file should be created.
**AND** the ASE database should remain empty, and no temporary directories should be created.