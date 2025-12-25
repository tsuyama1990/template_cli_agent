# Cycle 2 User Acceptance Tests: Structure Generation and Configuration System

**Version:** 1.0.0
**Status:** Final

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 2. The focus of this cycle is on validating the new user-centric features that represent the first major step towards true "fire-and-forget" automation. These tests are designed to ensure that the two primary deliverables of the cycle—the **Two-Tier Configuration System** and the **automated Structure Generator (Module A)**—not only function correctly but also provide a seamless and intuitive user experience. The goal is to verify, from an end-user's perspective, that the system can now be successfully operated with minimal, high-level input, and that the complex, expert-level decisions it makes under the hood are both correct and transparent. Passing these UATs will confirm that we have successfully removed the user from the tedious and error-prone loop of initial data preparation.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C2-001** | **Successful Run from Minimal Input:** This is the primary "golden path" scenario for this cycle. It tests the entire integrated workflow, starting from the simplest possible user input for a common material class (an alloy). Its success demonstrates that the Config Expander and Structure Generator are working together and are correctly integrated with the core engine from Cycle 1. This test is the ultimate validation that the system can now operate with a significantly higher level of autonomy. | **High** |
| **UAT-C2-002** | **Correct Heuristic-Based Configuration:** This scenario directly tests the "intelligence" of the `ConfigExpander`. The system's ability to make sound, expert-level decisions from minimal input is a core value proposition. This test ensures that the heuristics are not only functioning but are also correct from a scientific standpoint. By checking the generated `exec_config_dump.yaml`, the user can verify that the system correctly classified the material and chose appropriate, high-quality settings, building trust in the automation. | **High** |
| **UAT-C2-003** | **Material-Specific Structure Generation:** A key feature of the `StructureGenerator` is its versatility. It is not a one-size-fits-all tool but a collection of specialized algorithms tailored to different classes of materials. This test verifies that the `ConfigExpander`'s material classification correctly triggers the appropriate algorithm in the `StructureGenerator`. It confirms that the system applies the correct physical heuristics for each material type (e.g., SQS for alloys, NMS for molecules), which is crucial for generating a relevant and physically meaningful initial dataset. | **High** |
| **UAT-C2-004** | **Provenance of Generated Configuration:** Scientific reproducibility is a non-negotiable requirement. Since the system is now making many decisions automatically, it is essential that these decisions are recorded. This scenario tests the system's commitment to provenance. It verifies that the complete, unabridged `exec_config_dump.yaml` is saved as a permanent artifact of the run. This file acts as the ultimate audit trail, ensuring that the user can always go back and understand the exact parameters that were used, and it makes the entire calculation fully reproducible. | **Medium** |
| **UAT-C2-005** | **User Input Validation:** As the system's user interface is now the `input.yaml` file, it is critical that it handles user errors gracefully. A user will inevitably make typos or forget required fields. This scenario tests the system's robustness to such errors. It ensures that the Pydantic-based validation provides immediate, clear, and helpful feedback, guiding the user to correct their input file without allowing the pipeline to proceed with a faulty configuration and fail later in a more confusing way. | **Medium** |

---

## 2. Behavior Definitions

### **UAT-C2-001: Successful Run from Minimal Input**

**Scenario:** A materials science researcher wants to generate a new MLIP for a binary alloy, FePt. They know the composition but do not want to spend time creating dozens of slightly different atomic arrangements or looking up the optimal DFT settings. They want the pipeline to do all the heavy lifting.
> **GIVEN** a valid `input.yaml` file that is intentionally minimal, containing only the following lines:
> ```yaml
> system:
>   composition: "FePt"
> ```
> **AND** the user has a fully working installation of the MLIP-AutoPipe pipeline, including a configured Quantum Espresso executable.
>
> **WHEN** the user executes the main pipeline command, providing only the path to this minimal `input.yaml` file.
>
> **THEN** the system must execute the entire, multi-stage workflow (configuration expansion, structure generation, DFT labelling, and MLIP training) without crashing or requiring any further user interaction.
> **AND** the initial log output must clearly indicate that it is parsing the `input.yaml` and expanding it into a full configuration.
> **AND** subsequently, the log must show that the system is entering the "Structure Generation" phase and is creating a new set of initial structures before it proceeds to the labelling phase.
> **AND** the process must complete successfully with a zero exit code, and a final trained MLIP model file (e.g., `fept.model`) must be present in the designated output directory.

### **UAT-C2-002: Correct Heuristic-Based Configuration**

**Scenario:** A user is working with a covalent material, Silicon Carbide (SiC), and wants to verify that the system's automated decision-making aligns with their own expert knowledge of how such a system should be treated.
> **GIVEN** an `input.yaml` file containing only the composition of Silicon Carbide:
> ```yaml
> system:
>   composition: "SiC"
> ```
> **WHEN** the user runs the pipeline, which then generates and saves an `exec_config_dump.yaml` file in the output directory.
>
> **THEN** the user must be able to open and inspect this `exec_config_dump.yaml` file.
> **AND** within this file, the `system` section must contain a key `structure_type` with the correctly identified value `covalent`.
> **AND** the `dft_compute` section must be fully populated with the correct SSSP Precision pseudopotentials and energy cutoffs for both `Si` and `C`, demonstrating correct parameter retrieval.
> **AND** the `structure_generation` section must specify a strategy that is physically appropriate for a rigid, covalent material, such as `deep_rattling` or `melt_quench`, confirming that the heuristic engine made a sound scientific choice for the generation algorithm.

### **UAT-C2-003: Material-Specific Structure Generation**

**Scenario:** A researcher is working on several different projects and wants to be confident that the pipeline is not just a "black box" but is applying the correct, distinct, and physically-appropriate heuristic for generating initial structures for each of their unique material systems.
> **GIVEN** the user has prepared a set of four different, minimal `input.yaml` files, one for each material class:
> 1.  `alloy.yaml` containing `composition: "CuAu"`
> 2.  `molecule.yaml` containing `composition: "H2O"`
> 3.  `ionic.yaml` containing `composition: "NaCl"`
> 4.  `covalent.yaml` containing `composition: "C"`
>
> **WHEN** the user runs the pipeline four separate times, once for each of these input files.
>
> **THEN** the log output produced during the run of `alloy.yaml` must contain a clear and explicit message stating that it is using the **SQS (Special Quasirandom Structures)** method.
> **AND** the log output for the `molecule.yaml` run must explicitly state that it is using the **NMS (Normal Mode Sampling)** method.
> **AND** the log output for the `ionic.yaml` run must explicitly state that it is using the **AIRSS (Ab Initio Random Structure Searching)** method.
> **AND** the log output for the `covalent.yaml` run must explicitly state that it is using a method like **Deep Rattling** or **Melt-Quench**. This confirms that the internal logic correctly dispatches to the appropriate physical model for each case.

### **UAT-C2-004: Provenance of Generated Configuration**

**Scenario:** A user has completed a successful run that was initiated from a very simple `input.yaml`. For inclusion in a publication, they need to report the exact, detailed parameters that were automatically chosen and used by the system.
> **GIVEN** a pipeline run was initiated with a minimal `input.yaml` and has completed successfully, creating an output directory.
>
> **WHEN** the user navigates into and inspects the contents of this output directory.
>
> **THEN** the directory must contain the expected final products: the trained model file and the ASE database.
> **AND** crucially, it must also contain a file named exactly `exec_config_dump.yaml`.
> **AND** when the user opens this YAML file, its contents must be a complete, explicit, and verbose configuration. It should contain no implicit defaults.
> **AND** the user must be able to take this `exec_config_dump.yaml`, use it as an input to a new pipeline run, and have it reproduce the exact same calculation, thus confirming its role as a perfect, reproducible record.

### **UAT-C2-005: User Input Validation**

**Scenario:** A new user is creating their first `input.yaml` file. They make a common mistake and misspell one of the main keys. They expect the program to fail fast with a helpful error, rather than running for minutes before crashing with an obscure internal error.
> **GIVEN** an `input.yaml` file where the user has made a typo in a required key, for instance:
> ```yaml
> system:
>   composiiton: "FePt"  # Note the typo: "composiiton" instead of "composition"
> ```
> **WHEN** the user attempts to execute the pipeline command, pointing to this invalid file.
>
> **THEN** the program must not start any of the expensive computational steps (like structure generation or DFT). It must exit immediately, within seconds.
> **AND** it must print a user-friendly error message to the console, not a frighteningly long Python stack trace.
> **AND** the error message must be specific and helpful. For example, it should clearly state that the field `composiiton` is not a valid or expected key within the `system` section.
> **AND** ideally, the error message should also suggest the correct key, for example: "Error in input file: Field 'composiiton' is not recognized. Did you mean 'composition'?" This provides a positive user experience and helps the user learn the correct configuration format.
