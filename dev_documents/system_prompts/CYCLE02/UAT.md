# User Acceptance Testing (UAT): CYCLE02 - Configuration & Heuristics

## 1. Test Scenarios

This UAT focuses on verifying the new "Two-Tier Configuration" system. The "user" is now a regular scientist or engineer who wants to use the pipeline without needing to know the intricate details of DFT settings. The goal is to confirm that their simple input is intelligently expanded into a complete and correct configuration file that can drive the workflow.

| Scenario ID | Test Scenario                                                   | Priority |
| :---------- | :-------------------------------------------------------------- | :------- |
| UAT-C02-01  | Generate a full configuration from a minimal magnetic alloy input | High     |
| UAT-C02-02  | Generate a full configuration for a non-magnetic system         | High     |
| UAT-C02-03  | Verify user overrides are respected in the final configuration  | Medium   |

---

### **Scenario UAT-C02-01: Generate a full configuration from a minimal magnetic alloy input**

**(Min 300 words)**
This is the primary user story for this cycle. The user wants to generate a potential for a common magnetic alloy, FePt, but only knows the chemical composition. They will create the simplest possible `input.yaml` file and expect the system to make all the necessary expert decisions for them. The user will then meticulously inspect the auto-generated `exec_config_dump.yaml` to verify the correctness of the heuristic engine.

The user will start by creating an `input.yaml` containing only:
```yaml
system:
  elements: ["Fe", "Pt"]
  composition: "FePt"
```
After running `cdd run input.yaml`, they will open the generated `exec_config_dump.yaml`. The acceptance criteria are based on a checklist of expected values in this file. The user will first check the `system` section and confirm that `structure_type` has been correctly identified as `"alloy"`. Most importantly, in the `dft_compute` section, they will verify that `magnetism` is set to `"ferromagnetic"`, a critical setting for this material that a non-expert might forget. They will also verify that the `pseudopotentials` are set to a high-quality standard like `"SSSP_1.3_PBE_precision"`. Crucially, they must check that `ecutwfc` and `ecutrho` have been set to values appropriate for Fe and Pt (e.g., 90.0 and 720.0 Ry, respectively), demonstrating the SSSP heuristic is working. This successful test provides high confidence that the system can correctly configure a complex, magnetic material calculation from minimal user input.

---

### **Scenario UAT-C02-02: Generate a full configuration for a non-magnetic system**

**(Min 300 words)**
This scenario serves as a counterpoint to the first, ensuring the heuristic engine is not overly aggressive and can correctly handle simpler, non-magnetic materials. The user will test the system with a material like Silicon (Si), a canonical semiconductor. This tests the engine's ability to differentiate between material types and apply different sets of rules.

The user will create a minimal `input.yaml`:
```yaml
system:
  elements: ["Si"]
  composition: "Si"
```
Upon running the pipeline, they will again inspect the `exec_config_dump.yaml`. The acceptance criteria here are different from the FePt case. In the `system` section, they should see `structure_type` identified as `"covalent"`. In the `dft_compute` section, the critical check is that the `magnetism` parameter is either absent or explicitly set to `null` or `None`, confirming the system did not incorrectly enable magnetic calculations. Furthermore, the user will verify that the `ecutwfc` and `ecutrho` values are correct for Silicon based on the SSSP data (e.g., 40.0 and 320.0 Ry), which will be different from the values in the FePt test. They will also check that simulation parameters like `temperature_steps` are reasonable for Silicon, perhaps based on its high melting point. Successfully passing this test demonstrates the system's logic is not one-size-fits-all and can adapt its heuristics to different classes of materials.

---

### **Scenario UAT-C02-03: Verify user overrides are respected in the final configuration**

**(Min 300 words)**
This scenario tests the flexibility of the configuration system, ensuring that an expert user who wishes to provide more specific inputs can do so, and that their inputs are not ignored by the heuristic engine. The user will test this by providing a minimal input file but with one or two specific parameters overridden.

The user will create an `input.yaml` for FePt, but this time they will specify a custom DFT command and a non-default cutoff energy, simulating a situation where they have specific computational resources or wish to run a lower-precision test.
```yaml
system:
  elements: ["Fe", "Pt"]
  composition: "FePt"
dft_compute:
  command: "mpirun -np 8 pw.x"
  ecutwfc: 60.0
```
After running the pipeline, they will inspect the `exec_config_dump.yaml`. The acceptance criteria are twofold. First, they must verify that the overridden values are present. The `command` in the `dft_compute` section must be exactly `"mpirun -np 8 pw.x"`, and the `ecutwfc` must be `60.0`, not the SSSP-derived default of `90.0`. Second, they must verify that the *other* heuristics still worked correctly. For example, the `magnetism` should still be set to `"ferromagnetic"`, and the `ecutrho` should still be the SSSP default (or a value derived from the user's `ecutwfc`, e.g., `ecutrho: 480.0`). This test confirms that the `ConfigExpander` correctly layers user-provided values over its own defaults, providing a powerful balance between automation and expert control.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C02-01 - Generate a full configuration from a minimal magnetic alloy input**

*   **GIVEN** a clean project environment.
*   **AND** the user has created a file named `input.yaml` in the current directory.
*   **AND** `input.yaml` contains only the system's elemental composition: `system: {elements: [Fe, Pt], composition: FePt}`.
*   **WHEN** the user executes the command `cdd run input.yaml` from the command line.
*   **THEN** the system should create a new file named `exec_config_dump.yaml`.
*   **AND** this new file should be a valid YAML file.
*   **AND** when parsed, the `system.structure_type` key in the new file should have the value `"alloy"`.
*   **AND** the `dft_compute.magnetism` key should have the value `"ferromagnetic"`.
*   **AND** the `dft_compute.pseudopotentials` key should be set to a standard protocol like `"SSSP_1.3_PBE_precision"`.
*   **AND** the `dft_compute.ecutwfc` key should have a value of `90.0` (or the correct SSSP value for Pt).
*   **AND** the `dft_compute.command` key should be set to a reasonable default, like `"mpirun -np 4 pw.x"`.
*   **AND** the `simulation.temperature_steps` key should be populated with a default list of temperatures suitable for an alloy.
*   **AND** the pipeline should proceed to the (mocked) labelling and training steps using these generated parameters.

---

**Scenario: UAT-C02-02 & UAT-C02-03 - Combined Gherkin for a Non-Magnetic System with User Overrides**

*   **GIVEN** a clean project environment.
*   **AND** the user has created a file named `input.yaml` with the following content:
    ```yaml
    system:
      elements: ["Si"]
      composition: "Si"
    dft_compute:
      kpoints_density: 0.20
    ```
*   **AND** the user intends to test both the non-magnetic heuristic and the override mechanism.
*   **WHEN** the user executes the command `cdd run input.yaml` from the command line.
*   **THEN** the system should create a new file named `exec_config_dump.yaml`.
*   **AND** the `system.structure_type` key in the new file should have the value `"covalent"`.
*   **AND** the `dft_compute.magnetism` key should not be present, or its value should be `null`.
*   **AND** the `dft_compute.ecutwfc` key should have the correct SSSP value for Si (e.g., `40.0`).
*   **AND** the `dft_compute.kpoints_density` key **MUST** have the value `0.20`, reflecting the user's override, not a default value.
*   **AND** all other required keys, such as `dft_compute.command` and `mlip_training.model_type`, should be populated with their standard default values.
*   **AND** the rest of the pipeline should continue execution using this merged configuration of heuristics and user overrides.
