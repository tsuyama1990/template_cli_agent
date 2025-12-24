# CYCLE02 User Acceptance Testing (UAT)

## 1. Test Scenarios

This UAT plan focuses on verifying the new automation capabilities introduced in CYCLE02: the Config Expander heuristic engine and the automatic Structure Generator.

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C02-001 | End-to-End Workflow from Minimal Input (Alloy) | High     |
| UAT-C02-002 | Correct Material Classification and Strategy Selection | High     |
| UAT-C02-003 | Verification of Generated `exec_config_dump.yaml` | High     |
| UAT-C02-004 | Plausibility Check of Generated Structures | Medium   |
| UAT-C02-005 | Workflow with User-Provided Single Structure | Medium   |

---

**Scenario UAT-C02-001: End-to-End Workflow from Minimal Input (Alloy)**

*   **Description**: This is the primary "happy path" test for this cycle. It verifies that the system can run the full workflow, from automated configuration and structure generation through to model training, starting with a minimal user input for a binary alloy.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation from CYCLE01.
    *   A minimal `input.yaml` file containing only `system:\n  composition: "SiGe"`.
*   **Acceptance Criteria**:
    *   The command `uv run mlip-pipe run input.yaml` must execute and exit with a status code of 0.
    *   The system must create an `exec_config_dump.yaml` file.
    *   The log output must indicate that the system correctly classified the material as an "alloy".
    *   The log output must indicate that the SQS (Special Quasirandom Structures) generation method was used.
    *   The workflow must proceed to run DFT calculations on a set of generated SiGe structures.
    *   The workflow must successfully train and save an MLIP model file.
    *   No user intervention should be required after the initial command is executed.

---

**Scenario UAT-C02-02: Correct Material Classification and Strategy Selection**

*   **Description**: This test verifies that the Config Expander's heuristic engine correctly classifies different types of materials and selects the appropriate structure generation strategy for each.
*   **Preconditions**:
    *   A functional MLIP-AutoPipe installation.
    *   Four separate minimal `input.yaml` files for:
        1.  An alloy: `composition: "CuAu"`
        2.  A molecule: `composition: "H2O"`
        3.  An ionic solid: `composition: "NaCl"`
        4.  A covalent solid: `composition: "C"` (carbon)
*   **Acceptance Criteria**:
    *   When the workflow is run with the "CuAu" config, the log must show the material was classified as "alloy" and the generation strategy was "SQS".
    *   When the workflow is run with the "H2O" config, the log must show the material was classified as "molecule" and the generation strategy was "NMS" (Normal Mode Sampling).
    *   When the workflow is run with the "NaCl" config, the log must show the material was classified as "ionic" and the generation strategy was "AIRSS".
    *   When the workflow is run with the "C" config, the log must show the material was classified as "covalent" and the generation strategy was "Deep Rattling" or "Melt-Quench".

---

## 2. Behavior Definitions

**Behavior for UAT-C02-001: End-to-End Workflow from Minimal Input (Alloy)**

```gherkin
Feature: Fully Automated Workflow from Minimal Input
  As a user,
  I want to provide only the chemical composition of my material,
  So that the system automatically generates the data and trains a potential without any further input.

  Scenario: Generating a potential for a binary alloy
    GIVEN a minimal configuration file `input.yaml` with the composition "SiGe"
    AND a functional Quantum Espresso installation
    WHEN I execute the command `uv run mlip-pipe run input.yaml`
    THEN the process should complete successfully with an exit code of 0
    AND the system should log that it has classified the material as an "alloy"
    AND the system should log that it is using the "SQS" structure generator
    AND a complete `exec_config_dump.yaml` file should be created
    AND a final MLIP model file should be saved to the output directory.
```

---

**Behavior for UAT-C02-003: Verification of Generated `exec_config_dump.yaml`**

```gherkin
Feature: Heuristic Configuration Expansion
  As a developer and a user,
  I want the system to make intelligent, physically-based decisions when expanding the configuration,
  So that the resulting DFT and simulation parameters are sensible and robust.

  Scenario: Inspecting the auto-generated full configuration file
    GIVEN a minimal configuration file `input.yaml` with the composition "SiGe"
    WHEN I run the command `uv run mlip-pipe run input.yaml`
    THEN an `exec_config_dump.yaml` file must be generated
    AND this file must be a valid YAML file
    AND the `dft_compute` section must contain specific, non-null values for `ecutwfc`, `ecutrho`, and `pseudopotentials`
    AND the pseudopotentials listed must be appropriate for Si and Ge according to the SSSP protocol
    AND the `structure_generation` section must specify "sqs" as the chosen strategy.
```

---

**Behavior for UAT-C02-005: Workflow with User-Provided Single Structure**

```gherkin
Feature: Seeding Workflow with an Initial Structure
  As a user,
  I want to provide a known starting crystal structure for my material,
  So that the automatic generation process can create variations of that specific structure.

  Scenario: Starting the generation from a user-provided CIF file
    GIVEN a minimal configuration file `input.yaml` with `composition: "NaCl"`
    AND the config file also contains `initial_structure: "path/to/nacl.cif"`
    AND `nacl.cif` is a valid Crystallographic Information File for salt
    WHEN I execute the command `uv run mlip-pipe run input.yaml`
    THEN the system should use `nacl.cif` as the base for its structure generation (e.g., applying deformations or creating supercells)
    AND the workflow should complete successfully by generating a trained MLIP model.
```
