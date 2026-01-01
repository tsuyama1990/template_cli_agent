# CYCLE 02: UAT.md - Structure Generation and Configuration Expansion

## 1. Test Scenarios

This User Acceptance Test (UAT) for Cycle 02 is designed to validate the new user-centric features: the "Two-Tier Configuration" system and the automated `StructureGenerator`. The focus is on the user experience, demonstrating that the system can now be initiated with a vastly simplified input, removing the need for the user to understand complex DFT parameters or provide their own atomic structures. This UAT will showcase the system's new "intelligence" in interpreting minimal user requests and autonomously setting up a complete and valid workflow.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C02-01** | **Automated Workflow from Minimal Input for an Alloy**| **High** |

**Scenario UAT-C02-01: Automated Workflow from Minimal Input for an Alloy**

This scenario will test the complete, automated workflow for a representative binary alloy, Iron-Platinum (FePt), a material of significant scientific interest. The user experience begins with creating a very simple YAML file containing only the elements and a target composition. From this single, minimal input, the user will observe the system performing a series of automated actions:
1.  **Configuration Expansion:** The system will correctly identify the material as an 'alloy' and generate a comprehensive `exec_config_dump.yaml`, populating it with best-practice DFT parameters for Fe and Pt.
2.  **Structure Generation:** The system will automatically generate a set of diverse initial structures using the Special Quasirandom Structures (SQS) method, including strained variations.
3.  **Core Workflow Execution:** The system will then seamlessly feed these generated structures into the core labeling and training engines developed in Cycle 01.

A Jupyter Notebook, `C02_UAT.ipynb`, will be provided to make this process transparent and verifiable. The notebook will guide the user through creating the minimal `input.yaml`, running the single CLI command, and then inspecting all the generated artifacts—the expanded configuration, the new structures in the database, and the final trained model—confirming that the automation works end-to-end.

## 2. Behavior Definitions

**Feature:** Automated pipeline initiation from a minimal configuration file.

**Scenario:** Successful generation of a trained potential for an FePt alloy starting from a minimal user input.

*   **GIVEN** a clean working directory.
*   **GIVEN** a minimal YAML configuration file named `input.yaml` with the following content:
    ```yaml
    system:
      elements: ["Fe", "Pt"]
      composition: "FePt"

    simulation:
      temperature: [300, 1000]
    ```
*   **GIVEN** an empty ASE database file (`uat_c02.db`).
*   **GIVEN** a mock of the Quantum Espresso executable that produces valid, pre-computed output for FePt structures.

*   **WHEN** the user executes the command `uv run mlip-pipe --config-file input.yaml --database-file uat_c02.db`.

*   **THEN** the command should execute and exit with a success code (0).
*   **AND** a new file named `exec_config_dump.yaml` should be created in the working directory.
*   **AND** the `exec_config_dump.yaml` file should contain a `system` section with `structure_type: "alloy"`.
*   **AND** the `exec_config_dump.yaml` file should contain a `dft_compute` section with fully populated, physically reasonable parameters for FePt (e.g., `ecutwfc`, `pseudopotentials`).
*   **AND** the ASE database `uat_c02.db` should be populated with more than 5 new atomic structures.
*   **AND** all structures in the database should subsequently be updated to have the state 'labeled' and contain energy and force data.
*   **AND** a final trained model file named `model.ace` should be created.
*   **AND** the standard output should show log messages indicating the new initial stages: "Reading minimal config...", "Expanding configuration...", "Configuration expanded successfully.", "Starting Structure Generator...", "Generated XX initial structures.", before continuing with the familiar logs from Cycle 01.
