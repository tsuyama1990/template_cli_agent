# CYCLE 05: UAT.md - Advanced Features and Usability

## 1. Test Scenarios

This final User Acceptance Test for the initial project scope is designed to validate the advanced scientific capabilities and usability improvements introduced in Cycle 05. The UAT will demonstrate the system's expanded applicability to a wider and more complex range of materials, specifically focusing on the newly implemented features for magnetic systems. It will also showcase the enhanced user experience provided by the improved command-line interface.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C05-01**| **Automated Workflow for a Magnetic Material (Iron)** | **High** |

**Scenario UAT-C05-01: Automated Workflow for a Magnetic Material (Iron)**

This scenario will focus on a classic and important magnetic material: body-centered cubic (BCC) Iron (Fe). The user will provide a minimal `input.yaml` containing only `elements: ["Fe"]`. This will trigger the system's advanced heuristics. The user will observe that the system not only identifies the material and sets up the standard workflow but also automatically recognizes Iron as a magnetic element and enables the special DFT calculation parameters required to model it accurately.

The user experience is centered on the system's "expert knowledge." The user does not need to know the specific Quantum Espresso parameters for a spin-polarized calculation; the system handles it for them. The `C05_UAT.ipynb` Jupyter Notebook will guide the user through this process. It will involve:
1.  Creating the simple `input.yaml` for Iron.
2.  Running the CLI command.
3.  Inspecting the `exec_config_dump.yaml` to verify that the system automatically set `magnetism: "ferromagnetic"`.
4.  Observing the new, user-friendly progress bars in the CLI output during the labeling and training stages.
This UAT provides a powerful demonstration of the system's core philosophy: encoding expert knowledge to free the user from needing to manage complex, domain-specific details, thereby making advanced materials simulation more accessible.

## 2. Behavior Definitions

**Feature:** Automated handling of magnetic materials and improved user feedback.

**Scenario:** Successful generation of a trained potential for BCC Iron, initiated from a minimal config, with clear user feedback.

*   **GIVEN** a clean working directory.
*   **GIVEN** a minimal YAML configuration file named `input.yaml` with the content:
    ```yaml
    system:
      elements: ["Fe"]
    simulation:
      temperature: [300, 1000]
    ```
*   **GIVEN** an empty ASE database file (`uat_c05.db`).
*   **GIVEN** a mock of the Quantum Espresso executable that is programmed to inspect the input file it receives.

*   **WHEN** the user executes the command `uv run mlip-pipe --config-file input.yaml --database-file uat_c05.db`.

*   **THEN** the command should execute and exit with a success code (0).
*   **AND** a file `exec_config_dump.yaml` should be created containing `magnetism: "ferromagnetic"`.
*   **AND** the standard output of the command should display a progress bar (e.g., from the `rich` library) during the "Labeling Structures" phase.
*   **AND** the mock for the Quantum Espresso executable must confirm that the input file it was passed contains the line `nspin = 2` and a non-zero `starting_magnetization` value, which are required for a magnetic calculation.
*   **AND** a final trained model file named `model.ace` should be created.
