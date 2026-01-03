# User Acceptance Testing (UAT): Cycle 2

This document outlines the User Acceptance Testing (UAT) scenarios for the second cycle of the MLIP-AutoPipe project. The focus of this cycle is on the advanced scientific features, including the MD/MC Exploration Engine and Farthest Point Sampling, as well as the new Web User Interface. These tests are designed to be performed from a user's perspective, ensuring that the new features are not only functionally correct but also intuitive, robust, and deliver on their scientific promise.

## 1. Test Scenarios

These scenarios are designed to validate that the new, complex features are functioning correctly from a user's perspective and that the Web UI provides a seamless and intuitive experience for both novice and expert users.

### Scenario ID: UAT-C2-001
-   **Priority**: High
-   **Title**: End-to-End Pipeline Execution with MD Exploration and FPS Sampling via CLI
-   **Description**: This is the core scientific validation scenario for the entire project. It confirms that the main promise of the application—to automatically run a simulation and intelligently select a diverse dataset—is fulfilled. The user will configure a pipeline run for a simple, well-understood system like bulk Silicon. The configuration will specify using the new, advanced features: a short Molecular Dynamics simulation at an elevated temperature to ensure the atoms move significantly, followed by Farthest Point Sampling (FPS) to select the most diverse structures from the resulting trajectory. The success of this test demonstrates that the scientific heart of the application is working as expected. The user will verify this by observing tangible evidence of the pipeline's operation. The final structures in the output database should be visibly different from the initial, perfect crystal structure, and their calculated potential energies should show variation, proving they are distinct snapshots from a dynamic simulation. This confirms the system's capability to generate meaningful, evolved datasets, which is the foundational requirement for its use in training high-quality MLIPs. It provides confidence that the complex interplay between the MD engine, the MLIP calculator, and the FPS algorithm is correctly orchestrated.

### Scenario ID: UAT-C2-002
-   **Priority**: High
-   **Title**: Successful Pipeline Execution and Visualization via Web User Interface
-   **Description**: This scenario ensures that the new Web UI is a fully functional and user-friendly alternative to the Command Line Interface. The user will launch the web application and, using only their browser, replicate the setup from the previous scenario (UAT-C2-001). They will use the interactive form widgets to define the Silicon system, configure the MD simulation parameters, and select the FPS sampler. After launching the job, they will monitor its progress in real-time via the log viewer embedded in the UI. Upon successful completion, they will use the UI's built-in visualisation tool to load the resulting database and inspect one of the generated atomic structures in a 3D view. This scenario is critical because it validates the entire UI-to-backend workflow. It proves that the graphical interface is not merely a cosmetic feature but a complete control panel for the application. A successful test here means the application is now accessible to a much broader audience, including students, experimentalists, and scientists who may not be comfortable working exclusively on the command line. It lowers the barrier to entry and is a key step in making the tool a truly practical and widely adopted solution.

### Scenario ID: UAT-C2-003
-   **Priority**: Medium
-   **Title**: Verification of Automatic Ensemble Switching (Slab vs. Bulk)
-   **Description**: This scenario is designed to amaze a more expert user by demonstrating the system's built-in "intelligence" and scientific rigor. The user will perform two separate, comparable runs via the CLI. The first run will be for a standard bulk material, like the Silicon crystal from the previous scenarios. The second run will be for a surface slab, which the user can create by instructing the generator to add a vacuum layer along one of the cell axes. The key to this test is inspecting the detailed log files produced by both runs. The test is successful if the logs clearly show that the pipeline automatically detected the difference in geometry and, as a result, selected the correct thermodynamic ensemble for each simulation. The log for the bulk run should contain a message indicating that the NPT ensemble (constant pressure, constant temperature) was used, which is appropriate for allowing the cell volume to relax. The log for the slab run should, in contrast, show that the NVT ensemble (constant volume, constant temperature) was used. This is the physically correct choice for a slab system, as it prevents the vacuum layer from artificially collapsing during the simulation. This scenario proves that the `detect_vacuum` feature and the automatic switching logic are functioning correctly. It's a powerful demonstration that the tool isn't just a simple script but a sophisticated scientific instrument that helps users avoid common pitfalls and produce more accurate results.

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behavior for each UAT scenario in a structured, unambiguous format.

---

**Scenario**: `UAT-C2-001: End-to-End Pipeline Execution with MD Exploration and FPS Sampling via CLI`

*   **GIVEN** I am in a clean directory.
*   **AND** I have created a configuration file named `config_md.yaml` that specifies a short MD run using a fast classical potential (for testing purposes) and FPS sampling:
    ```yaml
    system:
      elements: [Si]
      composition: {Si: 1.0}
      lattice_constant: 5.43
      lattice: 'diamond'
    generation:
      num_structures: 1 # Start from a single perfect crystal
      supercell_size: [2, 2, 2]
    exploration:
      model_name: 'EMT' # Use a fast classical potential for this test
      temperature_k: 800
      num_steps: 500
      timestep_fs: 1.0
    sampling:
      method: 'fps'
      num_samples: 10
    database:
      path: Si_explored.db
    ```
*   **WHEN** I execute the command `mlip-autopipec run --config-path config_md.yaml` in my terminal.
*   **THEN** the process should complete successfully with an exit code of 0.
*   **AND** the terminal log should contain messages indicating that the `MDExplorer` and `FarthestPointSampler` components were executed.
*   **AND** a new file named `Si_explored.db` should be created in the current directory.
*   **AND** when I inspect the database with an ASE tool, it must contain exactly 10 structures, as requested in the sampling configuration.
*   **AND** when I calculate the potential energy of the first structure in the database, it must be significantly different from the potential energy of a pristine, unperturbed Silicon crystal, proving that the simulation has altered the structure.
*   **AND** the structures within the database should show a noticeable variance in their potential energies, indicating that FPS has selected diverse points from the simulation trajectory.

---

**Scenario**: `UAT-C2-002: Successful Pipeline Execution and Visualization via Web User Interface`

*   **GIVEN** the MLIP-AutoPipe package is installed and the `streamlit` dependency is available.
*   **WHEN** I launch the web application by running the command `mlip-autopipec web-ui` in my terminal.
*   **AND** I open the local URL (e.g., `http://localhost:8501`) provided by the command in my web browser.
*   **AND** I use the input forms, sliders, and dropdowns in the web page's sidebar to enter the same parameters as in the `config_md.yaml` file from scenario UAT-C2-001.
*   **AND** I click the "Run Pipeline" button in the UI.
*   **THEN** a log panel should appear on the main part of the page and begin displaying the real-time output of the backend pipeline execution.
*   **AND** after the pipeline execution finishes, the log panel should display a "Pipeline completed successfully" message.
*   **AND** a "Results Viewer" section should become visible on the page.
*   **AND** I should be able to select the `Si_explored.db` database from a dropdown or file selector in the Results Viewer.
*   **AND** upon selecting the database, an interactive 3D visualisation of one of the generated Silicon structures should appear, which I can rotate and zoom with my mouse.

---

**Scenario**: `UAT-C2-003: Verification of Automatic Ensemble Switching (Slab vs. Bulk)`

*   **GIVEN** I have two configuration files. The first, `config_bulk.yaml`, is identical to the file from `UAT-C2-001`.
*   **AND** the second configuration file, `config_slab.yaml`, is identical to `config_bulk.yaml` except for an added parameter in the `generation` section to create a vacuum layer (e.g., `vacuum_layer_angstroms: 10`).
*   **WHEN** I first execute the command `mlip-autopipec run --config-path config_bulk.yaml > bulk_run.log 2>&1` to capture all output.
*   **AND** I then execute the command `mlip-autopipec run --config-path config_slab.yaml > slab_run.log 2>&1` to capture all output.
*   **THEN** both commands should complete successfully with an exit code of 0.
*   **AND** when I search the contents of the `bulk_run.log` file, I must find a log message that explicitly states "Bulk system detected. Using NPT ensemble." or similar.
*   **AND** when I search the contents of the `slab_run.log` file, I must find a log message that explicitly states "Vacuum layer detected. Using NVT ensemble to preserve slab geometry." or similar.
*   **AND** the `bulk_run.log` must NOT contain the "NVT" message, and the `slab_run.log` must NOT contain the "NPT" message.
