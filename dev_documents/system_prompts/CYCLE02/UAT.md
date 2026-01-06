# User Acceptance Testing (UAT): MLIP-AutoPipe - Cycle 2

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for the second and final development cycle of the MLIP-AutoPipe project. This cycle introduces the advanced, science-focused features that form the core of the application's value proposition. The user persona remains the computational materials scientist, who will now rigorously test the system's ability to perform realistic, multi-stage simulations and intelligent data sampling. The UAT for this cycle is therefore designed to be more in-depth than that of Cycle 1, focusing on verifying the scientific correctness and robustness of the new Exploration and Sampling stages, as well as the successful integration of external Density Functional Theory (DFT) codes for the crucial Labeling stage. The goal is to give the user unshakable confidence that the tool is not just running, but is producing scientifically valid and high-quality results.

As in the previous cycle, a Jupyter Notebook (`UAT_CYCLE_02.ipynb`) will be the primary vehicle for this UAT. This interactive format is even more critical in this cycle, as it allows for the detailed, programmatic analysis and visualization required to verify the outputs of the complex new components. The notebook will guide the user through running a near-realistic workflow on a small, well-understood test system. It will then provide pre-written code snippets for interrogating the outputs in detail—for example, by plotting the energy evolution during a simulation or visualizing the structural differences highlighted by the Farthest Point Sampling (FPS) algorithm. This hands-on, guided analysis is essential for the user to accept the scientific validity of the tool's outputs.

| Scenario ID | Test Scenario | Priority |
| :--- | :--- | :--- |
| UAT-C2-01 | **Successful End-to-End Run with MD Exploration and FPS Sampling** | **High** |
| UAT-C2-02 | **Verification of External DFT Labeling Engine Integration** | **High** |
| UAT-C2-03 | **Verification of Interactive Web UI Functionality and Consistency** | **Medium** |

---

### **Scenario UAT-C2-01: Successful End-to-End Run with MD Exploration and FPS Sampling** (Min 300 words)

**Description:** This scenario is the flagship test for Cycle 2, designed to validate the core scientific functionality of the application's advanced workflow. The user will configure and execute a complete pipeline that leverages the newly implemented Molecular Dynamics (MD) exploration engine and the intelligent Farthest Point Sampling (FPS) stage. The primary goal is to provide the user with tangible proof that the system is not just executing a series of steps, but that these steps are producing a scientifically meaningful result. Success in this scenario means the user can verify that the MD simulation has adequately "shaken up" the initial structure, that the FPS algorithm has correctly identified and selected a structurally diverse subset of configurations from the resulting trajectory, and that the final database reflects this intelligent selection. This test is the most critical validation of the new features and proves that the tool can deliver on its promise of creating diverse and information-rich datasets.

**Jupyter Notebook Steps:**
1.  **System Setup:** The notebook will begin by defining a small, well-understood test system that is known to exhibit interesting dynamics, for example, a 4-atom Silicon Carbide (SiC) unit cell, which has distinct atomic elements and bonding.
2.  **Configuration Generation:** A Python cell will programmatically generate a `config_cycle2.yaml` file. This configuration will be carefully designed to test the new features. It will specify:
    *   **Generation:** Create just 1 initial, perfect SiC crystal structure.
    *   **Exploration:** Run a short but energetic MD simulation for 50 steps at a high temperature (e.g., 1500 K). This high temperature is crucial to ensure significant atomic movement and exploration of the potential energy surface, generating a trajectory with varied structures.
    *   **Sampling:** Use the `fps` method to select the 10 most diverse structures from the 50-frame trajectory.
    *   **Labeling:** Crucially, for this test, it will use the `mock` engine. This is a deliberate choice to isolate the validation to *only* the Exploration and Sampling logic, ensuring the test is fast and independent of any external DFT code installation.
3.  **Pipeline Execution:** The notebook will execute the pipeline using the `!mlip-autopipec run --config config_cycle2.yaml` command, displaying the live log output to the user, which should show the progress of the MD simulation steps.
4.  **Verification and Analysis:**
    *   **Successful Execution:** The user will first verify that the command completes with a zero exit code.
    *   **Database Content Analysis:** The notebook will then connect to the output database. The first and most important assertion is that the database contains *exactly 10 rows*. This confirms that the final dataset's size is determined by the `sampling` stage, not the `exploration` stage.
    *   **Diversity Verification (Qualitative and Quantitative):** This is the core of the UAT. To give the user tangible confidence in the FPS algorithm, the notebook will perform a comparative analysis. It will load the 10 structures selected by FPS from the database. It will also load the *full 50 structures* from the raw trajectory file (which the tool will be configured to preserve for this UAT). A simple Python function will then calculate a diversity metric for the FPS set vs. a randomly chosen set of 10 structures from the trajectory. This metric could be the standard deviation of the potential energies or the average pairwise structural distance. The user will be shown that the diversity score for the FPS set is significantly higher. To make this intuitive, the notebook will generate a plot of potential energy vs. simulation timestep for all 50 trajectory points, with the 10 FPS-selected points highlighted. The user should visually confirm that FPS has selected points from the peaks, troughs, and transition regions of the energy landscape, not just the stable low-energy regions.

---

### **Scenario UAT-C2-02: Verification of External DFT Labeling Engine Integration** (Min 300 words)

**Description:** This scenario is designed to test the `LabelingEngine`'s ability to correctly and robustly interface with a real, external DFT code. The success of this test is absolutely critical for the tool's practical utility in a real research environment. While the previous scenario verified the data generation workflow, this test verifies that the generated data can be accurately labeled with high-fidelity quantum mechanical calculations. This UAT will require the user to have a compatible DFT code (in this case, Quantum Espresso) installed and its executable (`pw.x`) available in their system's PATH. The test will run a minimal, focused workflow that generates a single, simple structure and then processes it with the real `dft` engine, allowing the user to scrutinize the entire process from input generation to output parsing.

**Jupyter Notebook Steps:**
1.  **Prerequisite Check:** The notebook will begin with a crucial setup cell. It will run the shell command `!which pw.x` to check if the Quantum Espresso executable is available. If the command fails, the notebook will print a clear, user-friendly message explaining that this test is being skipped and why, preventing user frustration.
2.  **Targeted Configuration:** A specific `config_dft.yaml` will be created. This configuration will be designed to be as simple and fast as possible to isolate the labeling functionality. It will:
    *   Generate only 1 simple, well-behaved structure (e.g., a 2-atom Silicon unit cell in its ground state).
    *   Explicitly disable and skip the Exploration and Sampling stages.
    *   **Labeling:** The `labeling` section will be configured to use `engine: dft`. It will also specify the necessary parameters for the calculation, such as the `dft_executable` path (if not in PATH) and settings for pseudopotentials and k-points.
3.  **Execution with Logging:** The notebook will run the CLI. The user will be instructed to observe the log output carefully, which should show the detailed log messages from the `DFTLabelingEngine` as it creates a temporary directory, writes the `qe_input.in` file, executes `pw.x`, and parses the results.
4.  **Multi-level Verification:**
    *   **Successful Run:** First, the user verifies the command completes successfully.
    *   **Database Verification:** The notebook will open the resulting database and retrieve the single entry. The user will be prompted to inspect the stored data. The key verification is that the `energy`, `forces`, and `stress` values are present, have the correct numerical types and shapes, and are physically plausible for the system (e.g., the energy for a 2-atom Si cell should be a specific, known negative value). This confirms that the output parsing logic is working correctly.
    *   **Calculation Artifact Inspection:** To provide maximum transparency, the tool will be configured to preserve the temporary working directory for this UAT run. The notebook will then contain a cell with `!ls -l temp_dft_run/` to show the user the files that were created. The user will be able to inspect the generated input file (`qe_input.in`) and the raw output file (`qe_output.out`) to manually confirm that the process was executed as expected and that the parsed values in the database match the values in the raw output.

---

### **Scenario UAT-C2-03: Verification of Interactive Web UI Functionality and Consistency** (Min 300 words)

**Description:** This scenario tests the usability, functionality, and, most importantly, the *consistency* of the alternative user interface: the Web UI. The goal is to ensure that a non-expert user can successfully configure and launch a pipeline run interactively through a browser, without needing to touch the command line or manually edit a YAML file. Furthermore, it aims to prove that a workflow launched from the Web UI produces bit-for-bit identical results to the same workflow launched from the CLI. This confirms that the UI is not a separate, divergent application but is instead a reliable and consistent front-end to the same robust, underlying core logic. Success in this scenario means the user can trust both interfaces to produce the same high-quality scientific results.

**Test Steps (Manual, Guided by Notebook):**
This UAT will be performed manually by the user, guided by detailed instructions provided in the Jupyter Notebook.

1.  **Launch the Web UI:** The notebook will instruct the user to open a new terminal window, navigate to the project directory, and run the new command `mlip-autopipec webui`. It will tell them the expected URL that will be printed to the terminal (e.g., `http://localhost:8501`) and instruct them to open this URL in their web browser.
2.  **Interactive Configuration:** The user will be presented with the web interface. The notebook will provide a clear set of instructions, asking them to replicate the exact configuration from the successful CLI test in Scenario UAT-C2-01. This will involve using the interactive widgets: selecting 'Cu' and 'Au' from a multi-select box, setting the composition, adjusting the MD temperature using a slider, and selecting 'fps' from a dropdown menu for the sampling method.
3.  **Execution and Monitoring:** The user will then be instructed to click the "Run Pipeline" button in the UI. The user should observe that the UI becomes active, providing real-time feedback. This feedback could be a progress bar that updates as the pipeline stages complete, or a text area where the application's log messages are streamed in real time.
4.  **Verification of Consistency:**
    *   Upon completion, the UI should display a clear "Success" message. The notebook will then instruct the user to return to the Jupyter environment.
    *   A special verification cell will be provided in the notebook. This cell will first connect to the database produced by the CLI run in UAT-C2-01 (`cycle_02_cli.db`). It will then connect to the database produced by the Web UI run (which will be named something like `cycle_02_webui.db`).
    *   The notebook will then programmatically compare the two databases. It will assert that they have the same number of rows. It will then iterate through the rows of both databases, comparing the energies, forces, and atomic positions. The test passes only if the contents of both databases are identical, proving that the Web UI is a consistent and reliable interface to the core application.

## 2. Behavior Definitions

**Feature: Advanced Scientific Workflow with Exploration and Sampling**

**Scenario: A user runs a full, complex pipeline using Molecular Dynamics and Farthest Point Sampling.**
*   **GIVEN** a valid YAML configuration file that specifies a 100-step Molecular Dynamics exploration at a temperature of 1000K.
*   **AND** the configuration explicitly requests that `20` structures be sampled from the resulting trajectory using the `'fps'` method.
*   **AND** the labeling engine is configured to use a real, external DFT code for high-fidelity calculations.
*   **WHEN** the user executes the `mlip-autopipec run` command with this configuration.
*   **THEN** the system must first execute the MD simulation for all 100 steps, generating a trajectory file containing 100 distinct atomic configurations.
*   **AND** after the simulation is complete, the system must apply the FPS algorithm to the full trajectory to select the 20 most structurally diverse configurations.
*   **AND** the system must then execute an external DFT calculation for **only** those 20 selected structures, not for all 100 structures in the trajectory.
*   **AND** the final output database file must contain exactly `20` entries.
*   **AND** each of these entries must contain valid, non-mock energy, force, and stress data as calculated by the external DFT code.

---

**Feature: Interactive and Consistent Web UI**

**Scenario: A user configures and runs a job via the Web UI and compares it to a CLI run.**
*   **GIVEN** the user has launched the Web UI via the `mlip-autopipec webui` command and has access to it in their browser.
*   **AND** the user also has a YAML file, `cli_config.yaml`, that defines a specific workflow.
*   **WHEN** the user first executes the workflow via the command line: `mlip-autopipec run --config cli_config.yaml`, resulting in an output file `cli_output.db`.
*   **AND** the user then uses the Web UI's interactive widgets to exactly replicate all the settings from `cli_config.yaml`.
*   **AND** the user clicks the "Run Pipeline" button in the Web UI, resulting in an output file `webui_output.db`.
*   **THEN** the application backend, triggered by the UI, must execute the exact same workflow logic as the CLI-triggered run.
*   **AND** the UI should provide visual feedback indicating that the job is running and has completed successfully.
*   **AND** upon completion, the contents of `webui_output.db` must be identical to the contents of `cli_output.db`.
