# User Acceptance Testing (UAT): Cycle 2

This document outlines the User Acceptance Testing (UAT) scenarios for the second and final development cycle of the MLIP-AutoPipe project. This cycle is pivotal as it introduces the most advanced scientific features—namely the hybrid MD/MC exploration engine and Farthest Point Sampling—and, crucially, the web-based user interface (UI). The UAT for this cycle is therefore twofold. First, it must verify that these new, sophisticated scientific features are not only functional from a software perspective but also provide a tangible, measurable benefit in the quality of the generated dataset. Second, it must ensure that the entire power of the pipeline is made accessible and intuitive through a well-designed graphical interface. These tests are designed from the perspective of the end-user and are focused on validating the overall experience, the scientific correctness of the output, and the seamless integration of the new UI with the powerful backend. Passing these UAT scenarios will signify that the project has successfully met its ambitious goal of creating a tool that is both powerful for experts and accessible for newcomers.

## 1. Test Scenarios

### Scenario ID: UAT-C2-001
*   **Priority:** High
*   **Title:** Interactive Pipeline Execution via Web UI for a Ternary Alloy System
*   **Description:** This scenario serves as the primary "happy path" and holistic validation for the web interface. It is the quintessential test of the GUI's effectiveness and its integration with the underlying pipeline. The scenario requires a user to configure and launch a moderately complex pipeline run—for a ternary alloy system—entirely through the graphical interface, without any reliance on the command line. This test will validate the complete end-to-end functionality of the GUI, starting from user input for a non-trivial chemical system, proceeding to the successful launch and monitoring of the backend process, and culminating in the verification of the final generated data. Its successful completion will demonstrate that the GUI is not merely a "demo" but a fully functional and robust entry point to the application's core capabilities, making the tool accessible to a much wider range of users, including those who are not experts in command-line environments.
*   **UAT Tool:** The primary interaction for this test will be with the Web UI itself, which is launched via the `mlip_autopipec-gui` command. To guide the user and to provide a rigorous, automated verification of the final output, a Jupyter Notebook (`C2_UAT_WebApp.ipynb`) will be provided. This notebook will serve a dual purpose: first, as a set of clear instructions, telling the user exactly what values to enter into the web form (e.g., elements `[Fe, Cr, Ni]`, specific compositions, and exploration parameters). Second, after the user has visually confirmed that the run has completed in the UI, the notebook will contain cells that programmatically access the resulting database file, and run a series of assertions to automatically and unambiguously verify its correctness (e.g., checking the number of structures, the elemental composition, and the presence of energy/force metadata).

### Scenario ID: UAT-C2-002
*   **Priority:** High
*   **Title:** Scientific Validation: Demonstrating the Superiority of Advanced Exploration with FPS Sampling
*   **Description:** This scenario is arguably the most important from a scientific standpoint. It is designed to provide a clear, quantitative demonstration that the key advanced features developed in Cycle 2—the hybrid MD/MC engine and Farthest Point Sampling—provide a tangible scientific benefit. The user will execute two parallel pipeline runs. The first will be a "baseline" run, using the simpler algorithms from Cycle 1 (standard high-temperature MD, random sampling). The second will be an "advanced" run, using the exact same number of steps and samples but enabling the new features. By directly comparing the structural diversity of the two resulting datasets, the user can empirically validate the project's central hypothesis: that intelligent exploration and sampling leads to higher quality data.
*   **UAT Tool:** A dedicated Jupyter Notebook (`C2_UAT_AdvancedFeatures.ipynb`) will be the exclusive tool for this highly analytical test. This notebook will be a self-contained experiment, performing the following steps:
    1.  It will programmatically define and write two distinct configuration files: `config_simple.yaml` and `config_advanced.yaml`.
    2.  It will then execute the `mlip_autopipec` command-line tool twice, once for each configuration, generating two separate databases: `simple_output.db` and `advanced_output.db`.
    3.  Next, it will load both of these databases and, for each one, it will compute high-dimensional SOAP (Smooth Overlap of Atomic Positions) descriptors for every structure.
    4.  The notebook will then calculate a quantitative metric of structural diversity for each dataset (e.g., the mean or variance of the pairwise distances in the SOAP feature space).
    5.  Finally, it will perform a critical assertion: that the diversity score for the `advanced_output.db` is statistically and significantly higher than that for the `simple_output.db`. To make the result intuitive, it will also generate plots (e.g., a 2D projection of the SOAP vectors using PCA) to visually demonstrate the broader and more uniform coverage of the configuration space achieved by the advanced method.

### Scenario ID: UAT-C2-003
*   **Priority:** Medium
*   **Title:** Verification of Real-time Progress Monitoring and Responsiveness in Web UI
*   **Description:** This scenario focuses on a critical aspect of user experience: feedback and responsiveness during a long-running task. When a user initiates a pipeline run, which could take minutes or even hours, it is essential that the UI provides clear and continuous feedback that the system is working. A static, unresponsive page can lead to user frustration and uncertainty, potentially causing them to kill the process prematurely. This test is designed to ensure that the monitoring features are implemented correctly and effectively, providing the user with confidence that their task is progressing as expected.
*   **UAT Tool:** For this test, the user will interact directly with the Web UI. They will be instructed by a markdown document to start a pipeline run with a configuration designed to take a noticeable amount of time (e.g., 20 initial structures and 100 MD steps). The user's task is to observe the UI's behavior during the execution of the pipeline. They will need to confirm several specific visual cues: that a progress bar appears and smoothly updates its value, that textual status messages are logged to the screen (e.g., "Now processing structure 5 of 20"), and that the "Run" button becomes disabled to prevent accidental multiple submissions. This validates the front-end's ability to receive and display real-time updates from the back-end process.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN for UAT-C2-001

**GIVEN** the user has successfully launched the web application by executing the `mlip_autopipec-gui` command in their terminal.
**AND** the user has navigated to the application's local URL in their web browser, presenting them with the main interface.
**WHEN** the user meticulously fills out the configuration form located in the web page's sidebar, entering parameters for a ternary alloy system as instructed.
**AND** the user clicks the main "Run Pipeline" button to initiate the process.
**THEN** the UI must provide immediate visual feedback, for instance, by displaying a progress bar or a status indicator message like "Pipeline starting...".
**AND** after the pipeline completes its execution, a clear and unambiguous "Success" message must be displayed prominently in the UI.
**AND** a new file named `final_dataset.db` must have been created in the application's designated working directory on the server's filesystem.
**AND** a programmatic inspection of this database file must confirm that it contains the correct number of structures, and that these structures are composed of the correct chemical elements, consistent with the parameters the user entered into the UI form.

### GIVEN/WHEN/THEN for UAT-C2-002

**GIVEN** a user has created two distinct and valid configuration files: `config_simple.yaml` (which specifies `"random"` for the sampling method) and `config_advanced.yaml` (which specifies `"fps"` for the sampling method and sets `use_mc_swaps: true`).
**WHEN** the user sequentially executes the command-line pipeline for both of these configuration files, resulting in the creation of two separate databases: `simple_output.db` and `advanced_output.db`.
**THEN** both of the pipeline runs must complete successfully, without any runtime errors.
**AND** a quantitative analysis of the structures within the two databases must show that the structural diversity of the configurations in `advanced_output.db` is measurably and significantly higher than the diversity of the configurations in `simple_output.db`.
**AND** as an expected scientific side-effect, the average potential energy of the structures in `advanced_output.db` should be higher than the average in `simple_output.db`, which serves as evidence that the hybrid MD/MC engine is successfully exploring higher-energy, less stable, and therefore more informative regions of the potential energy surface.

### GIVEN/WHEN/THEN for UAT-C2-003

**GIVEN** the user has launched the web application and has entered a valid configuration that will take at least one minute to run.
**WHEN** the user clicks the "Run Pipeline" button to start the process.
**THEN** the "Run Pipeline" button must immediately become disabled or hidden to prevent the user from accidentally starting multiple, concurrent runs.
**AND** a progress bar must appear on the UI, initialized to a value of 0%.
**AND** during the pipeline's execution, the user must be able to observe the progress bar's value smoothly and continuously increasing.
**AND** in addition to the progress bar, textual status updates that provide more context (e.g., "Step 1/4: Generating initial structures...", "Step 2/4: Running MD/MC exploration...") must be displayed and updated on the screen in a designated log area.
**AND** once the pipeline has finished its work, the progress bar must show a full 100%, and a final, clear "Process Completed Successfully" status message must be visible to the user.
