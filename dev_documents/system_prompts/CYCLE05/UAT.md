# Cycle 5 User Acceptance Testing (UAT): Web UI and User Experience

This UAT plan covers the final cycle of the project, which introduces a web-based graphical user interface (UI) for the MLIP-AutoPipe. The goal is to ensure that the UI provides a user-friendly, intuitive, and effective way to interact with the powerful backend pipeline developed in the previous cycles. The tests are designed from the perspective of a user who may not be comfortable with command-line tools.

## 1. Test Scenarios

The testing for this cycle will be hands-on and interactive. The primary "test artifact" will be the user's experience interacting with the live web application. The scenarios below will guide the user through the complete workflow, from building a configuration to analyzing the results, all within the browser.

| Scenario ID | Scenario Name                               | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C5-01   | Interactive Pipeline Configuration          | High     |
| UAT-C5-02   | Executing and Monitoring the Pipeline       | High     |
| UAT-C5-03   | Viewing and Analyzing Results               | High     |
| UAT-C5-04   | Responsive and Robust UI Behavior           | Medium   |

### UAT-C5-01: Interactive Pipeline Configuration

**Description:**
This scenario tests the core feature of the UI: the ability to create a valid pipeline configuration without editing a text file. The user will interact with the form widgets in the UI's sidebar to define a complete data generation workflow.

**Execution:**
1.  Launch the web application (e.g., by running `streamlit run src/mlip_autopipec/main_gui.py`).
2.  Navigate to the configuration panel in the sidebar.
3.  **Generator Settings:**
    *   Select "alloy" as the generator mode.
    *   Enter "Cu,Au" in the elements text box.
    *   Use the sliders and input boxes to set a lattice constant and supercell size.
4.  **Exploration Settings:**
    *   Enable the exploration stage.
    *   Use the slider to select a temperature of 600 K.
    *   Use the input box to set the number of steps to 500.
5.  **Sampling Settings:**
    *   Enable the sampling stage.
    *   Select "fps" from the mode dropdown.
    *   Use the slider to set the number of samples to 25.
6.  **Verification:** The user will visually confirm that all the widgets update correctly and that their selections are reflected on the page. There should be no errors or warnings displayed. The "Run Pipeline" button should be enabled, indicating that the UI considers the configuration to be complete and valid.

### UAT-C5-02: Executing and Monitoring the Pipeline

**Description:**
This scenario tests the UI's ability to launch the backend pipeline and provide real-time feedback to the user on its progress.

**Execution:**
1.  After completing the configuration in the previous scenario, click the "Run Pipeline" button.
2.  **Verification:** The user will immediately observe the following changes in the UI:
    *   The "Run Pipeline" button should become disabled to prevent multiple submissions.
    *   A status indicator on the page should change from "Idle" to "Running".
    *   A log panel should appear and begin to fill with the output from the backend process. The user should see messages indicating that the "Generation", "Exploration", and subsequent stages are starting.
    *   A progress bar should appear, starting near 0% and gradually increasing as the pipeline completes its stages.
    *   The UI should remain responsive; the user should be able to scroll the log panel while the backend is working.

### UAT-C5-03: Viewing and Analyzing Results

**Description:**
This scenario tests the final, and most rewarding, part of the UI: the results page. After the pipeline has finished, the UI should present a clear and interactive summary of the generated dataset.

**Execution:**
1.  Wait for the pipeline execution from the previous scenario to complete.
2.  **Verification:** Upon completion, the UI should update automatically:
    *   The status indicator should change to "Finished".
    *   A new "Results" section should appear.
    *   **Energy Plot:** The user will verify that a histogram of the potential energies of the final 25 structures is displayed. The plot should be clearly labeled.
    *   **Structure Viewer:** The user will verify that a gallery or list of 3D structure viewers is displayed. They should be able to interact with one of these viewers (rotate, zoom) to inspect one of the final, sampled atomic structures.
    *   The "Run Pipeline" button may become enabled again, allowing the user to start a new run.

### UAT-C5-04: Responsive and Robust UI Behavior

**Description:**
This scenario covers general UI robustness and user experience. It checks how the UI handles changes and potential issues.

**Execution:**
1.  Start a new session with the web application.
2.  **Configuration Changes:** Change a setting in the configuration panel (e.g., switch the generator from "alloy" to "ionic"). The user should verify that the relevant form fields dynamically update (e.g., the alloy-specific fields disappear and ionic-specific fields like "oxidation states" appear).
3.  **Stop/Cancel (if implemented):** If a "Cancel" button is available, the user will start a pipeline run and then click "Cancel". They will verify that the backend process is terminated and the UI returns to an "Idle" state.
4.  **Verification:** The user will confirm that the application feels smooth and responsive and that it provides clear guidance and feedback throughout the process.

## 2. Behavior Definitions

**Scenario: Interactive Pipeline Configuration and Execution**

*   **GIVEN** the user has opened the MLIP-AutoPipe web application.
*   **WHEN** the user fills out all the required fields in the configuration sidebar for all desired stages (Generation, Exploration, Sampling).
*   **AND** the user clicks the "Run Pipeline" button.
*   **THEN** the UI's status should immediately change to "Running" and the run button should be disabled.
*   **AND** a log panel should appear and display the real-time output from the backend process.
*   **AND** a progress bar should show the pipeline's progress.

**Scenario: Viewing and Analyzing Results**

*   **GIVEN** a pipeline has been successfully executed via the web UI.
*   **WHEN** the backend process completes.
*   **THEN** the UI status should automatically change to "Finished".
*   **AND** a "Results" section must become visible.
*   **AND** this section must contain a plot showing the energy distribution of the final dataset.
*   **AND** this section must contain at least one interactive 3D viewer displaying a structure from the final dataset.

**Scenario: Dynamic UI updates**

*   **GIVEN** the user is on the configuration page.
*   **WHEN** the user changes the "generator mode" dropdown from "alloy" to "ionic".
*   **THEN** the input fields relevant only to the alloy generator must disappear.
*   **AND** the input fields relevant to the ionic generator must appear.
