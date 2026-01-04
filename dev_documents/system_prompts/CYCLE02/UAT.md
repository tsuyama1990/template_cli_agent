# User Acceptance Test Plan: MLIP-AutoPipe Cycle 2

This document outlines the User Acceptance Testing (UAT) plan for Cycle 2 of the MLIP-AutoPipe project. The focus of this UAT is to verify the new advanced features, including the more sophisticated exploration and sampling algorithms, and to validate the functionality and usability of the new Web-based User Interface (Web UI).

## 1. Test Scenarios

This UAT will involve two distinct methods of interaction. Scenarios testing the new backend algorithms can still be performed via the Command-Line Interface (CLI) and will be documented in a Jupyter Notebook (`UAT_Cycle2_CLI.ipynb`). Scenarios testing the Web UI will be performed manually in a web browser, following a scripted set of user actions.

| Scenario ID | Scenario Name                           | Priority | Description                                                                                                                                                                                                     |
| :---------- | :-------------------------------------- | :------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| UAT-C2-01   | **Run Advanced Pipeline via CLI**       | High     | This scenario validates that the new backend features (Hybrid MD/MC and FPS Sampling) can be successfully configured and run from the command line, producing a valid final database.                              |
| UAT-C2-02   | **Interactive Run via Web UI**          | High     | This is the primary "happy path" for the new Web UI. The user will use their web browser to configure a pipeline, launch it, monitor its progress, and see the final result, all without touching the command line. |
| UAT-C2-03   | **Web UI Form Validation**              | High     | This scenario tests the responsiveness and error-handling of the Web UI's input form, ensuring a smooth user experience.                                                                                           |
| UAT-C2-04   | **Visualization of Results in Web UI**  | Medium   | This scenario verifies that after a successful run, the user can visually inspect the generated atomic structures directly within the web browser, providing immediate scientific feedback.                           |

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behaviour for each test scenario.

---

### **Scenario: UAT-C2-01 - Run Advanced Pipeline via CLI**

(Minimum 500 words)

While Cycle 2 introduces a Web UI, the command-line interface remains a critical tool for power users, scripting, and batch processing. This UAT scenario is designed to ensure that the new, sophisticated backend algorithms are not only functional but also fully configurable and accessible through the original CLI entry point. The user in this case is likely an experienced computational researcher who wants to leverage the advanced capabilities of the pipeline for a more complex materials system, such as a ternary alloy. They specifically want to enable the hybrid MD/MC exploration to enhance the sampling of different chemical orderings and use the Farthest Point Sampling (FPS) method to ensure the final dataset is as diverse and information-rich as possible. This test validates that the configuration schema extensions for `MCConfig` and `FPSConfig` are working correctly and that the `PipelineRunner` can correctly interpret these settings to activate the new logic in the `MDEngine` and select the `FPSSampler`. A successful run demonstrates that the core application remains a powerful, scriptable tool and that the new features have been seamlessly integrated into the existing architecture. The test will be considered a pass if a user can successfully complete an end-to-end run with these advanced features enabled and produce a valid ASE database.

**GIVEN** a directory containing valid configuration files for a ternary alloy system (e.g., Fe-Pt-Co).
**AND** the `exploration` configuration enables the hybrid MD/MC mode.
**AND** the `sampling` configuration specifies the `fps` method and provides a valid FPS configuration.

**WHEN** the user executes the `mlip-autopipec run` command from the terminal, pointing to the configuration directory.

**THEN** the command should execute without raising any errors and exit with a status code of 0.
**AND** the console log output should indicate that the Hybrid MD/MC engine is being used.
**AND** the console log output should indicate that Farthest Point Sampling is being performed.
**AND** a final ASE database file (`final_structures.db`) must be created and be non-empty.

---

### **Scenario: UAT-C2-02 - Interactive Run via Web UI**

(Minimum 500 words)

This scenario is the cornerstone of the user-facing improvements in Cycle 2. It tests the primary user journey for the new Web UI: configuring and running a complete data generation workflow interactively through a web browser. The target user for this scenario could be a student, an experimentalist, or anyone who is not a command-line expert but wants to leverage the power of MLIP-AutoPipe. The goal is to provide an experience that is intuitive, seamless, and visually engaging. The user will start by simply launching the web server with a single command and then perform all subsequent actions in their browser. They will navigate to the local web page, be presented with a clean and organized form, and fill in the parameters for their desired material system. Upon clicking "Run," they should receive immediate feedback that the job has started and be able to watch its progress. This test validates the entire Web UI stack, from the frontend JavaScript that captures the user's input and sends it to the backend, to the FastAPI server that receives the request, starts the pipeline in a background process, and provides status updates. A successful completion of this test demonstrates that the Web UI is not just a cosmetic addition but a fully functional and robust interface to the core power of the application. It validates that the complexity of the underlying pipeline is successfully abstracted away, lowering the barrier to entry and making the tool accessible to a wider audience.

**GIVEN** the MLIP-AutoPipe web server is running.
**AND** the user has opened the application's URL in their web browser (e.g., `http://127.0.0.1:8000`).

**WHEN** the user fills out the web form with valid parameters for a simple system (e.g., elemental Silicon).
**AND** the user clicks the "Run Pipeline" button.

**THEN** a "Job Started" notification or status message should immediately appear on the screen.
**AND** a progress indicator (e.g., a progress bar or text like "Status: Running Exploration...") should be displayed and update periodically.
**AND** after a few minutes, the status indicator should change to "Completed".
**AND** a "Success" message should be displayed, and a link or button to view the results should become available.
**AND** throughout the process, the web page should remain responsive and not freeze.

---

### **Scenario: UAT-C2-03 - Web UI Form Validation**

(Minimum 500 words)

A good user interface does more than just submit jobs; it actively helps the user avoid making mistakes. This UAT scenario focuses on the interactive validation and error-handling capabilities of the Web UI. It tests for a user experience that is smooth and forgiving. When a user is entering data into a web form, they should receive immediate feedback if their input is invalid, long before they attempt to submit the job. This "live validation" is a hallmark of a modern web application and is crucial for preventing user frustration. For example, if a user types a negative number into the "Temperature" field, the field's border should immediately turn red, and a helpful message like "Temperature must be positive" should appear next to it. This scenario will test several common input errors. This client-side validation, likely implemented in JavaScript, provides a fast and responsive feedback loop. Furthermore, the scenario also tests the backend validation. Even if a user could somehow bypass the client-side checks, the FastAPI backend should still validate the submitted data using the Pydantic models and return a structured, user-friendly error message. A successful outcome for this test means that it is virtually impossible for a user to submit a configuration that is syntactically invalid, and that when they make a mistake, the UI guides them gently toward correcting it.

**GIVEN** the user has the Web UI open in their browser.

**WHEN** the user types a non-numeric value (e.g., "abc") into the "Temperature" input field.
**THEN** the input field should be highlighted in red.
**AND** a message should appear next to the field stating "Please enter a valid number."
**AND** the "Run Pipeline" button may be disabled until the error is corrected.

**WHEN** the user clears the "Elements" input field, leaving it empty.
**THEN** the input field should be highlighted in red.
**AND** a message should appear stating "Please specify at least one element."

**WHEN** the user enters an invalid element symbol (e.g., "Xx") and attempts to click "Run Pipeline".
**THEN** the form submission should be prevented.
**AND** an error message should be displayed, either next to the field or in a general notification area, indicating that "Xx is not a valid element."

---

### **Scenario: UAT-C2-04 - Visualization of Results in Web UI**

(Minimum 500 words)

The ultimate goal of the pipeline is to produce atomic structures. Being able to see these structures is not just a cosmetic feature; it is a critical part of the scientific workflow. This UAT scenario verifies the user's ability to visually inspect the results of their calculation directly within the Web UI. After a pipeline run has successfully completed, the user needs to be able to assess the output. Did the MD simulation produce the kind of amorphous structure they expected? Are the generated surfaces correct? A simple "Completed" message is not enough. This test ensures that the UI provides a way to answer these questions visually. Upon completion of a job, a "View Results" button or link should appear. Clicking this should either take the user to a new results page or display the results dynamically on the current page. The core of this feature will be an embedded 3D molecular viewer, powered by a JavaScript library like `3Dmol.js` or `nglview`. The test will verify that one or more of the final, sampled structures are loaded and displayed in this viewer. The user should be able to interact with the 3D model, rotating, panning, and zooming to get a clear understanding of the atomic configuration. This feature closes the loop on the interactive user experience, taking the user from abstract configuration parameters to a tangible, visual representation of the final scientific product.

**GIVEN** a pipeline run has been successfully completed via the Web UI (as in UAT-C2-02).
**AND** the UI is displaying a "Completed" status.

**WHEN** the user clicks on the "View Results" or similar button.

**THEN** a new section of the page should become visible, or the user should be navigated to a results page.
**AND** a 3D viewer component should be displayed on the page.
**AND** within this component, an atomic structure from the final dataset should be rendered.
**AND** the user must be able to interact with the 3D structure using their mouse (e.g., click and drag to rotate, scroll to zoom).
