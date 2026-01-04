# CYCLE02/UAT.md

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 2 of the MLIP-AutoPipe project. The focus of this cycle is to validate the advanced, high-value features introduced: the hybrid MD/MC exploration engine, the Farthest Point Sampling (FPS) module, and the new web-based user interface. The UATs are designed to provide the user with a hands-on experience to confirm that these new capabilities are not only functional but also intuitive and effective. As in the previous cycle, the tests will be delivered via a well-documented Jupyter Notebook (`UAT_Cycle02.ipynb`), which will allow for interactive execution and visualization of the results.

| Scenario ID | Test Scenario | Priority |
| :--- | :--- | :--- |
| UAT-004 | Run a full pipeline (Explore, Sample, Store) via the CLI | High |
| UAT-005 | Interactively configure and run a pipeline via the Web UI | High |
| UAT-006 | Verify the superior structural diversity of Farthest Point Sampling | Medium |

### UAT-004: Run a full pipeline (Explore, Sample, Store) via the CLI

This test scenario validates the complete, end-to-end workflow, now incorporating the new exploration and sampling stages, executed from the command line. This is the primary UAT for power users who will rely on the CLI for scripting and automation. The Jupyter Notebook will provide a detailed walkthrough. First, it will guide the user in creating a comprehensive configuration file that includes the new `exploration` and `sampling` sections. This configuration will specify a short but meaningful MD simulation (e.g., heating a small crystal structure to a high temperature for a few dozen steps). The user will then execute the `run-pipeline` command from a notebook cell. The final and most important part of the scenario will be the analysis of the output. The notebook will contain code to load the final, sampled database and perform several checks. It will verify that the number of structures matches the requested `n_samples`. More qualitatively, it will plot the potential energy of the structures over the course of the simulation trajectory and highlight the points that were selected by the sampler, giving the user a visual confirmation that the process worked as expected.

### UAT-005: Interactively configure and run a pipeline via the Web UI

This scenario is the primary UAT for the new graphical user interface, a critical feature for making the tool accessible to a wider audience. The test is designed to assess the UI's intuitiveness, functionality, and integration with the backend pipeline. The user will be instructed to first launch the web server using the new `mlip-autopipec run-webui` command. They will then open their web browser and navigate to the provided local address. The UAT will then guide them through the process of using the web form to build a configuration for a simple run (e.g., generating a few initial structures and running a very short MD exploration). The user will click the "Start Run" button and will be asked to observe the UI's status feedback, verifying that it updates from "Idle" to "Running" and finally to "Completed". The test will be considered successful if the user can complete this entire workflow without needing to consult documentation and if the run produces the expected output database file on the filesystem. This provides a holistic validation of the entire web stack.

### UAT-006: Verify the superior structural diversity of Farthest Point Sampling

This scenario provides a focused and quantitative test of the Farthest Point Sampling algorithm, which is a key scientific feature of the project. The goal is to provide the user with tangible proof that FPS is more effective at creating a diverse dataset than simple random sampling. The Jupyter Notebook will orchestrate a compelling A/B test. First, it will run a single, moderately long exploration pipeline, saving the entire raw trajectory (e.g., 1000 structures) to a file. This trajectory will then serve as a common input for two different sampling runs. The first run will use the `random` sampling method to select 50 structures. The second run will use the `fps` method to select 50 structures from the *exact same trajectory*. The core of the UAT will be the analysis and visualization of the two resulting datasets. The notebook will calculate the SOAP descriptors for all structures in both datasets and then use a dimensionality reduction technique like Principal Component Analysis (PCA) to project the high-dimensional SOAP vectors into 2D space. It will then generate a scatter plot for both the random and FPS results. The user will be asked to visually compare the two plots. The expected and successful outcome is that the points on the FPS plot are visibly more spread out and cover the space more uniformly, while the points on the random plot are more clustered. This provides a clear, data-driven, and visually intuitive demonstration of the value added by the FPS algorithm.

## 2. Behavior Definitions

The following Gherkin-style definitions provide a formal, unambiguous description of the expected system behavior for the Cycle 2 test scenarios.

### Scenario: Successful End-to-End Run with Exploration and Sampling via CLI

**GIVEN** a valid YAML configuration file that specifies a `system` for Aluminum (Al).
**AND** the configuration specifies an `exploration` block with `run_md: true`, `temperature: 1200`, and `n_steps: 50`.
**AND** the configuration specifies a `sampling` block with `method: 'fps'` and `n_samples: 10`.
**AND** the configuration specifies a `db` path of `'./al_md_fps.db'`.

**WHEN** the user executes the `mlip-autopipec run-pipeline` command with this configuration.

**THEN** the application should successfully generate an initial structure, run a 50-step MD simulation at 1200K, sample the resulting trajectory, and save the final structures.
**AND** the application should exit with a status code of 0.
**AND** the final database file `'./al_md_fps.db'` should be created.
**AND** the database should contain exactly 10 structures.
**AND** the potential energies of the structures in the database should be different from each other, indicating that the MD simulation successfully explored the potential energy surface.

### Scenario: Launching and Running a Job from the Web UI

**GIVEN** the MLIP-AutoPipe web server has been started with the `mlip-autopipec run-webui` command.
**AND** the user has opened the web UI in their browser.

**WHEN** the user fills out the web form to configure a simple generation task for 5 structures of Silicon (Si).
**AND** the user clicks the "Start Run" button.

**THEN** the UI in the browser should immediately provide feedback, for example, by disabling the start button and showing a "Running..." status message.
**AND** a pipeline process should be initiated on the backend server.
**AND** after the backend process completes, the UI status should automatically update to "Completed".
**AND** a new database file, corresponding to the run, should be created on the server's filesystem and should contain exactly 5 structures of Silicon.

### Scenario: Verifying FPS vs. Random Sampling via Analysis

**GIVEN** a large trajectory of 1000 `ase.Atoms` objects has been generated and saved.
**AND** a first sampling run is performed on this trajectory using the `random` method to select 100 samples, resulting in `random_samples.db`.
**AND** a second sampling run is performed on the exact same trajectory using the `fps` method to select 100 samples, resulting in `fps_samples.db`.

**WHEN** the user runs the analysis script provided in the UAT Jupyter Notebook.

**THEN** the script should generate two 2D scatter plots from the PCA of the SOAP vectors for each dataset.
**AND** the plot corresponding to `fps_samples.db` should show points that are visibly more spread out and have fewer dense clusters compared to the plot for `random_samples.db`.
**AND** a quantitative metric, such as the average nearest-neighbor distance between points in the PCA space, should be higher for the FPS dataset than for the random dataset.
