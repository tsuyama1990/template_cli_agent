# User Acceptance Testing (UAT): MLIP-AutoPipe Cycle 2

This document outlines the User Acceptance Testing (UAT) scenarios for the second cycle of the MLIP-AutoPipe project. The primary focus of this UAT is to meticulously verify the newly implemented advanced features. This includes the sophisticated hybrid MD/MC exploration engine, the intelligent Farthest Point Sampling (FPS) algorithm, the critical support for Machine Learning Interatomic Potential (MLIP) models, and the highly anticipated new Web User Interface. The tests are designed not only to ensure that these features are functionally correct but also that they provide a tangible and significant benefit to the end-user and contribute to a smooth, intuitive, and powerful user experience.

## 1. Test Scenarios

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C2-01   | Successful Pipeline Run with Advanced Features | High     |
| UAT-C2-02   | Interactive Configuration and Visualization via Web UI | High     |
| UAT-C2-03   | Quantitative Verification of FPS Diversity  | Medium   |

### UAT-C2-01: Successful Pipeline Run with Advanced Features

**(Min 300 words)**
This scenario is designed to validate the successful integration and orchestration of all new scientific features in a single, comprehensive run executed via the command-line interface. The user persona for this test is an advanced materials scientist whose research focuses on solid-state ion conductors. Their immediate goal is to generate a high-quality training dataset for a complex ionic material, ceria (CeO2), using a pre-trained MACE model they have developed. The scientific objective is to explore potential oxygen vacancy formation and migration pathways, so the user will enable the new hybrid MD/MC feature with a defined probability for Ce-O atom swaps to simulate defect dynamics. To ensure the resulting dataset is as information-rich and non-redundant as possible, the user will also enable Farthest Point Sampling. The user will begin by authoring a new YAML configuration file that specifies the `ionic` generator type, the file path to their trained MACE model, enables ZBL mixing as a safety measure for high-temperature stability, and sets specific parameters for both the Monte Carlo moves and the FPS sampler's SOAP descriptors. The core expectation is that the `mlip-autopipec` CLI can correctly parse this complex, multi-faceted configuration and execute the entire pipeline from start to finish without errors. The user will be closely monitoring the console output, expecting to see explicit log messages that confirm the MACE model is being loaded, that the active calculator is a mixed potential, and that MC swap moves are being attempted and reported with an acceptance ratio during the exploration phase. Upon successful completion, the resulting output database should contain a set of atomic structures that are not only physically valid but also structurally diverse, reflecting the more complex and powerful exploration dynamics of the hybrid engine. This test is the cornerstone of Cycle 2's validation, as it provides definitive proof that all the sophisticated backend components are working together harmoniously to deliver a scientifically superior result that was impossible to achieve with the simpler pipeline from Cycle 1.

### UAT-C2-02: Interactive Configuration and Visualization via Web UI

**(Min 300 words)**
This scenario focuses on the usability, functionality, and overall user experience of the new Web User Interface. The user persona is a graduate student or a researcher from a more experimental background who is less comfortable with command-line tools and YAML files, but wants to leverage the power of MLIP-AutoPipe for their work. The test begins with the user launching the Web UI by running a simple, memorable command like `mlip-autopipec ui`. Their default web browser should open automatically to a clean, well-organized, and intuitive interface. The user will then interact with a variety of standard web form elements—dropdown menus to select the generator type, sliders to set temperature and simulation steps, and text fields to enter chemical elements—to build a configuration for a simple alloy system, similar to the one in Cycle 1. The UI must provide helpful tooltips for less obvious parameters and should provide clear, real-time validation feedback (e.g., highlighting a field in red and showing an error message if they enter a negative temperature). Once the configuration is fully specified, the user will click a prominent "Start Pipeline Run" button. The UI must then transition to a "running" state, showing a progress indicator and a streaming log of messages from the backend process, so the user is kept informed and knows the task is running. After the pipeline finishes, the UI should automatically detect this and transition to a "results" view. The user expects to see a summary of the run and a paginated table of the generated structures. Upon clicking a row in this table, a 3D interactive visualization of that specific atomic structure must appear in an adjacent panel, allowing the user to rotate, zoom, and inspect the atoms. This scenario is crucial for validating that the tool is accessible to a broader scientific audience and that the UI provides a complete, seamless, and satisfying "round-trip" experience from interactive configuration to graphical result visualization.

### UAT-C2-03: Quantitative Verification of FPS Diversity

**(Min 300 words)**
This scenario is designed to provide a quantitative, objective, and scientifically rigorous validation of the benefits of the Farthest Point Sampling (FPS) algorithm. The user persona is a methods developer or a particularly thorough scientist who wants to prove that FPS creates demonstrably "better" datasets than the standard random sampling approach. To achieve this, the user will perform a carefully controlled computational experiment. First, they will define a moderately complex system (e.g., a Mo-Ta-W high-entropy alloy) and run the full pipeline with the `RandomSampler` enabled, generating a baseline database named `random_dataset.db`. Second, they will run the exact same pipeline—using the same initial structures and identical exploration parameters—but this time with the `FPSSampler` enabled, generating a second database named `fps_dataset.db`. The user will then use the provided UAT Jupyter Notebook for Cycle 2. This notebook will contain pre-written Python functions to load both databases, calculate the SOAP descriptors for every structure in each database using the same SOAP parameters, and then compute an average diversity score for each dataset. This score will be defined as the average pairwise Euclidean distance between all structures in the SOAP feature space. The core expectation is that the calculated diversity score for `fps_dataset.db` will be measurably and significantly higher than the score for `random_dataset.db`. A successful result would provide concrete, quantitative evidence that the FPS implementation is working as intended and is successfully fulfilling its primary mission: to select a more diverse and non-redundant set of structures from the same underlying trajectory data. This is a key selling point for the advanced feature set, as it directly translates to more efficient MLIP training and potentially more robust models. The notebook may also include a visualization (e.g., a 2D projection using PCA or t-SNE) to qualitatively show that the FPS points are more spread out in the feature space.

## 2. Behavior Definitions

**(Min 500 words)**

### Gherkin-style Definitions

**Scenario: UAT-C2-01 - Successful Pipeline Run with Advanced Features**

*   **GIVEN** a user has a valid YAML configuration file named `advanced_config.yml`.
*   **AND** the file specifies `generator_type: 'ionic'` for CeO2.
*   **AND** the file provides a valid path to a MACE model file via `mlip_model_path`.
*   **AND** the file has `use_zbl_mixing: true`.
*   **AND** the file specifies a non-zero `swap_probability` under a `mc_moves` section.
*   **AND** the file selects the `sampler_type: 'fps'` and provides a valid `fps_config`.
*   **AND** the user is in a terminal session.
*   **WHEN** the user executes the command `mlip-autopipec run --config-file advanced_config.yml`.
*   **THEN** the application should print log messages to the console indicating that the MACE model is being successfully loaded.
*   **AND** the log messages should confirm that the active ASE calculator is a mixed potential (MLIP+ZBL).
*   **AND** the log output for the exploration stage should periodically show messages indicating that MC swap moves are being attempted, along with an acceptance ratio.
*   **AND** after the exploration stage, the application should print messages indicating the start and completion of "Calculating SOAP vectors" and "Running Farthest Point Sampling".
*   **AND** the application must exit with a success code (0).
*   **AND** the specified output database file must be created and contain the expected number of structures.

**Scenario: UAT-C2-02 - Interactive Configuration and Visualization via Web UI**

*   **GIVEN** a user runs the command `mlip-autopipec ui` in their terminal.
*   **AND** a web browser automatically opens, displaying the MLIP-AutoPipe user interface with a sidebar for configuration.
*   **WHEN** the user selects "alloy" from the `generator_type` dropdown menu in the sidebar.
*   **AND** the user enters `Cu, Au` into the "Elements" text field.
*   **AND** the user adjusts a slider for "Temperature" to `500` K.
*   **AND** the user clicks the "Run Pipeline" button.
*   **THEN** the main area of the UI should display a progress indicator and a text area showing real-time log messages from the backend pipeline.
*   **AND** upon completion, a success message "Pipeline finished successfully" must be displayed.
*   **AND** a data table should appear in the UI, populated with rows corresponding to the generated structures.
*   **AND** when the user clicks on any row in the table, a 3D visualization widget must appear and render the atomic structure for that row, allowing for rotation and zoom.

**Scenario: UAT-C2-03 - Quantitative Verification of FPS Diversity**

*   **GIVEN** a user has two ASE database files, `random.db` and `fps.db`.
*   **AND** these databases were generated from identical pipeline runs where the only difference was the `sampler_type` configuration.
*   **AND** the user has opened the UAT Jupyter Notebook for Cycle 2.
*   **WHEN** the user executes the notebook cell responsible for loading the databases and calculating their respective diversity scores.
*   **THEN** the notebook should successfully connect to both database files and load all structures into memory without errors.
*   **AND** the notebook must output two clearly labeled floating-point numbers: "Random Sampler Diversity Score" and "FPS Diversity Score".
*   **AND** the value of the "FPS Diversity Score" must be numerically greater than the value of the "Random Sampler Diversity Score".
*   **AND** to be considered a significant success, the FPS score should be at least 15% greater than the random score.
*   **AND** the notebook may optionally display a 2D plot showing two distinct clusters of points, with the FPS points visibly more spread out than the random points.
