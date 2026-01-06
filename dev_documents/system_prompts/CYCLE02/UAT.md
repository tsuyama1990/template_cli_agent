# User Acceptance Test (UAT) Plan: CYCLE 02

This document outlines the User Acceptance Testing scenarios for the advanced features implemented in Cycle 2 of the MLIP-AutoPipe project. The focus of this UAT is to validate, from a user's perspective, the functionality, scientific correctness, and usability of the newly introduced `Exploration` and `Sampling` stages, as well as the brand-new Web UI.

As in Cycle 1, a Jupyter Notebook is the recommended environment for performing the backend-focused parts of this UAT. This will allow the user to run a full, complex pipeline, programmatically inspect the intermediate trajectory files, and perform quantitative analysis on the final sampled structures. This provides a comprehensive and interactive validation experience that goes beyond simple pass/fail checks. For the UI part, manual interaction with the web application will be necessary.

## 1. Test Scenarios

| ID    | Priority | Test Scenario                                   |
| :---- | :------- | :---------------------------------------------- |
| UAT-4 | High     | Verify successful end-to-end run of the full (4-stage) pipeline via CLI. |
| UAT-5 | High     | Verify the superior structural diversity of the FPS-sampled dataset compared to random sampling. |
| UAT-6 | Medium   | Verify the core functionality and usability of the Web UI for running a full pipeline. |

### Scenario Details

#### **UAT-4: Successful end-to-end run of the full (4-stage) pipeline via CLI.** (High Priority)
This is the ultimate "happy path" test for the entire application's backend. Its purpose is to verify that all four stages (Generation, Exploration, Sampling, Storage) are now correctly integrated and can execute in sequence without errors for a realistic use case. The user will configure a complete workflow, including a short but meaningful molecular dynamics simulation and the advanced Farthest Point Sampling algorithm. The successful completion of this test is a critical milestone. It demonstrates that the complex computational backend is robust, that the data flows correctly between all stages (from memory, to disk, and back), and that the parallel processing of the exploration stage is functioning as expected. It provides the user with the confidence that they can reliably use the command-line tool for their serious research tasks, automating what was previously a complex, error-prone, and manual process. This test validates the core engine of the entire product.

#### **UAT-5: Verify the superior structural diversity of the FPS-sampled dataset.** (High Priority)
This scenario is arguably the most important for validating the *scientific value* of the software. The key promise of the `Exploration` and `Sampling` engines is not just to produce data, but to produce a *diverse and informative* dataset. This UAT provides a tangible demonstration of this intelligence. It will involve a direct comparison of the outputs from two full pipeline runs that are identical in every way except for the sampling method: one will use the basic Random Sampling, and the other will use the advanced Farthest Point Sampling (FPS). By visualizing the structures and, more importantly, by performing a quantitative analysis of the structural similarity within each dataset, the user should be able to prove that the FPS-selected dataset is superior. They should find that it contains a richer variety of configurations—for instance, a healthy mix of low-energy stable states, high-energy distorted states, and perhaps even some transition states. This test is designed to amaze the user by providing clear, undeniable evidence of the "intelligence" of the sampling algorithm and the direct scientific benefit it provides for training better MLIPs.

#### **UAT-6: Verify the core functionality and usability of the Web UI.** (Medium Priority)
This scenario focuses on the accessibility and user experience of the graphical interface. The user will step away from the command line and interact with the Web UI to build a full pipeline configuration from scratch, launch the pipeline, and observe the results. The goal is to ensure that the UI is intuitive, responsive, and provides helpful feedback throughout the process. The user will test the interactive input widgets (e.g., sliders for temperature, text boxes for elements), the "Run Pipeline" functionality, and the real-time status updates that should log the progress of the backend. A successful test means the user can perform a complete, complex data generation task without needing to touch the command line or manually edit a single line of a text file. This validation is critical for confirming that the power of the backend has been made accessible to a broader, non-expert audience. The structure visualization feature will be a key point of amazement, allowing for instant visual confirmation and exploration of the data that the user has just created.

## 2. Behavior Definitions

---

**Scenario: UAT-4 - Successful end-to-end run of the full (4-stage) pipeline via CLI**

*   **GIVEN** I have created a valid YAML configuration file named `config_full.yaml`.
*   **AND** this file specifies a complete, four-stage pipeline:
    ```yaml
    system:
      generator_type: "alloy"
      elements: ["Ni", "Ti"]
      composition: {"Ni": 0.5, "Ti": 0.5}
      num_structures: 4
    exploration:
      md_type: "nvt"
      temperature_k: 1000.0
      timestep_fs: 1.0
      num_steps: 20
      mlip_model: "emt" # Use the fast, built-in ASE potential for testing
    sampling:
      sampler_type: "fps"
      num_samples: 5
    database_path: "uat_full_pipeline.db"
    ```
*   **WHEN** I execute the command `mlip-autopipec run-pipeline --config-path config_full.yaml` from my terminal.
*   **THEN** the command should complete successfully with an exit code of 0.
*   **AND** during execution, the terminal should display informative progress updates, clearly indicating when it is running the `Exploration` and `Sampling` stages.
*   **AND** upon completion, a database file named `uat_full_pipeline.db` should be created.
*   **AND** when I inspect this database, it should contain two distinct groups of data, identifiable by a key or tag: the `initial` structures and the `final_dataset`.
*   **AND** the `final_dataset` group in the database should contain exactly `5` structures, as requested in the sampling configuration.

---

**Scenario: UAT-5 - Verify the superior structural diversity of the FPS-sampled dataset**

*   **GIVEN** I have successfully run the full pipeline twice, once with `sampler_type: "random"` and once with `sampler_type: "fps"`. All other exploration settings are identical.
*   **AND** as a result, I have two database files: `uat_random.db` and `uat_fps.db`, each containing a final dataset of 20 structures.
*   **WHEN** I execute a Jupyter Notebook that loads the final datasets from both databases.
*   **AND** the notebook calculates a structural fingerprint (e.g., SOAP vectors) for every structure in both datasets.
*   **AND** the notebook then calculates a similarity metric for each dataset (e.g., the average cosine similarity of the fingerprint vectors within the set).
*   **THEN** the structures loaded from `uat_fps.db`, when visualized, should be visibly more diverse—showing a wider range of atomic orderings, lattice distortions, and coordination environments—than the more uniform structures from `uat_random.db`.
*   **AND** the calculated average similarity for the `uat_fps.db` dataset should be quantitatively lower than the average similarity for the `uat_random.db` dataset, providing mathematical proof that the FPS algorithm successfully selected a more diverse set of structures.

---

**Scenario: UAT-6 - Verify the core functionality and usability of the Web UI**

*   **GIVEN** the MLIP-AutoPipe Web UI application has been started and is now accessible in my web browser.
*   **WHEN** I use the interactive UI widgets in the "System Configuration" section to enter the elements `["Cu", "Au"]` and specify `10` initial structures.
*   **AND** I use the UI widgets in the "Exploration" section to set the temperature to 800K and the number of MD steps to 50.
*   **AND** I use the UI widgets in the "Sampling" section to select the "FPS" sampler and request `15` final samples.
*   **AND** I click the main "Run Pipeline" button.
*   **THEN** a progress indicator or spinner should appear immediately, and status messages should begin to be displayed in a log area, showing the pipeline's progress through the "Generation", "Exploration", "Sampling", and "Storage" stages.
*   **AND** after the pipeline finishes, a "Success!" message should be clearly displayed.
*   **AND** a 3D structure viewer component in the UI should load and display one or more of the final, sampled structures.
*   **AND** I should be able to interactively rotate and zoom in on the visualized atomic structure using my mouse.
*   **AND** a "Download Results" button should be available, and clicking it should allow me to save the generated database file to my local machine.
