# CYCLE02 USER ACCEPTANCE TEST (UAT): Advanced Exploration and User Interface

This document outlines the User Acceptance Testing (UAT) plan for the second major release of the MLIP-AutoPipe tool. The focus of this UAT is to validate the new, advanced scientific features and the new web-based user interface. These tests will ensure that the tool is not only more powerful and versatile but also more accessible to a wider range of users.

## 1. Test Scenarios

The following scenarios will be executed by the user to verify the new capabilities. As with the previous cycle, a Jupyter Notebook (`UAT_CYCLE02.ipynb`) will be provided to facilitate the execution and verification of the command-line scenarios. The final scenario will be a manual test of the web interface.

---

**Scenario ID:** UAT-C02-001
**Priority:** High
**Title:** Validate Hybrid MD/MC Engine for Enhanced Alloy Sampling

**Description (Min 300 words):**
This scenario is designed to test the core scientific advancement in Cycle 2: the hybrid Molecular Dynamics/Monte Carlo (MD/MC) exploration engine. Standard MD can sometimes get trapped in local energy minima, especially in complex alloys. The addition of MC atom swaps is intended to overcome these energy barriers and allow for a much broader exploration of the material's phase space. The user will configure a pipeline run for a binary alloy, such as Nickel-Aluminum (NiAl), which is known to have multiple ordered phases. They will enable the hybrid engine with a reasonable probability for atom swap moves.

The test's success will be evaluated by analyzing the output trajectory. The user will write a short script within the provided Jupyter Notebook to track the local environment of a few selected atoms over the course of the simulation. In a standard MD run, atoms would largely stay close to their initial neighbors. In a successful hybrid MD/MC run, the script should be able to detect that atoms have swapped positions, indicating that the MC moves are working correctly. Furthermore, the user will compare the final sampled structures with those from a standard MD run (from Cycle 1). The expectation is that the hybrid run will produce a more diverse set of configurations, potentially including different local ordering motifs, which are critical for training a robust MLIP. This test confirms that the most important new scientific feature is not only functional but delivers its intended benefit.

---

**Scenario ID:** UAT-C02-002
**Priority:** High
**Title:** End-to-End Pipeline Run for an Ionic Crystal (NaCl)

**Description (Min 300 words):**
This scenario validates the tool's extended capability to handle different classes of materials, specifically ionic crystals. The user will configure the pipeline to generate a dataset for Sodium Chloride (NaCl). This requires using the new `IonicGenerator`. The configuration will specify the cation (`Na`) and anion (`Cl`) and other relevant crystal structure information. The user will then run the full end-to-end pipeline.

The primary verification will be to inspect the generated structures to ensure they adhere to the physical constraints of ionic materials. First, the user will check that all generated structures, both initial and final, are perfectly charge-neutral, meaning they contain an equal number of Sodium and Chlorine atoms. This is a critical invariant that the `IonicGenerator` must enforce. Second, the user will analyze the radial distribution function, `g(r)`, for the final structures. They should observe that Na ions are preferentially surrounded by Cl ions (and vice-versa), and that ions of the same charge (Na-Na, Cl-Cl) repel each other and are rarely found as nearest neighbors. This confirms that the simulation correctly captures the fundamental physics of ionic bonding. This scenario's success demonstrates that the framework's modular design works and that it can be easily extended to new, scientifically distinct material systems.

---

**Scenario ID:** UAT-C02-003
**Priority:** Medium
**Title:** Verify Intelligent Data Sampling with Farthest Point Sampling (FPS)

**Description (Min 300 words):**
The purpose of this scenario is to test the effectiveness of the new Farthest Point Sampling (FPS) algorithm. Random sampling, while simple, can lead to a final dataset where many structures are very similar to each other. FPS is designed to produce a maximally diverse subset, which is often more efficient for training MLIPs. The user will perform an experiment to compare the two methods.

First, the user will run a single, relatively long exploration for an alloy system to generate a large trajectory file. This file will serve as the common input for both samplers. Next, the user will run the `sampling` and `storage` stages twice. In the first run, they will configure the `RandomSampler` to select 100 structures and save them to `random.db`. In the second run, they will configure the `FarthestPointSampler` to select 100 structures from the *same* trajectory and save them to `fps.db`. The user will then use visualization tools within the notebook to compare the two datasets. For instance, they could plot the energy distribution of the structures in each database. A successful FPS run is expected to produce a much broader energy distribution, indicating that it has selected a wider variety of low- and high-energy states. This scenario validates that the FPS implementation provides a tangible benefit over the simpler baseline, leading to higher-quality datasets.

---

**Scenario ID:** UAT-C02-004
**Priority:** High
**Title:** Basic Workflow Execution via Web User Interface

**Description (Min 300 words):**
This scenario provides a manual test for the new Web UI. The goal is to verify that a non-expert user can successfully configure and launch a pipeline run and view the results without ever touching the command line or a YAML file. The user will first launch the web application. They will be presented with a graphical interface containing input fields, sliders, and dropdown menus for the key configuration parameters.

The user will fill out the form to replicate the simple CuAu alloy configuration from the Cycle 1 UAT. They will select the `AlloyGenerator`, input the elements, choose the `MDEngine`, and set the number of steps. Once the form is complete, the user will click the "Launch Pipeline" button. The UI should provide feedback that the job has started. After the job completes, the user will navigate to a "Results" tab or page in the UI. This page should allow them to select the output database. Upon selection, it should display a summary of the results, such as the number of structures generated, and present a 3D visualization of a few of the final atomic structures. The test is successful if the entire workflow can be completed through the browser and the final results are consistent with those produced by the CLI in Cycle 1.

## 2. Behavior Definitions

The following definitions describe the expected system behavior in a Gherkin-style format.

---

**Feature:** Advanced Material Generation and Sampling

**As a** Computational Chemist
**I want to** use advanced simulation and sampling techniques for various material types
**So that** I can generate more diverse and higher-quality datasets for my MLIPs.

---

**Scenario:** Generating a dataset for an ionic crystal

**GIVEN** I have a configuration file specifying:
  - A generator for an `ionic` system with cation `Na` and anion `Cl`.
  - An exploration phase using the hybrid `MD/MC` engine.
  - A sampler using the `FPS` algorithm.
  - An output database file named `NaCl_dataset.db`.

**WHEN** I execute the main pipeline command.

**THEN** the command should complete successfully.
**AND** a database file named `NaCl_dataset.db` must be created.
**AND** every structure in the database must be charge-neutral (equal numbers of Na and Cl atoms).
**AND** an analysis of the structure's nearest neighbors should show a strong preference for Na-Cl pairs.

---

**Feature:** Web-Based User Interface

**As a** Researcher
**I want to** use a graphical interface to configure and run the pipeline
**So that** I can generate datasets without needing to write YAML files or use the command line.

---

**Scenario:** Launching a job from the Web UI

**GIVEN** the Web UI application is running and I have opened it in my browser.

**WHEN** I fill the web form with parameters for a simple `CuAu` alloy.
**AND** I click the "Launch Pipeline" button.

**THEN** the UI should indicate that a job is running in the background.
**AND** after some time, the UI should indicate the job has completed successfully.
**AND** I should be able to navigate to a results page.
**AND** the results page should allow me to view the 3D structures generated by the pipeline I just configured.
