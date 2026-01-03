# User Acceptance Testing (UAT): CYCLE02

This document outlines the User Acceptance Testing scenarios for the advanced features implemented in CYCLE02 of the MLIP-AutoPipe project. The goal of this UAT plan is to provide a clear, user-centric framework for verifying that the new capabilities—the hybrid exploration engine, Farthest Point Sampling, ionic crystal support, and the Web UI—are not only functionally correct but also effective, usable, and provide a clear, demonstrable benefit to the end-user. These tests are designed to answer the question: "Does the new feature work, and does it help me achieve my scientific goals better than before?"

## 1. Test Scenarios

The following scenarios have been designed to provide comprehensive coverage of the new features in CYCLE02. They focus on validating the scientific effectiveness of the new algorithms and the usability of the new interfaces.

| Scenario ID | Title                                               | Priority | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|-------------|-----------------------------------------------------|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| UAT-C02-001 | Validate Hybrid MD/MC Engine for Alloy Exploration  | High     | This scenario is the primary test of the most significant new technical feature in this cycle. The user will run the pipeline for a binary alloy system that is known to benefit from atomic swaps, such as a disordered solid solution or a material with multiple ordered phases (e.g., SiGe or CuAu). The core goal is to verify that the hybrid engine explores the material's configuration space more effectively and efficiently than a standard MD simulation. Success is not just about the pipeline completing, but about finding clear, physical evidence in the output that the Monte Carlo moves are active and are successfully creating new, diverse atomic arrangements that would be inaccessible to a simple MD run in a similar amount of time. This test confirms the scientific validity of the new explorer. |
| UAT-C02-002 | Compare FPS and Random Sampling for Diversity       | High     | This scenario is designed to directly and quantitatively demonstrate the value of the new Farthest Point Sampling (FPS) algorithm. It provides a head-to-head comparison against the baseline random sampler from CYCLE01. The user will first generate a large pool of candidate structures by running the exploration stage. Then, they will use this same pool of structures as input for two separate sampling runs: one using the `random` sampler and another using the `fps` sampler. The user will then perform a data-driven analysis to compare the diversity of the two resulting datasets. The explicit goal is to confirm, with quantitative evidence, that FPS produces a more structurally varied and less redundant set of structures, thereby proving its superiority for creating efficient, high-quality training sets. |
| UAT-C02-003 | Generate a Dataset for an Ionic Crystal (NaCl)      | Medium   | This scenario validates the expanded material support of the pipeline, ensuring it is not just a tool for alloys. The user will configure and run the entire pipeline from start to finish to generate a dataset for a simple, classic ionic crystal: Sodium Chloride (NaCl). The test will verify that the new `IonicGenerator` can correctly create charge-neutral initial structures based on a known crystal prototype. Success is defined by the creation of a valid final database where all the structures within it are physically plausible for an ionic system, demonstrating that the entire pipeline can now handle this completely different class of materials. |
| UAT-C02-004 | End-to-End Run and Visualization via Web UI         | High     | This scenario tests the usability, functionality, and completeness of the new Web UI. The user will not use the command line at all for this test. Instead, they will use the browser-based interface to configure all the necessary parameters for a standard pipeline run (e.g., the binary alloy from CYCLE01's UAT), launch the job, monitor its progress, and inspect the results. The goal is to verify that the UI provides a complete, intuitive, and user-friendly alternative to the CLI for the entire workflow. A key part of this test is to confirm that the visualization component of the UI works correctly, allowing the user to view the generated structures directly in their browser. |

**UAT Approach using a Jupyter Notebook:**

To ensure these complex and sometimes data-intensive tests are easy for a user to perform and interpret, a `CYCLE02_UAT.ipynb` Jupyter Notebook will be the primary vehicle for this UAT. This interactive notebook will be especially critical for scenarios like `UAT-C02-002`, where quantitative analysis is required to prove the benefit of the new feature. The notebook will be meticulously structured to guide the user through each scenario step-by-step. It will contain:
1.  Code cells to programmatically create and save the distinct Hydra configuration files needed for each test run.
2.  Cells that execute the pipeline runs using shell commands, capturing the output logs for inspection.
3.  For `UAT-C02-002`, the notebook will contain a detailed analysis section. This section will include Python code to load the two generated databases, calculate SOAP descriptors for every structure, perform a dimensionality reduction (e.g., with PCA), and then create plots that visually show the distribution of the structures in descriptor space. This will provide a clear, visual demonstration of the greater "volume" of space covered by the FPS samples compared to the random samples. It will also calculate and display a quantitative diversity metric to support the visual evidence.
4.  For the Web UI test, the notebook will first provide clear instructions on how to launch the UI with a single command. It will then describe the steps the user should follow within the UI. Finally, it will contain a Python cell that programmatically checks the file system for the expected output from the UI-initiated run, thus providing an automated way to confirm its successful completion.

## 2. Behavior Definitions

The following Gherkin-style behavior definitions describe the expected outcomes for each high-priority test scenario in detail.

---

### **Scenario: UAT-C02-001 - Validate Hybrid MD/MC Engine for Alloy Exploration**

This test verifies that the hybrid engine is correctly performing its core function: Monte Carlo atom swaps.

**GIVEN**
The user has created a valid Hydra configuration for a Silicon-Germanium (SiGe) binary alloy.
And the `explorer` section is configured to use the `hybrid_md_mc` engine.
And the configuration specifies a reasonably high `swap_frequency` (e.g., every 10 MD steps) and a simulation temperature (e.g., 1000K) high enough to make swaps thermodynamically plausible.

**WHEN**
The user executes the full pipeline using the CLI with this configuration.
And the exploration stage completes successfully without any errors related to the MC moves.

**THEN**
A final database should be created successfully.
And when the user analyzes the intermediate trajectory files, there must be clear and unambiguous evidence of atomic swaps having occurred.
And this can be verified by tracking a specific atom: for instance, a Germanium atom that was initially in a local environment surrounded only by Silicon atoms should, at a later timestep in the trajectory, be found adjacent to another Germanium atom. This provides strong, undeniable proof that the MC swaps are active and are effectively altering the material's structure beyond what simple atomic vibration from MD would achieve.

---

### **Scenario: UAT-C02-002 - Compare FPS and Random Sampling for Diversity**

This test provides a clear, quantitative validation of the superiority of Farthest Point Sampling over random sampling.

**GIVEN**
A common set of simulation trajectory files has already been generated by a completed exploration run, representing a large pool of candidate structures.

**WHEN**
The user first runs only the sampling and storage stages of the pipeline with the `sampler` configured to use the `random` implementation, creating a database named `random_results.db`.
And the user then runs the sampling and storage stages a second time, starting from the *exact same* trajectory files, but this time with the `sampler` configured to use the `fps` implementation, creating a database named `fps_results.db`.

**THEN**
Two separate databases, `random_results.db` and `fps_results.db`, must be created, each containing the same number of sampled structures.
And when the user performs a quantitative comparison of the structures in both databases (e.g., by executing the analysis in the provided UAT Jupyter notebook), the results must show that the database generated by FPS has a significantly higher diversity score.
And this will be demonstrated visually by a plot showing the FPS samples are more spread out in descriptor space.
And it will be confirmed numerically by a metric (like the standard deviation of pairwise distances) that is higher for the FPS set, proving that FPS is a more effective strategy for creating a diverse dataset.

---

### **Scenario: UAT-C02-004 - End-to-End Run and Visualization via Web UI**

This test ensures that the Web UI is a complete, functional, and user-friendly interface for the entire pipeline.

**GIVEN**
The user has launched the MLIP-AutoPipe Web UI application using the `mlip_autopipec gui` command, and it is accessible in their browser.

**WHEN**
The user navigates to the Web UI in their browser.
And the user interacts with the various input forms, sliders, and dropdowns to visually configure all the parameters for a simple binary alloy system (similar to the test case from `UAT-C01-001`).
And the user clicks the "Run Pipeline" button in the UI.

**THEN**
The UI should provide immediate visual feedback that the pipeline process has been successfully started in the background.
And after some time, the UI's status display should indicate that the run has completed successfully.
And when the user checks the file system, the corresponding output files, including the final ASE database, must be present in the expected location.
And the contents of the database must be correct and consistent with the parameters that were set in the Web UI.
And the user should be able to navigate to a "Results Visualization" tab or section in the UI, which should display a 3D visualization of one or more of the structures from the final database, confirming that the entire workflow, from configuration to visualization, can be accomplished through the UI.
