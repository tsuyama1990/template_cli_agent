# User Acceptance Testing (UAT): Cycle 2

This document outlines the User Acceptance Testing scenarios for the deliverables of Cycle 2. The focus of this cycle is on the advanced scientific features and the new web-based user interface, which together represent a major leap in the capabilities and usability of the MLIP-AutoPipe tool. These UATs are designed to confirm that the new capabilities not only work correctly from a technical standpoint but also provide a tangible benefit to the end-user and a smooth, intuitive experience. The tests will validate the scientific correctness of the new exploration and sampling methods, the core functionality of the web UI, and the performance improvements from parallelization.

## 1. Test Scenarios

These scenarios cover the major new features introduced in this cycle: Farseed Point Sampling (FPS), hybrid MD/MC exploration, the Web UI, and parallel execution. Each test is designed from the perspective of a user trying to achieve a specific scientific or operational goal.

| Scenario ID | Title                                       | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C2-01   | Generate a Diverse Dataset using FPS        | High     |
| UAT-C2-02   | Run Pipeline via Interactive Web UI         | High     |
| UAT-C2-03   | Enhanced Exploration with Hybrid MD/MC      | Medium   |
| UAT-C2-04   | Verify Parallel Execution Speed-up          | Medium   |

### Scenario UAT-C2-01: Generate a Diverse Dataset using FPS

**Description:** This UAT is designed to validate the core scientific premise of the new Farthest Point Sampling (FPS) method. The central goal is to demonstrate to the user that choosing FPS results in a dataset that is measurably more diverse than one generated using the baseline random sampling from Cycle 1. This is a critical test of the feature's value, as the primary reason for implementing FPS is to improve the quality of the final dataset, which should lead to better-trained MLIP models. The test requires the user to generate two comparable datasets, one using random sampling and one using FPS, and then use a provided analysis tool (such as a Jupyter Notebook) to quantify and compare their diversity. This comparison is key; it provides an objective measure of success. The user experience here is about empowering the researcher with a more advanced scientific tool and giving them the confidence that it works as advertised. The successful completion of this test will confirm that the complex underlying machinery of SOAP descriptors and the FPS algorithm is correctly integrated and delivering its intended benefit, justifying its inclusion in the tool. It moves the application from being a simple automation script to an intelligent data generation framework.

**User Story:** As a researcher, I want to use an advanced sampling method like FPS to automatically select the most structurally unique configurations from a long simulation trajectory, ensuring my final dataset is as diverse as possible, so that I can train a more robust and accurate MLIP model with the same amount of data.

**Execution Steps:**
1.  The user creates two configuration files for the same physical system, `config_random.yaml` and `config_fps.yaml`. They are identical in every respect except for the `sampling.method` parameter, which is set to `'random'` and `'fps'` respectively. Both are configured to sample 10 final structures.
2.  The user runs the pipeline twice, once for each configuration, generating `random_dataset.db` and `fps_dataset.db`.
3.  The user is provided with a simple Jupyter Notebook (`analyse_diversity.ipynb`). This notebook contains Python code to load both databases, compute a structural fingerprint (e.g., average SOAP vector) for each structure, and then calculate an average pairwise similarity score for each dataset.
4.  The user executes the notebook. The output should clearly show that the average similarity score for `fps_dataset.db` is significantly lower than the score for `random_dataset.db`, providing quantitative proof that FPS has selected a more diverse set of structures.

### Scenario UAT-C2-02: Run Pipeline via Interactive Web UI

**Description:** This UAT validates the functionality, usability, and overall user experience of the new web interface. The primary goal is to ensure that a user, potentially one with limited command-line experience, can successfully configure and launch a pipeline run entirely from their web browser. This test covers the full user journey within the UI: starting the application, navigating the interface, filling in parameters using interactive widgets, initiating the run, monitoring its progress in real-time, and retrieving the results. The experience should be intuitive, requiring no prior knowledge of YAML or CLI syntax. The UI should provide clear labels, sensible default values, and helpful tooltips for the various parameters. A crucial aspect of this test is the real-time feedback mechanism; the user must be able to see the log output from the backend process directly in the browser, which provides reassurance that the task is running and allows for immediate diagnosis of any issues. A successful test will demonstrate that the web UI is a viable and user-friendly alternative to the CLI, making the tool accessible to a much broader audience and lowering the barrier to entry for complex materials science simulations.

**User Story:** As a scientist who is not an expert in command-line tools, I want to use a simple and intuitive graphical interface to configure my data generation pipeline, run it with a click of a button, and see the results, so that I can leverage the power of the tool without needing to learn how to write YAML files or operate in a terminal.

**Execution Steps:**
1.  The user starts the web application by executing a simple command in their terminal: `mlip-autopipec run-web-ui`.
2.  The user opens the local URL provided by the command (e.g., `http://localhost:8501`) in their web browser.
3.  The user is presented with a clean, well-structured form. They interact with the form widgets, such as typing element names into text boxes, moving sliders to set temperature, and selecting 'fps' from a dropdown menu.
4.  After filling out the form, the user clicks the prominent "Run Pipeline" button.
5.  Immediately, a log output area appears on the page and begins streaming the real-time output from the backend pipeline, exactly as it would appear in the terminal.
6.  The user waits for the pipeline to complete, observing the log messages. Upon successful completion, the UI displays a clear "Success!" message and a "Download Database" button becomes active.
7.  The user clicks the download button and their browser successfully downloads the resulting database file (e.g., `mlip_data.db`).

## 2. Behavior Definitions

These Gherkin-style definitions formalize the expected behavior for the UAT scenarios of Cycle 2, providing a detailed, unambiguous contract for the new features. They are written to be comprehensive, covering not only the technical correctness of the features but also the nuances of the user experience, such as the clarity of feedback and the scientific validity of the results. This ensures that the implementation fully aligns with the goals of making the tool more powerful, more intelligent, and easier to use.

**Scenario: Generate a Diverse Dataset using FPS**
```gherkin
GIVEN the user has two identical configuration files, "config_fps.yaml" and "config_random.yaml"
  AND "config_fps.yaml" has the 'sampling.method' parameter set to 'fps'
  AND "config_random.yaml" has the 'sampling.method' parameter set to 'random'
  AND both files request the same number of final samples (e.g., 20) from the same initial system
WHEN the user runs the CLI pipeline for both of these configurations
  AND two output databases, "fps_dataset.db" and "random_dataset.db", are successfully generated
  AND the user then executes a provided Jupyter Notebook designed to analyze and compare the diversity of the two databases
THEN the notebook must run without errors, successfully loading and processing both database files
  AND the final analysis output of the notebook must show a quantitative result
  AND this result must clearly indicate that the average structural similarity within "fps_dataset.db" is significantly lower than the average similarity within "random_dataset.db".
```

**Scenario: Run Pipeline via Interactive Web UI**
```gherkin
GIVEN the user has started the web UI application via the 'mlip-autopipec run-web-ui' command
  AND their web browser is open to the application's local URL
WHEN the user interacts with the web form to define a complete pipeline configuration, setting parameters for the system, exploration, and sampling stages
  AND the user clicks the "Run Pipeline" button on the web page
THEN the web page must not become unresponsive and should immediately show a visual indicator that the process has started
  AND a text area must appear on the page, displaying the live log output from the backend CLI process as it executes
  AND the log messages should show the pipeline progressing through the Generation, Exploration, Sampling, and Storage stages
  AND once the backend process is complete, the web page must display a clear and prominent success message
  AND a "Download Database" button must become visible and clickable
  AND when the user clicks the download button, the browser must initiate a download of the final, valid database file.
```

**Scenario: Enhanced Exploration with Hybrid MD/MC**
```gherkin
GIVEN the user has a starting structure file, "ordered_alloy.xyz", representing a chemically ordered alloy
  AND the user has two configuration files, "config_md.yaml" and "config_md_mc.yaml", both pointing to this starting structure
  AND "config_md.yaml" has the 'mc_config.enabled' parameter set to false
  AND "config_md_mc.yaml" has 'mc_config.enabled' set to true, with a specified 'swap_frequency'
WHEN the user runs the pipeline for both configurations to generate two different trajectories
  AND the user uses an external visualization tool (e.g., 'ase gui') to view the final structure from both the pure MD and the hybrid MD/MC runs
THEN the final structure from the pure MD run should visually retain a high degree of its initial chemical order
  AND the final structure from the hybrid MD/MC run should show significant chemical disorder, with atoms of different species clearly mixed throughout the crystal lattice, providing qualitative evidence that the Monte Carlo swaps were successfully performed.
```

**Scenario: Verify Parallel Execution Speed-up**
```gherkin
GIVEN a configuration file that specifies the generation of a non-trivial number of seed structures (e.g., 8)
  AND this requires running 8 independent exploration simulations
WHEN the user first runs the pipeline with the 'exploration.num_parallel' parameter set to 1
  AND the user records the total time taken for the "Exploration" stage as reported in the application's log output (Time T1)
  AND the user then runs the exact same pipeline but with 'exploration.num_parallel' set to a value greater than 1 that matches their available CPU cores (e.g., 4)
  AND the user records the new total time for the "Exploration" stage (Time T4)
THEN the recorded time T4 must be substantially less than the recorded time T1, demonstrating a clear and significant performance improvement.
```
