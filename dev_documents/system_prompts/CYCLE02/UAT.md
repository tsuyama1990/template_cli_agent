# User Acceptance Testing (UAT): Cycle 2 - Advanced Features and Web UI

This document outlines the User Acceptance Testing (UAT) plan for the second development cycle of the MLIP-AutoPipe project. This cycle is pivotal as it introduces the advanced scientific features and the graphical user interface (Web UI) that transform the project from a specialist's command-line tool into a more accessible and powerful platform. The UAT will focus on two key areas: first, verifying that the new, more complex exploration and sampling features produce physically meaningful and superior results compared to the baseline from Cycle 1. Second, it will rigorously assess the Web UI to ensure it provides a smooth, intuitive, and complete user experience, from initial configuration to final results visualization. The scenarios are designed from the perspective of a materials scientist who may not be an expert in command-line tools but has a deep understanding of the scientific problem. A combination of interactive browser-based testing and data analysis in a provided Jupyter Notebook (`UAT_Cycle2.ipynb`) will be used to comprehensively verify the outcomes and ensure that the new features deliver on their promised value.

## 1. Test Scenarios

These scenarios are designed to allow a user to perform a thorough validation of the new functionalities introduced in Cycle 2. They cover the end-to-end user experience of the Web UI and the scientific validity of the new backend features.

| Scenario ID | Priority | Description |
| :--- | :--- | :--- |
| **UAT-C2-001** | High | **Web UI Golden Path End-to-End Run**: This is the most critical scenario for this cycle. It requires the user to perform a complete, end-to-end workflow for a standard alloy system, but using *only* the new Web UI. This test will verify everything from graphical configuration to launching the run, monitoring its progress, and interactively visualizing the final results. Its success is the primary criterion for the UI's release readiness. |
| **UAT-C2-002** | High | **Scientific Verification of Automatic Ensemble Switching**: This test is focused on a key piece of scientific intelligence. It will verify the system's ability to correctly identify a surface/slab system, automatically switch to the appropriate NVT ensemble, and thus prevent the unphysical cell deformation that would otherwise occur. This confirms that the system can be trusted to produce physically realistic simulation results for non-bulk systems. |
| **UAT-C2-003** | Medium | **Scientific Verification of Hybrid MD/MC (Atom Swap)**: This scenario is designed to confirm the scientific value of the new hybrid MD/MC feature. The user will compare the results of a run with and without atom swaps enabled, verifying that the hybrid method generates a more diverse set of chemical orderings in a binary alloy, which is crucial for training robust potentials for such systems. |
| **UAT-C2-004** | High | **Scientific Verification of Farthest Point Sampling (FPS)**: This test validates the newly implemented intelligent sampling feature. The user will compare a dataset generated using FPS against one generated with random sampling. The goal is to quantitatively and qualitatively verify that the FPS dataset is more structurally diverse, demonstrating its superiority for creating efficient and comprehensive training sets. |
| **UAT-C2-005** | Medium | **Interactive Exploration and Visualization of Results**: This scenario focuses on the usability of the results page in the Web UI. It ensures that the user can effectively browse, search, filter, and visualize the generated structures after a run is complete, confirming that the UI is not just a launcher but also a powerful tool for data analysis and exploration. |

---

### **Scenario UAT-C2-001: Web UI Golden Path**

**(Min 300 words)**
This scenario provides a comprehensive test of the primary user workflow for the new Web UI, representing the "golden path" from start to finish. The entire test must be performed without any interaction with the command line. The user's goal is to create a small but complete dataset for a Copper-Gold (CuAu) alloy, a classic system for studying chemical ordering. The user will begin by opening their web browser and navigating to the web application's home page. They will be greeted by a clean, organized, and well-documented web form. They will proceed to fill out the form fields to define the system, entering "Cu" and "Au" as the elements, selecting a 50/50 composition, and choosing an FCC lattice. In the "Exploration" section of the form, they will specifically enable the new advanced features by ticking the checkboxes for "Enable Hybrid MC" and "Use ZBL Repulsion". In the "Sampling" section, they will select "FPS" from the dropdown menu.

Once the configuration is complete, the user will click the "Launch Run" button. The UI must provide immediate feedback, transitioning to a "Run in Progress" page. This page will display a unique run ID for reference and, crucially, will feature a log viewer that shows the real-time output from the backend pipeline. This allows the user to monitor the progress of the four stages. The user will wait for the run to complete. Upon successful completion, the page should automatically update to show a "Completed" status and provide a prominent, clickable link to the results page. On the results page, the user will first see a summary of the run's configuration. Below this, they will find a paginated table of the generated structures. The user will test the interactivity by sorting the table by energy and then clicking on the lowest-energy structure. This action should bring up an integrated 3D viewer that loads and displays a plausible, ordered L1<sub>0</sub> CuAu structure. This scenario is considered fully passed if the user can successfully and intuitively navigate this entire workflow without encountering any errors or confusing UI elements, and the final data is generated and displayed exactly as expected.

### **Scenario UAT-C2-002: Verify Automatic Ensemble Switching**

**(Min 300 words)**
This scenario is designed to test a crucial piece of scientific intelligence that has been built into the new exploration engine: the ability to automatically select the correct thermodynamic ensemble. In computational materials science, applying a constant pressure (NPT) simulation to a system that contains a vacuum, such as a surface slab, is a common but serious mistake. It can lead to the unphysical collapse of the vacuum layer as the simulation box shrinks to reach the target pressure, rendering the simulation results meaningless. This test will verify that the automatic ensemble switching feature robustly prevents this error. The user will configure a run to generate a dataset for a Platinum (111) surface. This will be done by using a generator that starts with a bulk Pt crystal, cleaves the (111) surface, and adds a 15 Å vacuum layer to separate the slab from its periodic images.

Critically, in the "Exploration" section of the configuration (either in the UI or a YAML file), the user will leave the "Ensemble" choice set to "auto" or will leave the field blank, relying on the system's default automatic detection. They will then launch the run. The first point of verification will be the log file, which should contain a clear, explicit message stating that a vacuum slab was detected and that the simulation ensemble has been automatically set to NVT. After the run completes, the user will perform a more rigorous scientific validation using the provided Jupyter notebook, `UAT_Cycle2.ipynb`. The notebook will load the main trajectory file from the simulation. The user will then execute a code block that plots the length of the `c` lattice vector (the dimension perpendicular to the slab) as a function of the simulation time. The expected and correct result is a flat line, indicating that the `c` vector remained constant throughout the simulation. This proves that the NVT ensemble was correctly chosen and the vacuum layer was preserved. For a compelling comparison, the user could then be guided to run a second, shorter simulation where they manually force the ensemble to NPT, and the resulting plot would show the `c` vector shrinking over time, visually demonstrating the very problem that the automatic feature successfully prevents. The UAT is passed if the default, automatic behavior correctly identifies the slab and preserves the vacuum layer's dimension.

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behavior for each test scenario in a structured, unambiguous format.

---

**Scenario: UAT-C2-001 - Web UI Golden Path**

```gherkin
GIVEN the user has successfully opened the MLIP-AutoPipe Web UI in their web browser.

WHEN the user interacts with the web form to input the following parameters:
  - System: A Copper-Gold (CuAu) alloy with a 50/50 composition.
  - Exploration: The "Enable Hybrid MC" feature is turned ON.
  - Sampling: The "FPS" method is selected to produce 20 final structures.
AND the user clicks the main "Launch Run" button.

THEN the user's browser should immediately navigate to a new status page.
AND this status page must display a unique run ID and show the run's current status as "In Progress".
AND after a period of time, the run should complete successfully, and the status on the page should update to "Completed".
AND the results page should contain a link or table displaying information for exactly 20 curated structures.
AND when the user clicks on any structure in the results table, an interactive 3D visualization of that CuAu structure must be rendered on the page.
```

---

**Scenario: UAT-C2-002 - Verify Automatic Ensemble Switching**

```gherkin
GIVEN the user has configured a pipeline run for a Platinum (111) surface slab which includes a 15 Å vacuum layer.
AND the configuration for the exploration stage explicitly leaves the 'ensemble' parameter unset to rely on the default automatic detection mechanism.

WHEN the user launches this pipeline run.

THEN the application's log file (visible in the UI or on disk) must contain a specific log message indicating that a vacuum slab was detected and that the simulation ensemble was consequently set to NVT.
AND the MD simulation must complete the run successfully without any errors.
AND a post-run analysis of the output trajectory file, performed using the provided Jupyter Notebook, must show that the simulation cell dimension perpendicular to the slab surface remains constant (within small thermal fluctuations) throughout the entire run.
```

---

**Scenario: UAT-C2-003 - Verify Hybrid MD/MC (Atom Swap)**

```gherkin
GIVEN the user has prepared two separate configurations for a Nickel-Iron (NiFe) alloy, which are identical in every way except:
  - `config_md.yml` has `enable_hybrid_mc: false`.
  - `config_mc.yml` has `enable_hybrid_mc: true` and specifies a non-zero swap frequency.

WHEN the user runs the pipeline twice, once with each of these two configurations.

THEN two distinct databases, "md_run.db" and "mc_run.db", must be successfully created.
AND the user opens the provided analysis Jupyter Notebook, `UAT_Cycle2.ipynb`.
AND the notebook is used to calculate a chemical short-range order parameter for every structure in both databases.
AND a histogram plotting the distribution of this order parameter must show a significantly wider and more varied distribution for the "mc_run.db" dataset when compared to the narrower distribution from the "md_run.db" dataset.
AND this result provides clear, quantitative evidence that the atom swap moves have enabled the exploration of a more diverse range of chemical configurations.
```

---

**Scenario: UAT-C2-004 - Verification of Farthest Point Sampling (FPS)**

```gherkin
GIVEN the user has prepared two configurations for the same chemical system, identical except for the sampling method:
  - One configuration specifies `method: 'random'`.
  - The other configuration specifies `method: 'fps'`.

WHEN the user executes two separate runs, one for each configuration.

THEN two distinct databases, "random_dataset.db" and "fps_dataset.db", must be created.
AND the user opens the provided Jupyter Notebook, `UAT_Cycle2.ipynb`.
AND the notebook is used to compute a structural similarity metric (e.g., based on SOAP kernels) for all pairs of structures within each database.
AND the analysis must quantitatively show that the average similarity score between structures in the "fps_dataset.db" is lower than the average similarity in the "random_dataset.db".
AND a visual plot of the similarity matrices for both databases must qualitatively confirm that the FPS-sampled set has lower off-diagonal correlations, indicating higher structural diversity.
```

---

**Scenario: UAT-C2-005 - Interactive Visualization of Results**

```gherkin
GIVEN a pipeline run has completed successfully and the user has navigated to its results page in the Web UI.

WHEN the user interacts with the results table on the page.

THEN the table must be paginated, allowing the user to easily navigate through a large set of structures (e.g., pages 1, 2, 3...).
AND the table's columns, such as 'Structure ID' and 'Potential Energy', must be sortable by clicking on the headers.
AND the user must be able to type a filter query, such as "Energy < -500 eV", into a filter box.
AND upon entering the filter query, the table must dynamically update to show only the structures that satisfy this criterion.
AND when the user clicks on any structure row in the filtered or unfiltered table, the 3D viewer component on the page must load and display the correct corresponding atomic structure.
```
