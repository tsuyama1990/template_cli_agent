# CYCLE 02: USER ACCEPTANCE TESTING (UAT) - Advanced Features and Web UI

This document outlines the User Acceptance Testing (UAT) scenarios for the advanced features and the new Web User Interface delivered in Cycle 2 of the MLIP-AutoPipe project. The goal is to ensure these new, complex features are not only functional but also provide a tangible benefit and a positive experience for the user.

## 1. Test Scenarios

| Scenario ID | Title                                          | Priority |
| :---------- | :--------------------------------------------- | :------- |
| UAT-C2-001  | Interactive Pipeline Run via Web UI            | High     |
| UAT-C2-002  | Verifying the Effect of Farthest Point Sampling (FPS) | High     |
| UAT-C2-003  | Generating an Ionic Crystal Structure        | Medium   |
| UAT-C2-004  | Monitoring a Running Job from the Web UI       | High     |

### Scenario UAT-C2-001: Interactive Pipeline Run via Web UI

**(Min 300 words)**
This is the primary "happy path" scenario for the new graphical interface. It validates the entire user journey, from opening the web page to successfully generating and visualizing a dataset. The user experience should be seamless and intuitive, abstracting away the command-line complexity. This UAT will be performed directly in a web browser, with the server running locally.

The user will be guided through the following steps, perhaps by a short tutorial video or a markdown guide:
1.  **Launch:** The user starts the web server by running a simple command like `mlip-autopipec web-ui`.
2.  **Navigation:** The user opens their web browser and navigates to the provided local URL (e.g., `http://127.0.0.1:8000`).
3.  **Configuration:** The user is presented with a clean, form-based interface. They will interactively fill in the parameters for a simple binary alloy (e.g., SiGe). The form will have clear labels, sensible defaults, and perhaps tooltips explaining what each parameter does. For example, they will select "alloy" as the system type, enter `Si` and `Ge` as elements, and specify a composition.
4.  **Execution:** The user clicks a prominent "Run Pipeline" button. The UI should immediately respond, indicating that the job has been submitted and is now running. A unique job ID should be displayed.
5.  **Completion & Visualization:** Once the job completes (which will be quick for this test case), a results section will become visible. This section will feature an embedded 3D molecular viewer. The generated structures (e.g., the 5-10 structures sampled from the run) will be automatically loaded and displayed. The user can then rotate, zoom, and inspect each atomic structure individually directly within the browser. This provides immediate, visual confirmation that the workflow was successful and the output is physically sensible. This instant visual feedback is a core part of the amazing user experience this UAT is designed to validate.

### Scenario UAT-C2-002: Verifying the Effect of Farthest Point Sampling (FPS)

**(Min 300 words)**
This UAT is designed to demonstrate the scientific value of one of the key new features: Farthest Point Sampling. It aims to provide the user with a clear, visual understanding of *why* FPS is superior to random sampling for creating diverse datasets. This will be presented as a comparative Jupyter Notebook, `UAT-C2-002.ipynb`, which will be a powerful educational tool.

The notebook will guide the user through this comparison:
1.  **Setup:** A baseline dataset is needed. The notebook will first run a longer, more complex MD simulation (the "Exploration" phase) on a system and save the entire trajectory, which contains thousands of highly correlated structures.
2.  **Run A: Random Sampling:** The notebook will then programmatically create a configuration to run only the "Sampling" part of the pipeline on this trajectory data, using the `random` method to select 50 structures.
3.  **Run B: FPS Sampling:** Next, it will create another configuration to run the sampling, but this time using the `fps` method to select 50 structures.
4.  **Analysis and Visualization:** This is the core of the UAT. The notebook will perform a dimensionality reduction technique (like PCA or t-SNE) on the SOAP descriptors of the *entire* trajectory, creating a 2D map of the "structure space".
    *   First, it will plot all the points from the full trajectory as a grey background, showing the overall shape of the explored space.
    *   Then, it will overlay the 50 points selected by random sampling in blue. The user will visually observe that these points are likely clumped together in the most densely populated regions of the map.
    *   Finally, it will overlay the 50 points selected by FPS in red. The user will see a dramatic difference: the red points will be spread out much more evenly across the entire map, covering the space far more effectively and including outliers. This visualization provides an immediate, intuitive "aha!" moment, proving the value of the advanced sampling feature.

### Scenario UAT-C2-003: Generating an Ionic Crystal Structure

**(Min 300 words)**
This scenario validates the expanded capabilities of the structure generator module. It confirms that the system can now handle material classes beyond simple alloys, specifically ionic crystals like Sodium Chloride (NaCl), which require adherence to charge-balance rules. This demonstrates the framework's extensibility. The test will be conducted via the Web UI to reinforce the unified user experience.

1.  **Configuration:** The user will navigate to the Web UI. In the configuration form, they will select `ionic` as the System Type.
2.  **Input:** The form will dynamically update to show fields relevant to ionic crystals. The user will input `Na` and `Cl` as the elements. The system might have an input for oxidation states (e.g., `Na: 1`, `Cl: -1`) or infer them. The user will request a standard rock salt crystal structure.
3.  **Execution:** The user will click "Run Pipeline".
4.  **Verification:** The primary verification is visual. When the results are displayed in the 3D viewer, the user will see the classic, alternating checkerboard pattern of Na and Cl ions characteristic of a rock salt structure. The notebook accompanying this UAT, `UAT-C2-003.ipynb`, will provide a more rigorous check. It will connect to the output database, retrieve the generated structure, and use a chemical analysis library (like `pymatgen`) to calculate the electrostatic energy or formal oxidation states of the structure, asserting that the generated crystal is charge-neutral. This confirms that the underlying physical constraints have been correctly implemented and satisfied.

### Scenario UAT-C2-004: Monitoring a Running Job from the Web UI

**(Min 300 words)**
This UAT focuses on the user experience for long-running jobs, where clear feedback on the status is essential. It tests the asynchronous nature of the web backend and the frontend's ability to poll for and display progress updates.

1.  **Setup:** The user will configure a job in the Web UI that is designed to take a noticeable amount of time (e.g., 30-60 seconds). This could involve requesting a larger number of initial structures or a longer exploration phase.
2.  **Execution:** The user clicks "Run Pipeline".
3.  **Monitoring:** The user will observe the status display area on the web page. The UI should not freeze or become unresponsive. Instead, it should provide a series of updates, demonstrating that it is actively polling the backend server. The user should see status messages like:
    *   "Job submitted. Status: Pending"
    *   "Status: Generating initial structures..."
    *   "Status: Exploring structures (5 / 20 complete)..."
    *   "Status: Sampling final structures..."
    *   "Status: Finished."
4.  **Verification:** The user verifies that the progress updates are logical and correspond to the actual stages of the pipeline. The progress counter (e.g., "5 / 20") should increment as the job runs. The final "Finished" status should coincide with the results appearing in the 3D viewer. This confirms that the entire monitoring system is working correctly, providing the user with the transparency and peace of mind needed when running computationally expensive tasks.

## 2. Behavior Definitions

**GIVEN** the user has the web server running and has opened the application in their browser
**AND** they have filled out the web form to configure a 'SiGe' alloy with 5 initial structures
**WHEN** they click the "Run Pipeline" button
**THEN** a job ID is immediately displayed
**AND** the status text shows "Pending" or "Running"
**AND** after a short time, the status text shows "Finished"
**AND** a 3D viewer appears, populated with 5 atomic structures
**AND** inspecting the structures in the viewer shows they contain both 'Si' and 'Ge' atoms.

---

**GIVEN** a user has run two sampling processes, one 'random' and one 'fps', on the same large trajectory file
**WHEN** they view the comparative analysis in the `UAT-C2-002.ipynb` notebook
**THEN** the 2D plot of the structure space should show the 'random' samples clustered in dense areas
**AND** the 'fps' samples should be visibly more spread out across the entire space than the 'random' samples.

---

**GIVEN** the user has configured a long-running job in the Web UI
**WHEN** they submit the job
**THEN** the status display on the page should update its text at least 3 times without the user refreshing the page
**AND** the sequence of status messages should follow the logical progression of the pipeline (e.g., "Generating" -> "Exploring" -> "Finished").
