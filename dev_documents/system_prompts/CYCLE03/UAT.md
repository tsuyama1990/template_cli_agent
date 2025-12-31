# CYCLE03/UAT.md

## 1. Test Scenarios

User Acceptance Testing (UAT) for Cycle 3 is designed to verify the effectiveness and efficiency of the new exploration and sampling workflow. The key user story for this cycle is: "As a computational materials scientist, I want the system to intelligently explore a vast range of possible atomic structures using a fast surrogate model, and then automatically select a small, information-rich subset for expensive DFT calculations, ensuring I get the most accurate potential for my budget." The UAT will be conducted via a Jupyter Notebook (`uat_cycle03.ipynb`) that allows the user to visualize the results of the exploration and sampling process, building confidence in this critical optimization step.

| Scenario ID | Scenario Description                                                                | Priority |
| :---------- | :---------------------------------------------------------------------------------- | :------- |
| UAT-C03-01  | **Verify Surrogate Model MD Exploration**                                           | High     |
|             | The user will run the full `Generate -> Explore -> Label -> Train` workflow. They will verify that the system successfully uses the MACE surrogate model to run a high-temperature Molecular Dynamics (MD) simulation, which serves as the exploration phase. The goal is to confirm that this computationally intensive step executes correctly. |          |
| UAT-C03-02  | **Validate Intelligent Sampling and Data Reduction**                                | High     |
|             | The user will inspect the state of the ASE database before and after the exploration module runs. They will verify that the number of structures selected for DFT labeling is exactly the number requested in the configuration and is significantly smaller than the total number of frames in the MD trajectory. This demonstrates the data reduction and intelligent filtering capability of the DIRECT sampler. |          |
| UAT-C03-03  | **Visualize Sampled Structure Diversity**                                           | Medium   |
|             | The user will use the Jupyter notebook to load the structures selected by the sampler. The notebook will generate plots (e.g., a histogram of potential energies, a 2D projection of the structural descriptors using PCA/t-SNE) to visually confirm that the selected structures are diverse and cover a wide range of configurations, from low-energy stable states to high-energy distorted states. This provides qualitative evidence that the sampling is not just random, but intelligent. |          |

## 2. Behavior Definitions

The behavior of the system will be defined using the Gherkin-style (GIVEN/WHEN/THEN) syntax. These definitions will be implemented as interactive steps within the `uat_cycle03.ipynb` notebook.

### Scenario: UAT-C03-01 - Verify Surrogate Model MD Exploration

*   **GIVEN** a minimal `input.yaml` file for a simple system like Aluminum ("Al").
*   **AND** a configuration file that enables the exploration stage and specifies a short MD run (e.g., 1000 steps).
*   **WHEN** the user executes the main CLI command `uv run mlip-pipe run --input input.yaml`.
*   **THEN** the command should execute and finish without errors.
*   **AND** the system's log output should clearly show messages indicating that it is "Starting surrogate model MD exploration," running a specific number of steps, and then "Finished exploration."
*   **AND** the log should also indicate that it is "Calculating descriptors for trajectory" and subsequently "Performing DIRECT sampling."

### Scenario: UAT-C03-02 - Validate Intelligent Sampling and Data Reduction

*   **GIVEN** the successful completion of the workflow from scenario UAT-C03-01.
*   **AND** the configuration file specified that the sampler should select a specific number of structures (e.g., `num_samples_to_select: 50`).
*   **AND** the MD simulation was configured to run for 1000 steps, generating 1000 candidate structures.
*   **WHEN** the user connects to the output ASE database (`mlip_autopipec.db`) within the Jupyter Notebook.
*   **AND** the user queries the database for all structures that have been marked as `status='labeled'`.
*   **THEN** the total number of labeled structures in the database should be exactly 50.
*   **AND** this number (50) should be significantly less than the total number of MD steps (1000), demonstrating a data reduction ratio of 20:1.

### Scenario: UAT-C03-03 - Visualize Sampled Structure Diversity

*   **GIVEN** the successful completion of the workflow from scenario UAT-C03-01.
*   **AND** the user has loaded the 50 sampled `Atoms` objects from the database into the Jupyter Notebook.
*   **AND** the user has also loaded the full MD trajectory file produced by the exploration step.
*   **WHEN** the user generates a histogram of the potential energies of all 1000 structures from the full MD trajectory.
*   **AND** overlays a scatter plot showing the potential energies of the 50 selected samples on the same histogram.
*   **THEN** the scatter points should be distributed across the entire range of the histogram, not just clustered at the bottom. This visually confirms that the sampler is selecting not only stable, low-energy structures but also high-energy, distorted structures that are crucial for learning repulsive interactions.
*   **AND** when the user generates a 2D t-SNE plot of the SOAP descriptors for all 1000 trajectory frames.
*   **AND** overlays the points corresponding to the 50 selected samples.
*   **THEN** the points for the selected samples should be spread out across the entire t-SNE map, indicating they are structurally diverse and not redundant.
