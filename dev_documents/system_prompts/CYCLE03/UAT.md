# Cycle 3 User Acceptance Testing (UAT): Advanced Sampling and Hybrid Exploration

This document provides the User Acceptance Testing (UAT) plan for Cycle 3. This cycle introduced two major scientific enhancements: the ability to perform hybrid Molecular Dynamics / Monte Carlo (MD/MC) simulations and the implementation of advanced data sampling techniques (Random and Farthest Point Sampling). This UAT will allow a user to verify these new capabilities and appreciate the increased intelligence of the pipeline.

## 1. Test Scenarios

The UAT will again be presented as a Jupyter Notebook (`UAT_Cycle3.ipynb`) to provide an interactive and explanatory testing experience. The user will be able to clearly see the effect of the new features on the generated dataset.

| Scenario ID | Scenario Name                               | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C3-01   | Verification of Hybrid MD/MC Exploration    | High     |
| UAT-C3-02   | Verification of the Sampling Engine         | High     |
| UAT-C3-03   | Comparing FPS and Random Sampling           | High     |

### UAT-C3-01: Verification of Hybrid MD/MC Exploration

**Description:**
This scenario allows the user to verify the effect of the new hybrid MD/MC exploration mode. The goal is to confirm that enabling the Monte Carlo "atom swap" move leads to a more rapid exploration of the chemical space in an alloy, compared to a pure MD simulation.

**Execution via Jupyter Notebook (`UAT_Cycle3.ipynb`):**
1.  **Setup (MD only):** The first set of cells will create a configuration file (`config_md_only.yaml`) for a binary alloy (e.g., Ni-Al). The exploration will be configured to run a standard MD simulation at a moderate temperature for a set number of steps. The pipeline will be run, generating a trajectory.
2.  **Setup (MD/MC):** The next set of cells will create a second configuration file (`config_md_mc.yaml`) that is *identical* to the first, except it will enable the `mc_swap` feature with a reasonably high frequency. The pipeline will be run again with this new configuration, generating a second trajectory from the *same* starting structure.
3.  **Verification:** The notebook will then load both trajectories. The key verification is to show that the chemical ordering of the atoms has changed more significantly in the MD/MC run. A simple metric, like counting the number of Ni-Al bonds at the start and end of each simulation, will be calculated. The user will verify that the change in this bond count is substantially larger for the MD/MC trajectory, providing clear evidence that the atom swaps are successfully altering the material's chemical structure. A visualization of the final structures from both runs will also highlight the difference.

### UAT-C3-02: Verification of the Sampling Engine

**Description:**
This scenario is a straightforward test of the new `SamplingEngine`. The user will run the full pipeline and verify that the final number of structures stored in the database matches the number requested in the sampling configuration. This confirms that the engine is correctly selecting a subset of the candidates generated during exploration.

**Execution via Jupyter Notebook (`UAT_Cycle3.ipynb`):**
1.  **Setup:** A cell will create a configuration file (`config_sampling.yaml`).
    *   `exploration`: Configure to run for enough steps to generate a large number of candidate structures (e.g., >500).
    *   `sampling`: Enable this stage, set `mode: "random"`, and configure `num_samples: 50`.
    *   `db_path`: `sampling_test.db`.
2.  **Execution:** The pipeline will be run with this configuration.
3.  **Verification:** The user will connect to the output database `sampling_test.db`. A code cell will perform a query to count the number of structures in the final, curated dataset. The user will verify that this count is exactly 50, matching the `num_samples` parameter. This confirms the engine is correctly truncating the dataset to the desired size.

### UAT-C3-03: Comparing FPS and Random Sampling

**Description:**
This is the most important UAT for this cycle. It is designed to provide the user with a tangible and intuitive understanding of the benefit of Farthest Point Sampling (FPS) over a naive random selection. The user will generate two datasets—one with FPS, one with random sampling—and visually and quantitatively compare their diversity.

**Execution via Jupyter Notebook (`UAT_Cycle3.ipynb`):**
1.  **Setup:** The user will first run the pipeline (Generate -> Explore) to create a large pool of candidate structures from a single, long trajectory. This ensures both samplers start with the exact same set of candidates.
2.  **Run FPS:** A configuration file will be used with `mode: "fps"` and `num_samples: 20`. The pipeline's sampling and storage stages will be run, saving the results to `fps_dataset.db`.
3.  **Run Random:** A second configuration file will be used with `mode: "random"` and `num_samples: 20`. The sampling and storage stages will be run again on the *same* initial pool of candidates, saving the results to `random_dataset.db`.
4.  **Verification:** The notebook will provide a powerful comparison:
    *   **Quantitative Analysis:** A code cell will load both datasets. For each dataset, it will calculate the potential energy of every structure. It will then plot two histograms of the energy distributions on the same axes. The user will be able to see that the FPS-selected energies are typically much more spread out, covering a wider range of the potential energy surface, while the random samples are more clustered around the average energy. The standard deviation of the energy for the FPS set should be significantly higher.
    *   **Visual Analysis:** The notebook will display a grid of 5-6 structures from each dataset. The user will be able to visually inspect them. The structures from the random set will likely look very similar to each other, while the FPS set will contain more visually distinct and diverse configurations (e.g., some more ordered, some more disordered, some with defects). This provides compelling, intuitive proof of the superiority of FPS.

## 2. Behavior Definitions

**Scenario: Verification of Hybrid MD/MC Exploration**

*   **GIVEN** an initial binary alloy structure.
*   **WHEN** the user runs two separate exploration simulations from this structure: one with pure MD, and one with hybrid MD/MC.
*   **THEN** the trajectory from the hybrid MD/MC simulation should show evidence of atom swaps.
*   **AND** a quantitative measure of chemical ordering (like the number of bonds between different atom types) should show a significantly larger change over the course of the hybrid MD/MC simulation compared to the pure MD simulation.

**Scenario: Verification of the Sampling Engine**

*   **GIVEN** a configuration file that specifies `sampling.num_samples: 50`.
*   **AND** the exploration stage generates more than 50 candidate structures.
*   **WHEN** the user runs the full `generate -> explore -> sample -> store` pipeline.
*   **THEN** the final output database should be created successfully.
*   **AND** the database should contain exactly 50 structures in its main table.

**Scenario: Comparing FPS and Random Sampling**

*   **GIVEN** a large, fixed pool of candidate structures.
*   **WHEN** the user creates two final datasets from this pool, one using `sampling.mode: "fps"` and the other using `sampling.mode: "random"`, both selecting 20 structures.
*   **THEN** both `fps_dataset.db` and `random_dataset.db` should be created, each containing 20 structures.
*   **AND** the distribution of potential energies of the structures in `fps_dataset.db` should be visibly wider and have a higher standard deviation than the distribution of energies in `random_dataset.db`.
*   **AND** a visual inspection of the structures should show greater structural diversity in the FPS-selected set.
