# User Acceptance Testing (UAT): CYCLE05 - Advanced Sampling (DIRECT)

## 1. Test Scenarios

This UAT focuses on verifying the intelligent data selection capabilities of the new `Sampler` module. The user, acting as a research scientist, wants to ensure that the pipeline is not just generating a lot of data, but is wisely selecting only the most valuable data points for expensive DFT calculations. The success of this cycle is measured by the dramatic reduction in the number of candidate structures after the sampling step.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C05-01  | Verify reduction of structures after sampling a stable trajectory | High     |
| UAT-C05-02  | Verify sampler selects diverse structures from a mixed trajectory | High     |
| UAT-C05-03  | User can control the size of the sampled dataset                  | Medium   |

---

### **Scenario UAT-C05-01: Verify reduction of structures after sampling a stable trajectory**

**(Min 300 words)**
This scenario tests the sampler's ability to recognize and discard redundant data. A long MD simulation of a stable crystal at low temperature will produce thousands of frames, but most of them are very similar, representing small thermal vibrations. A brute-force approach would waste huge amounts of resources by calculating nearly identical structures. This UAT will confirm that the sampler correctly identifies this situation and selects only a few representative frames.

The user will set up a workflow to run a long MACE-driven MD simulation (e.g., 20,000 steps) for a stable crystal structure at a low temperature (e.g., Si at 300 K). They will configure the system to save all 20,000 frames from the `Explorer`'s trajectory. When the `Sampler` module runs, the user will monitor the logs. The primary acceptance criterion is a clear log message indicating the input and output sizes. The user expects to see a message like: "Sampler received 20,000 structures. After clustering and sampling, 50 structures were selected."

The user will then inspect the structures that were added to the database for labelling. They will verify that the number of new `'unlabeled'` entries is indeed the small number reported in the logs (e.g., 50), not the original 20,000. This demonstrates the massive data reduction that is the key benefit of this module. By visualizing a few of the selected structures, the user should see that they are all very similar, confirming that the sampler correctly identified that this trajectory explored only one main basin on the potential energy surface.

---

### **Scenario UAT-C05-02: Verify sampler selects diverse structures from a mixed trajectory**

**(Min 300 words)**
This scenario tests the sampler's intelligence in a more complex situation. The user will generate a trajectory that contains distinct and diverse atomic environments, such as a simulation of a melting crystal or a phase transformation. The goal is to verify that the sampler does not just reduce the data, but that it intelligently preserves examples of all the important structural variations present in the raw data.

To create a suitable test case, the user will run a MACE-driven MD simulation at a high temperature, designed to induce a phase change or create defects. For example, they could simulate a small Si crystal being heated until it becomes amorphous. The raw trajectory from the `Explorer` will thus contain frames of the perfect crystal, distorted intermediate structures, and fully amorphous configurations.

After the `Sampler` processes this trajectory, the user will inspect the small set of selected structures. The acceptance criteria are qualitative. The user will visualize all of the selected structures. They must be able to identify representatives from each stage of the simulation. They should find at least one structure that looks like the initial, stable crystal. They should find several that are highly distorted or contain defects. And they should find at least one that is clearly amorphous. If the selected set contains this structural diversity, it proves that the descriptor-based clustering is working. It successfully identified the different "types" of structures in the trajectory and the stratified sampling ensured that representatives from all types, even rare ones, were preserved for DFT labelling.

---

### **Scenario UAT-C05-03: User can control the size of the sampled dataset**

**(Min 300 words)**
This scenario tests the configurability of the sampler. An expert user might want to control the trade-off between the cost of the DFT calculations and the quality of the resulting potential. They can do this by changing the number of clusters the sampler uses, which directly influences the size of the final dataset. This UAT verifies that the system respects this user-defined parameter.

The user will run the same workflow twice on an identical trajectory. In the first run, they will set the number of clusters to a low value in their configuration file: `sampler: {num_clusters: 50}`. They will run the sampling step and record the number of structures that are selected for labelling (e.g., the system logs "35 structures were selected").

In the second run, they will change the configuration to a higher value: `sampler: {num_clusters: 200}`. They will re-run the sampling step on the exact same trajectory data. The acceptance criterion is that the number of selected structures in the second run is significantly larger than in the first run (e.g., the system now logs "150 structures were selected"). This confirms that the `num_clusters` parameter provides a direct and effective "knob" for the user to control the size and diversity of the training set, allowing them to tailor the workflow to their specific computational budget and accuracy requirements.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C05-01 & UAT-C05-02 - Combined Gherkin for Intelligent Sampling**

*   **GIVEN** the pipeline has run a long MACE-driven MD simulation of a Si crystal at high temperature.
*   **AND** the `Explorer` module has produced a trajectory containing 10,000 frames.
*   **AND** this trajectory contains frames corresponding to both crystalline and amorphous Si.
*   **AND** the `Sampler` is configured to select from this trajectory.
*   **WHEN** the `WorkflowOrchestrator` passes the 10,000 frames to the `Sampler.sample` method.
*   **THEN** the system should log that it is starting the sampling process on 10,000 input structures.
*   **AND** the sampler should calculate a descriptor (SOAP) for each of the 10,000 frames.
*   **AND** the sampler should perform a clustering analysis on the 10,000 descriptor vectors.
*   **AND** the system should log a message indicating the final number of selected structures, which must be much smaller than 10,000 (e.g., around 100-200).
*   **AND** only this small number of selected structures should be added to the database with the `state` set to `'unlabeled'`.
*   **AND** when a user visualizes the set of selected structures, they must find at least one structure that is clearly crystalline.
*   **AND** the user must also find at least one structure that is clearly amorphous.
*   **AND** the subsequent `LabelingEngine` step should only be triggered for this small, diverse subset of structures.

---

**Scenario: UAT-C05-03 - User can control the size of the sampled dataset**

*   **GIVEN** the pipeline has generated a fixed trajectory of 5,000 frames.
*   **AND** the user has configured the `input.yaml` with `sampler: {num_clusters: 20}`.
*   **WHEN** the `Sampler` module processes the trajectory.
*   **THEN** the system should select a certain number of structures, `N1`.
*
*   **GIVEN** the user has modified the `input.yaml` to `sampler: {num_clusters: 100}`.
*   **WHEN** the `Sampler` module processes the *exact same* 5,000-frame trajectory again.
*   **THEN** the system should select a new, larger number of structures, `N2`.
*   **AND** `N2` must be significantly greater than `N1`.
*   **AND** the user can thereby conclude that the `num_clusters` parameter directly and effectively controls the size of the final training dataset.
