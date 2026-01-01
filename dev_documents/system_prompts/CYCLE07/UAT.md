# User Acceptance Testing (UAT): CYCLE07 - Advanced Active Learning & Boundary Treatment

## 1. Test Scenarios

This UAT is designed for an expert user to validate the sophisticated refinements made to the active learning loop. The focus is on ensuring the data being generated is of the highest possible physical quality and that the learning process itself is efficient and adaptive.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C07-01  | Verify correct boundary treatment for a trapped covalent structure| High     |
| UAT-C07-02  | Observe adaptation of the dynamic uncertainty threshold           | High     |
| UAT-C07-03  | Verify force masking during the training process                  | Medium   |

---

### **Scenario UAT-C07-01: Verify correct boundary treatment for a trapped covalent structure**

**(Min 300 words)**
This scenario tests the critical new feature of ensuring that structures sent for DFT labelling are physically sound. The user will intentionally trigger the trapping of a structure from a covalent material (like silicon) and then meticulously inspect the "cleaned" structure that is sent to the `LabelingEngine`.

The user will set up an active learning run for Silicon. They will use a model trained only on the perfect crystal structure and run a high-temperature simulation to quickly generate amorphous, uncertain configurations. When the simulation traps a structure, the user will pause the workflow before the DFT calculation is submitted. They will retrieve this trapped structure from the database or a temporary file.

The primary acceptance criterion is a visual and structural analysis of this atom cluster. The user will open the structure in a visualization tool. They should observe that it is a non-periodic cluster of atoms. At the outer edge of the cluster, they must be able to identify Hydrogen atoms that have been automatically added. They will verify that these hydrogens are bonded to Silicon atoms that would otherwise have been "dangling" (i.e., they were at the edge of the cutoff radius). The presence and correct placement of these passivating hydrogen atoms is the key evidence that the boundary treatment is working. The user will also check the `atoms.info` dictionary of the structure to confirm the presence of a `force_mask` array, which identifies the original "core" atoms.

---

### **Scenario UAT-C07-02: Observe adaptation of the dynamic uncertainty threshold**

**(Min 300 words)**
This scenario validates the new dynamic uncertainty threshold. The user wants to confirm that the pipeline becomes "smarter" as it learns, adjusting its sensitivity to uncertainty based on the maturity of the model. A good dynamic threshold should be low at the beginning (when the model knows little) and get higher as the model's overall confidence increases.

The user will set up and run a multi-generation active learning workflow (e.g., for 4-5 generations). The configuration will be set to use a dynamic threshold, e.g., `"dynamic_95percentile"`. The user's main task is to monitor the logs at the beginning of each generation.

The acceptance criteria are based on the logged values of the uncertainty threshold.
*   **Generation 1:** The user should see a log message like: "Calculating dynamic threshold from 50 training points... New threshold set to 1.85."
*   **Generation 2:** After a new structure has been added and the model retrained, the user should see a new calculation: "Calculating dynamic threshold from 51 training points... New threshold set to 2.15."
*   **Generation 3:** "Calculating dynamic threshold from 52 training points... New threshold set to 2.30."

The user must observe a clear, monotonically increasing trend in the calculated threshold across the generations. This trend is direct proof that the system is correctly assessing the model's growing confidence (the distribution of uncertainties on the training set is widening) and is adapting its criterion for what it considers a "novel" or "surprising" new structure.

---

### **Scenario UAT-C07-03: Verify force masking during the training process**

**(Min 300 words)**
This scenario tests the final step of the boundary treatment process: ensuring that the unphysical forces on the buffer atoms of a trapped cluster are ignored during model training. This is a subtle but critical point for model accuracy.

The user will perform the test from UAT-C07-01 to generate a trapped, passivated cluster of silicon with a `force_mask` applied. They will allow this structure to be labelled by the `LabelingEngine`. The key part of this UAT is to inspect the training process itself. The user will need to enable verbose logging for the `TrainingEngine`.

When the `TrainingEngine` starts its next run, its training set will include this new, masked structure. The user must find log messages that indicate the force mask is being used. For example, the `TrainingEngine` should log: "Structure ID 53 contains a force mask. Applying mask to loss function." During the training output, where the loss values are printed, the user should be able to confirm (if the library supports it) that the contribution to the force error from this structure is lower than it would be otherwise, or that the number of forces being included in the loss is only for the "core" atoms. A simpler, more practical criterion is to check that the final trained model has not been corrupted. The user could check, for example, that the forces on a perfect, bulk silicon crystal are still zero, which might not be the case if the model had been trained on large, unphysical forces from the buffer atoms of the trapped cluster.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C07-01 - Verify correct boundary treatment for a trapped covalent structure**

*   **GIVEN** an active learning simulation is running for a covalent material like Silicon.
*   **AND** the simulation encounters a high-uncertainty atomic configuration.
*   **AND** the `SimulationEngine` traps the frame containing this configuration.
*   **AND** the configuration specifies that boundary treatment with passivation is enabled.
*   **WHEN** the `SimulationEngine` processes this trapped frame before adding it to the database.
*   **THEN** a new, smaller, non-periodic `ase.Atoms` object should be created.
*   **AND** this new cluster should contain the high-uncertainty atom and its neighbors within a defined radius.
*   **AND** if any Si atoms in the "core" region had their bonds cut at the edge of the cluster, a Hydrogen atom should be automatically added to passivate that dangling bond.
*   **AND** the final `ase.Atoms` object that is added to the database for labelling must contain these passivating Hydrogen atoms.
*   **AND** the `info` dictionary of this `ase.Atoms` object must contain a boolean array called `force_mask` that distinguishes core atoms from buffer/passivation atoms.

---

**Scenario: UAT-C07-02 - Observe adaptation of the dynamic uncertainty threshold**

*   **GIVEN** an active learning workflow configured to run for multiple generations.
*   **AND** the `uncertainty_threshold` parameter is set to `"dynamic_95percentile"`.
*   **WHEN** the `WorkflowOrchestrator` starts the first generation.
*   **THEN** it should calculate a threshold, `T1`, based on the initial dataset's uncertainty distribution.
*   **AND** it should use `T1` for the first simulation run.
*   **AND** after this run adds a new data point and the model is retrained.
*   **WHEN** the `WorkflowOrchestrator` starts the second generation.
*   **THEN** it must re-calculate a new threshold, `T2`, based on the now-larger dataset.
*   **AND** `T2` should be greater than `T1`, reflecting the model's increased confidence.
*   **AND** the second simulation run must use the new, higher threshold `T2`.
*   **AND** this process of re-calculating an increasing threshold should repeat for all subsequent generations.
