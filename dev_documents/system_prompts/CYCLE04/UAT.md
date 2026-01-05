# Cycle 4 User Acceptance Testing (UAT): Robustness and Intelligence

This UAT plan for Cycle 4 is designed to allow a user to verify the major new features related to pipeline robustness, advanced material generation, and automated, physics-based decision-making. The tests will provide tangible proof that the application is now more resilient, capable, and "smarter" than in previous versions.

## 1. Test Scenarios

The tests will be conducted using a Jupyter Notebook (`UAT_Cycle4.ipynb`), enabling interactive execution and clear visualization of the results.

| Scenario ID | Scenario Name                               | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C4-01   | Verification of Simulation Stability with ZBL | High     |
| UAT-C4-02   | Verification of `IonicGenerator`            | High     |
| UAT-C4-03   | Verification of Auto Ensemble Switching     | High     |
| UAT-C4-04   | Verification of Fault Tolerance             | High     |

### UAT-C4-01: Verification of Simulation Stability with ZBL

**Description:**
This scenario provides a direct, dramatic demonstration of the benefit of ZBL potential mixing. The user will run a very high-temperature simulation, which is prone to crashing, with and without the ZBL repulsion enabled. This allows for a clear before-and-after comparison.

**Execution via Jupyter Notebook (`UAT_Cycle4.ipynb`):**
1.  **Setup (ZBL Off):** A cell will create a configuration file (`config_zbl_off.yaml`).
    *   `generator`: A simple alloy or ionic structure.
    *   `exploration`: Set a very high temperature (e.g., 2500 K) and set `enable_zbl_repulsion: false`.
2.  **Execution (ZBL Off):** The pipeline is run. The user should observe the log output, expecting to see messages indicating that one or more simulations failed and were caught.
3.  **Setup (ZBL On):** An identical configuration file (`config_zbl_on.yaml`) is created, but this time with `enable_zbl_repulsion: true`.
4.  **Execution (ZBL On):** The pipeline is run again.
5.  **Verification:** The user will compare the log outputs from both runs. The log from the ZBL-off run should show several "Caught exception in worker" messages. The log from the ZBL-on run should show significantly fewer, or ideally zero, such messages. The notebook will also check for the presence of files in the `quarantine/` directory for the first run, and its absence in the second. This provides clear proof that the ZBL feature dramatically improves simulation stability.

### UAT-C4-02: Verification of `IonicGenerator`

**Description:**
This scenario allows the user to test the new `IonicGenerator`. The goal is to verify that it can be configured to generate a valid ionic compound and that the resulting structure is charge-neutral.

**Execution via Jupyter Notebook (`UAT_Cycle4.ipynb`):**
1.  **Setup:** A cell will create a configuration file (`config_ionic.yaml`) that uses the new generator.
    *   `generator`: Set `mode: "ionic"` and provide parameters for a common ionic compound like Magnesium Chloride (MgCl2), including elements and their oxidation states (`Mg: 2`, `Cl: -1`).
2.  **Execution:** The pipeline's generation and storage stages are run.
3.  **Verification:** The user will connect to the output database. A code cell will retrieve one of the generated structures.
    *   **Composition Check:** The notebook will count the number of Mg and Cl atoms and assert that the ratio is 1:2, as expected for MgCl2.
    *   **Charge Neutrality Check:** The code will sum the oxidation states of all atoms in the structure. The user will verify that the total charge is exactly zero.
    *   **Visualization:** The structure will be visualized, allowing the user to see the generated ionic crystal.

### UAT-C4-03: Verification of Auto Ensemble Switching

**Description:**
This scenario demonstrates the pipeline's new "intelligence." The user will provide the system with two very different types of structures—a bulk crystal and a surface slab—and verify that the system automatically chooses the correct simulation ensemble (NPT for bulk, NVT for slab) for each.

**Execution via Jupyter Notebook (`UAT_Cycle4.ipynb`):**
1.  **Setup:** The user will not run the `generation` stage. Instead, two `xyz` files will be provided: `bulk.xyz` (a standard crystal) and `slab.xyz` (a crystal with a vacuum layer). The database will be pre-populated with these two structures.
2.  **Configuration:** A configuration file (`config_auto_ensemble.yaml`) will be created. The `exploration` section will be enabled with `ensemble: "auto"`.
3.  **Execution:** The pipeline's `exploration` stage is run. The user needs to be able to see which ensemble was chosen for each structure. This requires the `HybridMDMC` worker to be modified to log the ensemble it is using.
4.  **Verification:** The user will inspect the log output from the pipeline. They will verify two key lines in the log:
    *   "Running simulation for structure 1 using NPT ensemble."
    *   "Running simulation for structure 2 using NVT ensemble."
    This provides direct confirmation that the `detect_vacuum` function and the engine's logic are working together correctly.

### UAT-C4-04: Verification of Fault Tolerance

**Description:**
This scenario specifically tests the `ExplorationEngine`'s ability to handle a worker crash without halting the entire pipeline. While this is implicitly tested in UAT-C4-01, this test makes it the primary focus.

**Execution via Jupyter Notebook (`UAT_Cycle4.ipynb`):**
1.  **Setup:** The user will configure a pipeline to run exploration on 4 structures. To guarantee a failure, one of the input `xyz` files will be deliberately corrupted (e.g., by placing two atoms at the exact same coordinate), which will cause the MLIP calculator to fail.
2.  **Execution:** The pipeline is run.
3.  **Verification:** The user will inspect the output.
    *   **Graceful Handling:** The user will verify that the log shows an error being caught for the corrupted structure (e.g., "Worker for structure 3 failed. Saving to quarantine.").
    *   **Pipeline Completion:** Crucially, the user will verify that the simulations for the other 3 (valid) structures ran to completion and that the pipeline finished with a success message.
    *   **Quarantine Check:** The user will check that a `quarantine/` directory was created and that it contains the single corrupted structure file. This confirms the system's resilience.

## 2. Behavior Definitions

**Scenario: Verification of Simulation Stability with ZBL**

*   **GIVEN** a configuration for a very high-temperature MD simulation.
*   **WHEN** the user runs the exploration pipeline with `enable_zbl_repulsion: false`.
*   **THEN** the logs should indicate that one or more simulation workers failed and their errors were caught.
*   **AND** when the user re-runs the same pipeline with `enable_zbl_repulsion: true`.
*   **THEN** the logs should indicate that all simulation workers completed successfully.

**Scenario: Verification of `IonicGenerator`**

*   **GIVEN** a configuration to generate MgCl2, where Mg has an oxidation state of +2 and Cl has an oxidation state of -1.
*   **WHEN** the user runs the generation pipeline.
*   **THEN** the generated structures in the database must contain Mg and Cl atoms in a 1-to-2 ratio.
*   **AND** the sum of the oxidation states of all atoms in any generated structure must be 0.

**Scenario: Verification of Auto Ensemble Switching**

*   **GIVEN** a database containing one bulk crystal structure and one surface slab structure.
*   **AND** a configuration with `exploration.ensemble: "auto"`.
*   **WHEN** the user runs the exploration pipeline.
*   **THEN** the pipeline's log output must clearly state that it is running the bulk structure with the NPT ensemble.
*   **AND** the pipeline's log output must clearly state that it is running the slab structure with the NVT ensemble.

**Scenario: Verification of Fault Tolerance**

*   **GIVEN** a set of 4 input structures, where one is deliberately corrupted to cause a simulation failure.
*   **WHEN** the user runs the exploration pipeline.
*   **THEN** the pipeline must not crash and should run to completion.
*   **AND** the log output must show that 3 simulations succeeded and 1 simulation failed and was handled gracefully.
*   **AND** the database must be updated with the results of the 3 successful simulations.
*   **AND** the single corrupted structure file must be saved to a `quarantine/` directory.
