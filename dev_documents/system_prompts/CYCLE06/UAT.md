# User Acceptance Testing (UAT): CYCLE06 - On-the-Fly Simulation & Active Learning

## 1. Test Scenarios

This UAT is designed to validate the core "self-improving" feature of the MLIP-AutoPipe. The user, a materials scientist, will test the system's ability to run a simulation with a trained MLIP, detect when the model is uncertain, and automatically capture the uncertain structure for re-training. The focus is on the successful execution and verification of this feedback loop.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C06-01  | Successfully trap and label a high-uncertainty structure          | High     |
| UAT-C06-02  | Complete a multi-generation active learning run                   | High     |
| UAT-C06-03  | Verify simulation convergence when no uncertain structures are found | Medium   |

---

### **Scenario UAT-C06-01: Successfully trap and label a high-uncertainty structure**

**(Min 300 words)**
This scenario is the fundamental test of the active learning mechanism. The user wants to see the system in action: identify a weakness in its own model and capture the relevant data point. To achieve this, the user will start with a deliberately "imperfect" MLIP, trained on a very small and uniform dataset (e.g., only low-temperature crystal structures). They will then use this model to run a high-temperature simulation where new, unseen types of structures are likely to appear.

The user will start the `cdd active-learn` workflow. The system will first train the initial, imperfect model. Then, it will launch the `SimulationEngine`. The user will monitor the logs closely. The key acceptance criterion is to see the simulation stop prematurely. The system should log a message like, "Uncertainty of 6.7 exceeds threshold of 5.0 at step 1234. Trapping structure." This message is the primary evidence that the detection mechanism is working.

Following this, the user will see the `WorkflowOrchestrator` automatically send this newly trapped structure to the `LabelingEngine`. They should see logs from the Quantum Espresso calculation for this specific structure. Finally, the user will inspect the database. They will find a new entry corresponding to the trapped structure (frame 1234 from the simulation). They will verify that this new entry's `state` is now `'labeled'` and it has a valid `dft_result`. This successful test provides a complete, traceable record of the system identifying a model weakness and autonomously acquiring the necessary DFT data to correct it.

---

### **Scenario UAT-C06-02: Complete a multi-generation active learning run**

**(Min 300 words)**
This scenario tests the full, closed-loop nature of the active learning process over multiple iterations. The user wants to see the system not just find one uncertain structure, but to iteratively improve its model. The goal is to observe the system's behavior over several generations of training, simulating, trapping, and re-training.

The user will configure the pipeline for a multi-generation run, for example, setting `max_generations: 3`. They will start with a small initial dataset, similar to the previous scenario. When they launch the `cdd active-learn` command, they will monitor the overall process. The acceptance criteria are based on observing a clear, repeating pattern in the logs.

*   **Generation 1:** The user should see the system train an initial model (`Model_Gen1`), run a simulation, trap a new structure (`Structure_A`), and label it.
*   **Generation 2:** The user must see the system start a *new* training run. The logs should indicate that the training set now includes the initial data *plus* `Structure_A`. A new model, `Model_Gen2`, is created. The system then starts a new simulation using this improved model. It might run for longer before trapping another structure (`Structure_B`).
*   **Generation 3:** The user must again see a new training run, this time with a dataset containing the initial data, `Structure_A`, and `Structure_B`. `Model_Gen3` is created and a new simulation begins.

The successful observation of this distinct train->simulate->trap->label->retrain cycle over multiple generations is the core of this UAT. It proves that the orchestrator is correctly managing the state and data flow of the loop, leading to progressively more robust models.

---

### **Scenario UAT-C06-03: Verify simulation convergence when no uncertain structures are found**

**(Min 300 words)**
This scenario tests the "exit condition" of the active learning loop. A robust model should eventually be able to complete a full simulation without encountering any configurations that it deems uncertain. The user wants to verify that the pipeline correctly detects this state of "convergence" and terminates gracefully.

To test this, the user will start with a very good, pre-trained MLIP (perhaps one that has already gone through several active learning cycles). They will start the `cdd active-learn` workflow with this model. The system will launch the `SimulationEngine` to perform the OTF-MD simulation. The user will monitor the logs.

The key acceptance criterion is to see the simulation run to its full configured length (e.g., all 100,000 steps) without the uncertainty ever breaching the threshold. The system should print a log message at the end, such as "Simulation completed for 100,000 steps without finding any uncertain structures. The model appears to be converged for these conditions." Following this, the `WorkflowOrchestrator` should recognize that the simulation returned no new structure. It should then log a message like "Convergence reached. Terminating active learning loop." and the program should exit cleanly, without proceeding to another generation. This test is crucial as it demonstrates that the system has a stable stopping point and will not run indefinitely.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C06-01 & UAT-C06-02 - A Multi-Generation Active Learning Run**

*   **GIVEN** a database containing a small, initial set of labelled data points.
*   **AND** the active learning workflow is configured to run for a maximum of 3 generations.
*   **AND** the uncertainty threshold is set to a reasonable value.
*   **WHEN** the user starts the active learning workflow.
*   **THEN** **(Generation 1)** the `TrainingEngine` should train a model (`Model_Gen1`) using only the initial data.
*   **AND** the `SimulationEngine` should start an MD simulation using `Model_Gen1`.
*   **AND** the simulation should stop prematurely, logging that the uncertainty threshold was exceeded and a new structure (`Structure_A`) was trapped.
*   **AND** the `LabelingEngine` should be called to perform a DFT calculation on `Structure_A` and update it in the database.
*
*   **THEN** **(Generation 2)** the `TrainingEngine` should start a new training job.
*   **AND** the training data for this job must consist of the initial data plus the newly-labeled `Structure_A`.
*   **AND** a new model, `Model_Gen2`, should be created.
*   **AND** the `SimulationEngine` should start a new MD simulation using the improved `Model_Gen2`.
*   **AND** this simulation should run for more steps than the first one before potentially trapping another structure (`Structure_B`).
*   **AND** the `LabelingEngine` should be called for `Structure_B`.
*
*   **THEN** **(Generation 3)** the workflow should repeat, training `Model_Gen3` on all data collected so far.
*   **AND** the process should demonstrate a clear, iterative refinement loop.

---

**Scenario: UAT-C06-03 - Verifying Simulation Convergence**

*   **GIVEN** a database with a comprehensive set of labelled data, resulting in a robust MLIP.
*   **AND** the active learning workflow is configured to run.
*   **WHEN** the `TrainingEngine` trains the robust model.
*   **AND** the `SimulationEngine` starts an MD simulation with this model.
*   **THEN** the MD simulation should run for the entire number of configured steps (e.g., 100,000).
*   **AND** the system logs should show that the uncertainty value never exceeded the threshold during the run.
*   **AND** the `SimulationEngine` should terminate gracefully and return no new structure to the orchestrator.
*   **AND** the `WorkflowOrchestrator` should log a message indicating that convergence has been reached.
*   **AND** the active learning loop should terminate, and the program should exit, without starting a new generation.
