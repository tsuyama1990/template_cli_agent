# CYCLE 02: USER ACCEPTANCE TESTING (UAT)

This User Acceptance Testing (UAT) plan for Cycle 02 focuses on verifying the advanced, intelligent, and autonomous capabilities of the MLIP-AutoPipe system. The central goal is to provide a clear and compelling demonstration that the active learning loop, the core feature of this cycle, is not only functional but also effective. The UAT will test the system's ability to perform its most critical task: to autonomously identify the weaknesses in its own model, intelligently gather new data to address those weaknesses, and, as a result, systematically and measurably improve its own predictive accuracy over several generations. This test is the ultimate validation of the project's core promise.

Similar to the previous cycle, this UAT will be facilitated through an interactive Jupyter Notebook, `C02_UAT_Notebook.ipynb`. This format is chosen to make the complex, multi-stage process of active learning transparent and understandable to the user. The notebook will provide a guided walkthrough of a more scientifically interesting test case than the one used in Cycle 01. It will not only allow the user to execute the full active learning pipeline but also provide them with the analytical tools to dissect the results. The user will be able to see *why* the model chose to learn certain structures and to quantify the improvement in the model's accuracy, providing a tangible and satisfying user experience. This also serves as an advanced tutorial for how a user might analyze the results of their own future research runs.

## 1. Test Scenarios

### Scenario ID: UAT-C02-001
**Priority:** High
**Title:** Autonomous Potential Refinement via Active Learning for a Binary Alloy (FePt)

**Description:**
This UAT scenario is designed to be a comprehensive validation of the on-the-fly (OTF) active learning loop, which is the main deliverable of Cycle 02. To provide a meaningful test, the chosen material is Iron-Platinum (FePt). This system is significantly more complex than the Silicon used in Cycle 01; it is a binary alloy with a rich phase diagram, complex magnetic properties, and important technological applications, making it an excellent showcase for the power of the pipeline. The goal of this UAT is not to produce a state-of-the-art, publication-quality potential, but to provide undeniable proof that the system can take a rudimentary initial model and autonomously improve it through iterative, uncertainty-driven learning.

The entire user experience for this test will be managed within the `C02_UAT_Notebook.ipynb`. This notebook will guide the user through the setup, execution, and, most importantly, the detailed analysis of an active learning run.
1.  **Advanced Configuration**: The notebook will start by instructing the user on how to create the `input.yaml` file for FePt. This time, the file will not only contain the system definition but also a new `active_learning` section. The user will be guided to enable the active learning process and set key parameters, such as `max_generations: 3`, to control the duration of the run. This step demonstrates the new level of user control over the workflow.
2.  **Pipeline Execution**: The user will execute a cell to launch the full active learning pipeline. Because this process involves multiple generations of DFT calculations and model training, it will be a longer process than the test in Cycle 01. The notebook will be configured to stream the live log output from the pipeline directly into the cell's output. This provides a transparent view of the system's operations, allowing the user to see messages indicating the start and end of each generation, the detection of high-uncertainty structures, and the triggering of retraining events.
3.  **Introspective Verification**: This is the core of the UAT. After the pipeline has completed its run, the notebook will provide a series of analytical tools to dissect the results and verify that the active learning process worked as intended.
    *   **Generational Data Growth**: The user will execute a code cell that queries the database and plots a bar chart showing the number of new structures that were added in each generation. The expected result is to see a larger number of structures in "Generation 0" (the initial dataset) and a smaller, more focused number of structures in the subsequent "Generation 1" and "Generation 2". This visually demonstrates that the system is selectively adding data, not just generating it randomly.
    *   **Provenance of Learned Structures**: The notebook will guide the user to query the database for the specific structures that were added during the active learning generations. It will then provide code to visualize one of these structures and its associated uncertainty score, showing tangible proof of the uncertainty-driven selection process.
    *   **Quantifiable Accuracy Improvement**: This is the ultimate validation step. The notebook will contain a small, independent hold-out test set of FePt configurations that were not included in any stage of the training process. This test set will be designed to probe aspects of the potential energy surface that the initial model would likely struggle with. The user will execute a code cell that loads the initial "Generation 0" model and the final, refined model. It will then calculate and print the Root Mean Squared Error (RMSE) on the force predictions for both models against the hold-out set. A successful pass for this UAT requires the RMSE of the final model to be significantly and demonstrably lower than the RMSE of the initial model.

Passing this comprehensive UAT scenario provides clear, quantitative evidence that the system is not just an automation tool, but a true intelligent agent capable of autonomous learning and self-improvement.

## 2. Behavior Definitions

The following Gherkin-style behavior definitions describe the key observable outcomes of the new features introduced in Cycle 02. These definitions form the basis for the assertions and verification steps that will be carried out within the `C02_UAT_Notebook.ipynb`.

---

**Feature: Surrogate-based Exploration for Initial Seeding**
As a user, I want the system to leverage a fast, general-purpose surrogate model to intelligently explore the material's configuration space, so that the initial dataset is more diverse and informative than one based on simple lattice perturbations.

**Scenario:** The pipeline uses the MACE surrogate model and DIRECT sampling to generate an initial dataset for FePt.
**GIVEN** the user's `input.yaml` file specifies the active learning strategy as `direct_then_active`.
**WHEN** the user executes the full active learning pipeline command.
**THEN** the initial log messages displayed in the notebook must contain strings indicating that the `ExplorerSampler` is active, such as "Running surrogate MD with MACE model" or "Performing DIRECT sampling".
**AND** after the "Generation 0" training is complete, a query to the database must show that the initial dataset contains structures with varied stoichiometries or lattice types, beyond simple deformations of a single input structure.

---

**Feature: On-the-Fly Uncertainty Detection during Simulation**
As a user, I want the system to be aware of its own limitations and to automatically detect when it is encountering atomic configurations that it cannot reliably predict, so that it knows when and where it needs to learn more.

**Scenario:** An active learning MD simulation, using a partially trained MLIP, encounters a novel atomic configuration that is outside its training distribution.
**GIVEN** the system is executing an MD simulation within an active learning generation (e.g., using the "Generation 1" model).
**AND** the simulation evolves to a state containing an atomic environment (e.g., a defect structure) that is significantly different from any structure in the current training set.
**WHEN** the `SimulationEngine` evaluates the forces for this novel structure.
**THEN** the system's log output, streamed to the notebook, must contain a specific, clear message indicating the event, such as "Uncertainty of [some value] exceeds dynamic threshold of [some value]. Yielding structure for DFT calculation."
**AND** the main MD simulation for that generation must be paused.
**AND** a new entry for this novel structure must be added to the database with the status "uncalculated", marking it for future labeling.

---

**Feature: Autonomous Retraining and Generational Model Improvement**
As a user, I want the system to automatically take the new, high-uncertainty structures it has found, calculate their true properties with DFT, and use this new knowledge to retrain and improve its own model, creating a virtuous cycle of refinement.

**Scenario:** The system successfully completes a full active learning generation cycle after identifying one or more uncertain structures.
**GIVEN** the `SimulationEngine` has yielded at least one new high-uncertainty structure to the `WorkflowOrchestrator`.
**WHEN** the `WorkflowOrchestrator` processes this new structure and completes the generation.
**THEN** the system's logs must show that the `LabelingEngine` was invoked for the new structure.
**AND** the logs must subsequently show a message indicating that the `TrainingEngine` is starting a new training run, such as "Training Generation 2 model with [N+1] total structures...".
**AND** a new, distinct model file, such as `model_gen_2.ace`, must be created in the working directory, and it must be different from the previous generation's model file.

---

**Feature: Demonstrable and Quantifiable Accuracy Improvement**
As a user, I need to be confident that the computationally expensive active learning process is providing a real benefit. I want to see quantitative proof that the final model is more accurate than the initial one.

**Scenario:** A user compares the predictive accuracy of the first and the last models produced by a multi-generation active learning pipeline.
**GIVEN** a completed active learning run has produced a sequence of models, including `model_gen_0.ace` and `model_gen_3.ace`.
**AND** a separate, hold-out set of validation structures, which the models have never seen, is available in a file named `validation_set.xyz`.
**WHEN** the user executes the validation script in the UAT notebook, which calculates the Root Mean Squared Error (RMSE) of the force predictions for both models against the validation set.
**THEN** the printed output in the notebook must show the RMSE for both models.
**AND** the calculated RMSE for `model_gen_3.ace` must be numerically smaller than the RMSE for `model_gen_0.ace`, providing concrete evidence of the model's improvement.
