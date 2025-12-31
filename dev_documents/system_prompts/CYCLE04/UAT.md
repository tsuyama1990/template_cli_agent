# UAT.md: Cycle 04 - The Active Learning Loop

## 1. Test Scenarios

User Acceptance Testing for Cycle 04 is focused on verifying the fulfillment of the core promise of the MLIP-AutoPipe system: its ability to learn and improve itself autonomously through active learning. The user needs to gain a high degree of confidence that the on-the-fly (OTF) active learning loop is not just a theoretical concept but a practical, robust, and functioning mechanism. This UAT is designed to confirm that the system can reliably identify weaknesses in the MLIP during a simulation, and then systematically take the correct actions to remedy those weaknesses. The tests are designed to provide tangible, observable proof that the system is becoming "smarter," more robust, and more accurate over time, without any human intervention.

| Scenario ID | Scenario Description                                       | Priority |
| :---------- | :--------------------------------------------------------- | :------- |
| UAT-C04-01  | **End-to-End Verification of a Single, Complete Active Learning Cycle**         | High     |
|             | This is the most critical and comprehensive User Acceptance Test for this cycle. A user will start a simulation with a deliberately under-trained MLIP that has a known, significant "blind spot" (e.g., it has only been trained on the solid phase of a material and knows nothing about the liquid phase). The user's expectation is to observe the entire feedback loop in action: 1) the system should run the simulation successfully in the known domain, 2) it must detect high uncertainty when the simulation enters the unknown domain (the blind spot), 3) it must automatically pause the simulation and trigger a new DFT calculation on the novel structure, 4) it must then automatically retrain the MLIP with this new data, and finally, 5) it must be able to resume the simulation with the newly improved and more robust potential. This test provides a complete, end-to-end demonstration of the core value proposition of the entire project. |          |
| UAT-C04-02  | **Verification of the Dynamic and Adaptive Uncertainty Threshold**                           | Medium   |
|             | This scenario is designed to build trust in the "intelligence" of the learning mechanism. A user will inspect the system's log files across multiple, sequential active learning iterations. The user's expectation is to see clear evidence that the "uncertainty threshold" value is being automatically recalculated and adjusted at the start of each new simulation run. Specifically, the threshold should generally trend upwards as more data is added to the training set. This reflects the model's growing confidence and prevents the system from becoming overly sensitive and retraining on minor structural variations. This test confirms that the adaptive part of the learning mechanism is working as designed and that the system's behavior is stabilising as the model matures. |          |
| UAT-C04-03  | **Validation of the Correct and Physically Sound Extraction of Periodic Structures for Labelling**| High     |
|             | This scenario tests a critical, and technically challenging, aspect of the active learning loop. A user will run an OTF simulation of a bulk, periodic crystal. When an uncertain structure is detected and extracted, the user will manually inspect the specific atomic structure that is saved to the database for the subsequent DFT labelling. The user's expectation is that this extracted structure is a non-periodic, charge-neutral cluster that correctly includes not just the atom of highest uncertainty but also its complete surrounding environment (the buffer region). The user will verify that the structure is physically sound and free of common artefacts like unphysical dangling bonds. This UAT is essential for validating the critical, non-trivial structure extraction logic, which is a prerequisite for obtaining meaningful DFT labels. |          |

## 2. Behavior Definitions

These behaviors will be tested by a user running the active learning workflow from the command line, actively monitoring the log files for specific messages, and inspecting the state of the database and the filesystem at various key stages of the process.

### Scenario: UAT-C04-01 - End-to-End Verification of a Single, Complete Active Learning Cycle

*   **GIVEN** an MLIP has been previously trained on a limited dataset that is known to be incomplete, for example, a dataset containing only the solid, crystalline phase of a material like Silicon.
*   **AND** the user has configured and starts an on-the-fly (OTF) Molecular Dynamics simulation with this specific MLIP, setting the simulation temperature to a value high enough to induce melting.
*   **WHEN** the simulation runs, and the atomic system, driven by the high temperature, starts to lose its crystalline order and explore liquid-like, disordered configurations that were explicitly not part of the original training set.
*   **THEN** the system's console log must show a clear, explicit, and high-priority message indicating the detection of an extrapolation event, such as: "High uncertainty detected (Value: X exceeds Threshold: Y). Pausing simulation for retraining."
*   **AND** the log must subsequently show messages indicating that the `LabellingEngine` is being invoked to perform a new, expensive DFT calculation on the newly discovered structure.
*   **AND** following the successful completion of the DFT calculation, the log must show further messages indicating that the `TrainingEngine` is being invoked to retrain or update the MLIP model with the newly acquired data point.
*   **AND** the user must be able to verify on the filesystem that a new, updated MLIP model file has been saved, with a timestamp corresponding to the retraining event.
*   **AND** upon completion of the retraining, the system should be capable of automatically and seamlessly resuming the MD simulation from the exact point where it was paused, now using the newly improved and more robust potential.

### Scenario: UAT-C04-02 - Verification of the Dynamic and Adaptive Uncertainty Threshold

*   **GIVEN** a system has been configured to run for a minimum of three consecutive active learning iterations.
*   **WHEN** the user starts the main `run-active-learning-loop` command from the terminal.
*   **THEN** at the very beginning of the first simulation run, the console log must print a clear message stating the initial uncertainty threshold being used, for example: "Dynamic uncertainty threshold for iteration 1 set to: 0.5 eV/A".
*   **AND** after the first retraining cycle has successfully completed and the second simulation run is about to begin, the log must again print a *new* uncertainty threshold value, for example: "Dynamic uncertainty threshold for iteration 2 set to: 0.6 eV/A".
*   **AND** the user should be able to observe by inspecting the logs for all the iterations that this threshold value is generally trending upwards. While it may fluctuate slightly, the overall trend should be positive, reflecting the model's expanding training set and its consequently increasing baseline confidence. This provides clear evidence of the adaptive learning mechanism at work.

### Scenario: UAT-C04-03 - Validation of the Correct and Physically Sound Extraction of Periodic Structures for Labelling

*   **GIVEN** an on-the-fly (OTF) simulation of a bulk, periodic crystal (e.g., a supercell of FCC Aluminum) is running.
*   **AND** the system's log indicates that a high-uncertainty event has been detected and a structure has been extracted for relabelling.
*   **WHEN** the user connects to the `mlip.db` database and uses a script or client to extract the specific `ase.Atoms` object that was just added with the state 'active_learning_candidate'.
*   **THEN** upon opening this extracted structure in a molecular visualisation tool (like `ase gui`), the user must rigorously verify the following critical properties:
    *   The structure's periodic boundary conditions flag (`pbc`) must be set to `[False, False, False]`. The extracted structure must be an isolated, non-periodic cluster.
    *   The total number of atoms in the extracted structure must be significantly larger than the number of atoms in the simulation's unit cell. This provides evidence that a buffer region has been correctly included around a central core region.
    *   The visual representation of the structure must look physically sensible. There should be no atoms unnaturally sliced by a cell boundary or other unphysical artefacts that would result from a naive extraction method.
    *   If the material being simulated were an ionic crystal, the user would also verify that the total charge of the extracted cluster sums to zero, ensuring charge neutrality.
