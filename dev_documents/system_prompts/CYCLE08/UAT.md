# User Acceptance Testing (UAT): CYCLE08 - Advanced Simulation (kMC) & UI/UX

## 1. Test Scenarios

This final UAT is focused on two distinct areas: the new scientific capability of Kinetic Monte Carlo (kMC) and the overall polish and usability of the application. One user persona is the research scientist, who will validate the kMC engine. The other persona is a new user, who will evaluate the clarity and helpfulness of the command-line interface.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C08-01  | Successfully run a kMC simulation of vacancy diffusion            | High     |
| UAT-C08-02  | Experience clear progress indicators during a long workflow       | High     |
| UAT-C08-03  | Receive helpful and well-formatted output and error messages    | Medium   |

---

### **Scenario UAT-C08-01: Successfully run a kMC simulation of vacancy diffusion**

**(Min 300 words)**
This scenario tests the new kMC simulation engine. The user, a materials scientist, wants to study a classic rare-event problem: the diffusion of a single vacancy in a crystal lattice over a long timescale. This is something that is impossible to simulate with standard MD.

The user will start with a trained MLIP for a simple crystal (e.g., Aluminum) that is known to be accurate. They will create a large supercell of this crystal and manually remove one atom to create a vacancy. They will then create a configuration file for a kMC simulation, specifying the temperature and the number of events to run. The key input will be the `event_catalog`, which will be set to `["vacancy_hop"]`.

The user will launch the `cdd run-kmc` command. The primary acceptance criterion is the successful completion of the kMC run for the specified number of events. After the run, the system should output a final structure and a log of the vacancy's position over time. The user will visualize the initial and final structures. In the final structure, the vacancy must have moved from its original position. By plotting the vacancy's position over simulation time, the user should observe a random walk trajectory. From this trajectory, they can calculate a diffusion coefficient and compare it to known experimental or theoretical values for Aluminum at that temperature. A reasonable agreement would provide strong validation that the kMC engine is correctly simulating the physics of rare events using the MLIP.

---

### **Scenario UAT-C08-02: Experience clear progress indicators during a long workflow**

**(Min 300 words)**
This scenario focuses entirely on the user experience during a long, multi-stage process. A new user is starting a full active learning workflow from scratch. In previous versions, this might have involved long periods of waiting with no feedback. This UAT will verify that the new UI/UX improvements provide a clear, informative, and reassuring experience.

The user will launch a full `cdd active-learn` run on a new system. They will not focus on the scientific output, but on the terminal output itself. The acceptance criteria are based on a checklist of expected UI elements.
1.  **Labelling:** When the `LabelingEngine` is processing a batch of, say, 20 structures, the user **must** see a progress bar. The bar should clearly show the progress, like `Labeling Structures [██████----] 10/20`.
2.  **Training:** When the `TrainingEngine` is running, which can be slow, the user **must** see a spinner or an animated progress bar with a message like `Training MLIP model...`.
3.  **Simulation:** During the main MD or kMC simulation, which is the longest step, the user **must** see a progress bar indicating the number of simulation steps or events completed, e.g., `Running Simulation [███-------] 30,000 / 100,000 steps`.

The user should feel constantly informed about what the system is doing and how far along it is. The absence of long, silent pauses in the terminal output is the primary metric of success for this scenario. A positive user experience here is crucial for the adoption and usability of the tool.

---

### **Scenario UAT-C08-03: Receive helpful and well-formatted output and error messages**

**(Min 300 words)**
This scenario tests the clarity and professionalism of the system's logging and error reporting. The user will evaluate both a successful run and a deliberately triggered failure.

First, the user will perform a successful run and inspect the log file or terminal output. The acceptance criteria are related to the quality of the output. The logs should be well-structured, with clear timestamps and log levels (e.g., INFO, WARNING). Important summary information, like the final trained model accuracy or the number of trapped structures, should be highlighted (e.g., in bold or a different color) to be easily visible. The overall output should look clean and professional, not like a raw dump of text.

Second, the user will deliberately introduce an error. For example, they could point the configuration to a non-existent Quantum Espresso executable. They will then run the pipeline. The acceptance criterion is how the system handles the error. It **must not** crash with a long, cryptic Python traceback. Instead, the user should see a clearly formatted, user-friendly error message printed to the console, perhaps highlighted in red. The message should be helpful, for example: `[ERROR] Quantum Espresso executable not found at path: '/usr/bin/nonexistent/pw.x'. Please check your configuration file and ensure the path is correct.` This graceful handling of errors is a hallmark of a mature and user-friendly application.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C08-01 - Successfully run a kMC simulation of vacancy diffusion**

*   **GIVEN** a high-quality, pre-trained MLIP for Aluminum.
*   **AND** an `ase.Atoms` object representing a large Al supercell with a single vacancy.
*   **AND** a configuration file specifying a kMC simulation at 500 K for 10,000 events, with the event type `vacancy_hop`.
*   **WHEN** the user executes the command `cdd run-kmc`.
*   **THEN** the `SimulationEngine` should load the MLIP and the initial structure.
*   **AND** it should identify all possible vacancy hop events from the initial state.
*   **AND** it should start the main kMC loop, showing a progress bar for the 10,000 events.
*   **AND** the kMC loop should complete successfully.
*   **AND** the system should output a final `ase.Atoms` object.
*   **AND** when the initial and final structures are compared, the position of the vacancy must be different.
*   **AND** the system should produce a log file tracing the vacancy's position as a function of simulation time, allowing for analysis of its diffusive random walk.

---

**Scenario: UAT-C08-02 & UAT-C08-03 - A polished and informative user experience**

*   **GIVEN** a user is about to start a long, multi-stage workflow like `cdd active-learn`.
*   **WHEN** the workflow begins the labelling stage for 50 structures.
*   **THEN** a progress bar from the `rich` library must be displayed in the terminal, showing the fraction of structures that have been labelled.
*
*   **GIVEN** the workflow proceeds to the training stage.
*   **WHEN** the `TrainingEngine` is fitting the model.
*   **THEN** a spinner animation must be displayed with an informative message like "Fitting ACE model...".
*
*   **GIVEN** the user has made a mistake in their configuration file, such as an invalid path.
*   **WHEN** the user runs the application.
*   **THEN** the application must not crash with a raw Python stack trace.
*   **AND** it must instead print a clear, color-coded error message that explains the problem in simple terms.
*   **AND** it must guide the user on how to fix the problem (e.g., "Please check the 'command' field in your config file.").
*   **AND** all informational logs printed during a successful run should be cleanly formatted, with timestamps and log levels.
