# User Acceptance Test (UAT): Cycle 5

**Version:** 1.0.0
**Status:** Ready for Testing

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 5. The focus is on the final features and polish of the MLIP-AutoPipe system. This includes verifying the advanced Adaptive Kinetic Monte Carlo (akMC) simulation mode, ensuring the performance optimizations are effective, and confirming that the user interface and documentation provide a high-quality user experience.

| Scenario ID | Test Scenario Description                                                                | Priority |
|-------------|------------------------------------------------------------------------------------------|----------|
| UAT-C5-001  | **Verify Successful akMC Run for a Rare Event:** The user configures the system to simulate vacancy diffusion in a crystal using the `akMC` mode. The system must successfully identify the saddle point for the diffusion hop, calculate a reasonable energy barrier, and execute a KMC step. | High     |
| UAT-C5-002  | **Verify Active Learning During Saddle Point Search:** The user provides an MLIP that is known to be poor at describing the transition state for an event. During the akMC saddle point search, the system must trigger the uncertainty-based re-training loop, refine the potential, and then successfully converge on the correct saddle point. | High     |
| UAT-C5-003  | **Verify KMC Performance with Numba:** The user runs an akMC simulation with a large number of possible events. The KMC step, which involves selecting an event from the list, should execute very quickly (milliseconds), demonstrating the effectiveness of the Numba JIT compilation. | Medium   |
| UAT-C5-004  | **Verify Polished CLI Experience:** The user runs the main pipeline command. The command-line interface should display informative, real-time feedback, such as progress bars for long tasks and clearly formatted summary panels, using the `rich` library. | High     |
| UAT-C5-005  | **Verify Documentation and Tutorials:** The user follows the tutorial in the `docs` directory to run a complete workflow for an example system. The instructions must be clear, correct, and sufficient for a new user to successfully run the software and understand the output. | High     |


## 2. Behavior Definitions

This section provides detailed Gherkin-style behavior definitions for the test scenarios, outlining the exact steps and expected outcomes.

---

### **Scenario: UAT-C5-001 - Verify Successful akMC Run for a Rare Event**

*   **GIVEN** an `input.yaml` configured with `simulation_mode: "akmc"`.
*   **AND** the system is a crystal structure with a single vacancy.
*   **AND** a well-trained MLIP for this system is provided.
*   **WHEN** the user launches the simulation.
*   **THEN** the system should initiate a saddle point search from the initial state.
*   **AND** the logs should show that the saddle point search has converged on a transition state.
*   **AND** the system should calculate and log the energy barrier for the event (e.g., "Found event: Vacancy Hop with barrier 0.65 eV"). The barrier should be scientifically plausible.
*   **AND** the system should then perform a KMC step, and the logs should state that the simulation time has been advanced (e.g., "KMC step taken. Time advanced by 1.2 microseconds.").
*   **AND** the final atomic structure should show that the vacancy has moved to an adjacent lattice site.

---

### **Scenario: UAT-C5-002 - Verify Active Learning During Saddle Point Search**

*   **GIVEN** an `input.yaml` configured for an `akmc` simulation.
*   **AND** the provided initial MLIP was trained only on perfect crystal structures, with no information about vacancies or saddle points.
*   **WHEN** the user launches the simulation.
*   **THEN** the system will start a saddle point search.
*   **AND** during the search, as the atoms move towards the transition state, the uncertainty trigger must be activated.
*   **AND** the system must pause the saddle point search, extract the high-uncertainty structure, run DFT labeling, and re-train the potential.
*   **AND** after re-training, the system must automatically resume the saddle point search from where it left off, now with the improved potential.
*   **AND** the search should then converge to the correct saddle point.
*   **AND** the user should be able to see the new saddle point structure added to the training set in the log output.

---

### **Scenario: UAT-C5-003 - Verify KMC Performance with Numba**

*   **GIVEN** an akMC simulation is running and has identified a state with 100+ different possible escape events.
*   **WHEN** the system needs to perform the KMC step to select one of these events.
*   **THEN** the log message for the KMC step selection should appear almost instantaneously after the rates are calculated.
*   **AND** profiling the code should confirm that the Numba-jitted function for the KMC selection loop is responsible for this high performance.
*   **AND** the user should be able to see a progress bar for the KMC simulation, indicating that the process is running and not stalled.

---

### **Scenario: UAT-C5-004 - Verify Polished CLI Experience**

*   **GIVEN** any standard workflow run from the command line.
*   **WHEN** the `LabelingEngine` is running its loop over structures to calculate with DFT.
*   **THEN** the user must see a progress bar in their terminal showing the progress (e.g., `[ 10/50 ]`).
*   **WHEN** the configuration is first expanded at the beginning of a run.
*   **THEN** the user must see a formatted panel in the terminal that neatly summarizes the key configuration parameters being used for the run.
*   **AND** all log messages should be clearly legible and well-formatted.
*   **AND** the CLI should provide helpful error messages if the user provides an invalid input file or command-line argument.

---

### **Scenario: UAT-C5-005 - Verify Documentation and Tutorials**

*   **GIVEN** a new user who has just cloned the repository.
*   **WHEN** the user opens the `docs/user_guide.md` file.
*   **THEN** the document should provide clear, step-by-step instructions on how to install the software and its dependencies.
*   **WHEN** the user navigates to the `docs/tutorials/si_bulk_example.md` file.
*   **AND** the user follows the instructions in the tutorial exactly, copying and pasting the provided `input.yaml` content and running the specified commands.
*   **THEN** the software must run successfully and produce results that match the expected output described in the tutorial.
*   **AND** the tutorial must explain the meaning of the results and how to interpret them.
*   **AND** the user should be able to successfully run all the examples in the documentation without encountering any errors or inconsistencies.
