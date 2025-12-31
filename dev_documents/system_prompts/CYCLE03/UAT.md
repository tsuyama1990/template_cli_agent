# UAT.md: Cycle 03 - Efficient Exploration & Optimisation

## 1. Test Scenarios

User Acceptance Testing for Cycle 03 is designed to build the user's confidence in the "intelligence" and efficiency of the pipeline. The user, who is likely a technically-astute researcher or scientist, needs to be certain that the system is not operating as a brute-force tool that naively runs expensive DFT calculations on random or redundant atomic structures. Instead, they expect to see clear evidence that the system is strategically using the fast surrogate model to explore a wide range of possibilities and then intelligently "down-selecting" only the most unique and scientifically valuable structures for labelling. The computational performance of this crucial exploration and selection stage is also a key user concern and a primary focus of this UAT.

| Scenario ID | Scenario Description                                       | Priority |
| :---------- | :--------------------------------------------------------- | :------- |
| UAT-C03-01  | **Clear Verification of Surrogate-Model Driven Exploration Phase**     | High     |
|             | The user's primary goal in this scenario is to confirm that the system correctly executes the exploration phase using the fast surrogate model *before* any expensive calculations are performed. When initiating a workflow for a simple system, the user expects to see clear, unambiguous log messages indicating that a Molecular Dynamics (MD) simulation is being run with the MACE surrogate model. This simulation should generate a large number of candidate structures. Crucially, the user must be able to verify that no slow, expensive DFT calculations are launched during this phase. This test confirms that the exploration phase is correctly integrated into the workflow and is functioning as the primary, low-cost mechanism for generating a diverse pool of candidate data. |          |
| UAT-C03-02  | **Tangible Evidence of Diverse and Intelligent Structure Sampling (DIRECT)**        | High     |
|             | This scenario focuses on the outcome of the DIRECT sampling algorithm. After the exploration phase is complete, the user will inspect the specific set of structures that the system has selected for DFT labelling. The user's expectation is that this set will be both small (a fraction of the total trajectory) and structurally diverse. For a simulation designed to model the melting of a crystal, for example, the user expects to see not just perfect or near-perfect crystal structures, but also highly distorted, high-energy, and disordered liquid-like configurations among the selected candidates. This provides tangible proof that the DIRECT sampling algorithm is effective at identifying and capturing more than just the obvious, low-energy equilibrium states, which is essential for building a robust and transferable potential. |          |
| UAT-C03-03  | **Confirmation of Acceptable Performance of the Exploration Phase**        | Medium   |
|             | In this scenario, the user runs an exploration phase on a moderately sized system to validate its real-world performance. While it is understood that this is a computationally intensive step, the user expects that the descriptor calculation and clustering for a large trajectory (e.g., thousands of frames) will complete in a reasonable amount of time (e.g., minutes, not hours) on a standard workstation. The user also expects to see efficient use of their machine's resources during this time. This UAT is designed to validate that the Numba optimisations are effective in practice and that the exploration phase does not represent an intractable bottleneck that would make the entire pipeline impractical for daily use. |          |

## 2. Behavior Definitions

These behaviors will be tested by a user running the pipeline's exploration phase from the command line and then carefully inspecting the system's log files for specific messages, querying and visualising the contents of the database, and monitoring system resource usage.

### Scenario: UAT-C03-01 - Clear Verification of Surrogate-Model Driven Exploration Phase

*   **GIVEN** a project that has already successfully completed the initialization phase (as in Cycle 02), resulting in a database populated with a few seed structures.
*   **AND** the `exec_config_dump.yaml` file for the project explicitly specifies the use of the `mace_mp` surrogate model for the exploration phase.
*   **WHEN** the user launches the next stage of the pipeline from the command line (e.g., via a command like `uv run mlip-pipe run-exploration ...`).
*   **THEN** the user must be able to monitor the system's log files in real-time or inspect them after the run.
*   **AND** the logs must contain a clear, unambiguous message indicating the start of the surrogate model simulation, such as "Starting surrogate MD simulation with MACE calculator".
*   **AND** the logs should subsequently show the real-time progress of the MD simulation, for example by printing out the current timestep, temperature, and potential energy at regular intervals.
*   **AND** most critically, the user must be able to verify by inspecting the entire log file for this phase that there are absolutely no messages indicating the start or completion of a Quantum Espresso calculation. This confirms that the exploration is being performed without any DFT cost.

### Scenario: UAT-C03-02 - Tangible Evidence of Diverse and Intelligent Structure Sampling (DIRECT)

*   **GIVEN** a workflow has been configured and run for a simple crystal material (like Silicon) with a simulation temperature set to be just above its known melting point, encouraging a phase transition.
*   **AND** the exploration phase (as verified in UAT-C03-01) has completed successfully.
*   **WHEN** the user connects to the `mlip.db` database using a database client and extracts all the structures whose state is marked as 'selected_for_labelling'.
*   **AND** the user then visualises this specific subset of structures using a standard molecular visualisation tool like `ase gui`.
*   **THEN** the user must observe that the set of structures is not uniform or redundant.
*   **AND** the user should be able to clearly identify distinct structural motifs within the selection, including:
    *   Some structures that are still close to the perfect, original diamond-cubic crystal lattice of Silicon.
    *   Some structures that show significant thermal distortion, with atoms displaced far from their equilibrium positions, representing high-energy configurations.
    *   Some structures that are clearly non-crystalline, disordered, or amorphous, representing the molten, liquid state of the material.
*   **AND** the user will verify that the total number of these selected structures is a small and reasonable fraction (e.g., less than 1%) of the total number of frames that were generated in the initial surrogate MD trajectory, confirming the "funnelling" effect.

### Scenario: UAT-C03-03 - Confirmation of Acceptable Performance of the Exploration Phase

*   **GIVEN** a surrogate MD simulation has been run and has produced a trajectory file containing a substantial number of frames (e.g., 10,000 frames for a 64-atom system).
*   **WHEN** the user initiates the DIRECT sampling part of the exploration phase, which involves descriptor calculation and clustering.
*   **THEN** the system log must print clear, timestamped messages indicating the start and end of the most intensive steps, for example: "Starting SOAP descriptor calculation..." and "SOAP descriptor calculation finished. Elapsed time: X minutes".
*   **AND** the user will observe that the elapsed time for these steps is reasonable and practical for the given input size (e.g., for a 10,000-frame, 64-atom system, the user might expect this process to complete in under 30 minutes on a modern multi-core workstation).
*   **AND** while the descriptor calculation is running, the user will use standard system monitoring tools (like `top`, `htop`, or Windows Task Manager) to observe the resource usage of the Python process. The user should see that the process is efficiently utilising the available CPU resources, for example, by observing a CPU usage percentage significantly greater than 100% on a multi-core machine. This provides direct evidence that the Numba JIT compilation and parallelisation are working effectively.
