# UAT.md: Cycle 05 - Advanced Simulations & User Interface

## 1. Test Scenarios

User Acceptance Testing for the final development cycle, Cycle 05, is designed to evaluate the project as a complete, finished product. The UAT focuses on the total user experience and the delivery of the promised advanced scientific capabilities. For this cycle, the user will interact with the system exclusively through the new, formal command-line interface. Their expectation is that this interface will be professional, intuitive, robust, and self-documenting. Furthermore, they expect that the new, advanced Kinetic Monte Carlo (kMC) simulation feature will work as advertised, enabling them to confidently study the long-timescale phenomena that motivated the feature's development.

| Scenario ID | Scenario Description                                       | Priority |
| :---------- | :--------------------------------------------------------- | :------- |
| UAT-C05-01  | **Verification of a Professional CLI Experience for the Standard End-to-End Workflow**      | High     |
|             | A user's primary goal is to run a standard, end-to-end active learning workflow, from a minimal config file to a final trained model, using only the new command-line interface. The user expects a smooth, professional, and transparent experience. This includes being able to provide a simple, intuitive command, seeing clear and informative log messages that indicate the pipeline's progress, and receiving a final, trained model as the successful output. The user should not need to interact with any Python scripts directly. This UAT is critical as it validates the primary user-facing entry point and the overall usability of the entire application. It confirms that the complexity of the underlying pipeline has been successfully abstracted away from the user. |          |
| UAT-C05-02  | **Successful Execution and Verification of a Long-Timescale kMC Simulation**           | High     |
|             | A user, having already generated a robust MLIP using the main `run` command, now wants to leverage this potential to investigate a specific long-timescale phenomenon, such as a diffusion process in a crystal. They will use the new, dedicated `run-kmc` command. The user expects the system to successfully launch a kinetic Monte Carlo simulation, to be able to observe atomic "hops" or other rare events in the output trajectory, and, critically, to see the simulation time advance in the non-linear, event-based manner that is the hallmark of a kMC simulation. This UAT validates the delivery of the key advanced scientific feature of this cycle. |          |
| UAT-C05-03  | **Confirmation of CLI Robustness, User-Friendliness, and Help System**                         | Medium   |
|             | A new user, who is unfamiliar with the tool, wants to learn how to use it directly from the command line. They expect the CLI to have a comprehensive, built-in help system that provides clear and sufficient instructions for its use. They also expect the tool to be robust against common user errors, providing helpful and informative error messages if they misuse a command, provide an invalid argument, or point to a non-existent file. This UAT is essential for ensuring that the tool is accessible and user-friendly for beginners, which is crucial for its adoption by the wider scientific community. |          |

## 2. Behavior Definitions

These behaviors will be tested entirely from the user's command line. The user will execute the `mlip-pipe` command with various arguments and flags, and will verify the outcomes by observing the console output and inspecting the generated files.

### Scenario: UAT-C05-01 - Verification of a Professional CLI Experience for the Standard End-to-End Workflow

*   **GIVEN** a clean directory containing only a single, minimal `input.yaml` file for a simple system like Silicon.
*   **WHEN** the user opens their terminal and types the single, simple command: `mlip-pipe run input.yaml`.
*   **THEN** the application must start immediately and display a clear, well-formatted startup message that includes the application name and version number (e.g., "Starting MLIP-AutoPipe v1.0.0").
*   **AND** as the pipeline executes, the console output must provide a continuous stream of informative, human-readable status updates, clearly delineating the major stages of the workflow. The user should see messages like: "Configuration expanded successfully," "Starting surrogate model exploration phase," and "Entering active learning loop, iteration 1 of 5."
*   **AND** the process must complete successfully after the configured number of learning iterations.
*   **AND** upon successful completion, a final, clear message must be displayed in the console, for example: "Workflow complete. Final MLIP saved to `output/final_model.pt`".
*   **AND** the user must be able to navigate to the `output/` directory and find the correctly named, trained model file, confirming the successful completion of the entire process.

### Scenario: UAT-C05-02 - Successful Execution and Verification of a Long-Timescale kMC Simulation

*   **GIVEN** a previous, completed active learning run has produced a robust, well-trained MLIP for a system with a known diffusion mechanism (e.g., a crystal containing a single atomic vacancy).
*   **AND** the user has an `input.yaml` file that is correctly configured for a kMC simulation, pointing to the pre-trained MLIP.
*   **WHEN** the user executes the specific command for this task: `mlip-pipe run-kmc input.yaml`.
*   **THEN** the system must load the specified pre-trained MLIP and successfully start the Kinetic Monte Carlo simulation.
*   **AND** the log output generated by this run must be clearly and fundamentally distinct from a normal MD run. The user should observe a series of messages that reflect the event-based nature of the kMC algorithm, such as: "kMC step 1: Event search started...", "Found 5 potential events.", "Executing event: atom 12 hops to vacancy site.", and critically, a message indicating the non-linear advancement of time: "System time advanced by 1.2 nanoseconds."
*   **AND** upon completion of the simulation, the user must be able to load the output trajectory file into a visualiser and confirm that one or more rare events (like the vacancy moving through the crystal lattice) have indeed occurred. This provides definitive proof that the simulation has explored the long-timescale dynamics that would not have been accessible in a short MD run.

### Scenario: UAT-C05-03 - Confirmation of CLI Robustness, User-Friendliness, and Help System

*   **GIVEN** the user has the `mlip-pipe` tool correctly installed in their environment.
*   **WHEN** the user, acting as a beginner, types the command `mlip-pipe --help` in their terminal.
*   **THEN** the CLI must immediately print a well-structured and comprehensive help message. This message must list all available subcommands (e.g., `run`, `run-kmc`) and provide a brief, clear description of what each one does.
*   **AND WHEN** the user requests help for a specific subcommand, for example, `mlip-pipe run --help`.
*   **THEN** the CLI must print a more detailed help message specifically for the `run` command. This message must clearly explain the required arguments (e.g., `CONFIG_FILE`) and list any available command-line options with a description of their function.
*   **AND WHEN** the user deliberately tries to run a command with a typographical error, for example, `mlip-pipe rnu input.yaml`.
*   **THEN** the CLI must not crash with an unhandled exception. Instead, it must display a helpful, user-friendly error message that suggests the correct command, for example: "Error: No such command 'rnu'. Did you mean 'run'?".
*   **AND WHEN** the user tries to execute the `run` command but provides a path to a configuration file that does not exist.
*   **THEN** the CLI must again fail gracefully and provide a clear, actionable error message, such as: "Error: Invalid value for 'CONFIG_FILE': Path 'nonexistent_file.yaml' does not exist."
