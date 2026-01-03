# User Acceptance Testing (UAT): CYCLE01

This document outlines the User Acceptance Testing scenarios for the features implemented in CYCLE01 of the MLIP-AutoPipe project. The purpose of UAT is to validate the system from an end-user's perspective, ensuring that the software is not only technically functional but also usable, robust, and fit for its intended purpose. For CYCLE01, the goal is to verify that the core command-line pipeline provides a reliable and effective tool for a materials scientist to generate a high-quality training dataset for a multi-component alloy. These tests are designed to simulate real-world usage and to confirm that the key features of the system meet the user's requirements and expectations.

## 1. Test Scenarios

The following scenarios have been designed to cover the critical functionality delivered in CYCLE01. They test the primary success path, the system's robustness to interruption, its handling of user error, and its performance characteristics.

| Scenario ID | Title                                       | Priority | Description                                                                                                                                                                                                                                                                                                                                                                                                                                    |
|-------------|---------------------------------------------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| UAT-C01-001 | End-to-End Pipeline Run for a Binary Alloy  | High     | This scenario represents the single most important validation of the core functionality. It is the "golden path" test that simulates a standard, successful user interaction from start to finish. The user will configure the system for a simple, well-understood binary alloy (e.g., Copper-Gold) and execute the entire four-stage pipeline using a single command. The primary goal is to verify that all the individual components—generation, exploration, sampling, and storage—integrate seamlessly and correctly to produce the desired final output. Success is defined by the creation of a valid, richly populated ASE database containing the final sampled structures, each with its associated energy and force data. This test provides the ultimate confidence that the system as a whole is working as designed. |
| UAT-C01-022 | Verify Pipeline Resumption (Checkpointing)  | Medium   | This scenario is designed to test the pipeline's robustness, intelligence, and efficiency in the face of interruption, which is a common occurrence in long-running scientific workflows. A user will initiate a standard pipeline run and then deliberately interrupt it after the first stage (Generation) has completed. They will then restart the exact same command. The expected outcome is that the pipeline intelligently detects the pre-existing output from the generation stage and correctly skips it, resuming the workflow from the beginning of the exploration stage. This test is crucial for verifying that the checkpointing logic is functional, as this feature is what saves significant computation time and prevents the catastrophic loss of data in case of an unexpected shutdown. |
| UAT-C01-003 | Configuration Validation and Error Handling | High     | This scenario focuses on the system's robustness and user-friendliness when confronted with user error. The user will intentionally try to run the pipeline with a syntactically correct but logically invalid configuration file. Examples of such errors include specifying a composition where the fractions do not sum to 1.0, or misspelling a key parameter. The system must not crash or produce cryptic error messages. Instead, it is expected to gracefully fail at the very beginning of the process, printing a clear, informative, and human-readable error message that precisely identifies the specific problem in the configuration file and explains why it is wrong. This demonstrates the effectiveness of the Pydantic-based validation and is a critical feature for a positive user experience. |
| UAT-C01-004 | Verify Parallel Execution (`max_workers`)   | Medium   | This scenario validates the performance and scalability aspect of the computationally intensive exploration stage. The user will configure and run the same simulation workflow twice. The first run will be in serial mode, with `max_workers` set to 1. The second run will be in parallel mode, with `max_workers` set to a higher number, such as 4. The expectation is that the second run completes in a significantly shorter amount of time. Furthermore, the user should be able to verify (using system monitoring tools like `top` or `htop`) that multiple Python processes are indeed running concurrently during the exploration stage of the second run. This test confirms that the parallel processing logic is functional and that the tool can effectively leverage multi-core hardware. |

**UAT Approach using a Jupyter Notebook:**

To ensure these tests are easy for the user to execute, understand, and verify, a comprehensive Jupyter Notebook named `CYCLE01_UAT.ipynb` will be provided as the primary testing artifact. This approach moves beyond simple command-line instructions to create an interactive, self-documenting tutorial and testing environment. The notebook will be structured with clear headings for each test scenario and will contain a mixture of markdown explanations and executable code cells that will:
1.  Programmatically define and write the necessary Hydra configuration files (`config.yaml`) to disk for each specific scenario. This ensures the tests are perfectly reproducible.
2.  Use the `!` shell command syntax (e.g., `!mlip_autopipec run-pipeline --config-path ...`) to execute the pipeline directly from within the notebook, capturing the command's output in the cell.
3.  Contain Python code cells that use the `ase.db` library to connect to the generated `results.db` file, query its contents, perform assertions (e.g., `assert db.count() == 100`), and print detailed information about the database content.
4.  Display the results in a clear and understandable format, printing a bold "PASS" or "FAIL" message for each scenario based on the assertions.

This notebook-driven approach makes the UAT process transparent, reproducible, and highly interactive. It empowers the user to not just verify the results but also to easily inspect the outputs and even tweak the configurations to experiment further, turning the UAT process from a simple check into a valuable learning and exploration tool.

## 2. Behavior Definitions

The following are the detailed, Gherkin-style behavior definitions for the user acceptance test scenarios. They describe the system's expected behavior from the user's point of view in a clear, unambiguous format.

---

### **Scenario: UAT-C01-001 - End-to-End Pipeline Run for a Binary Alloy**

This test verifies the successful, "happy path" execution of the entire pipeline, confirming that all stages are correctly integrated and produce a valid final product.

**GIVEN**
The user is in a clean directory that contains no previous output files.
And the user has created a valid Hydra configuration file for a Copper-Gold (CuAu) binary alloy system.
And this configuration specifies:
- In the `generator` section: create 10 initial structures of a 50/50 CuAu alloy in an FCC lattice.
- In the `explorer` section: run a short 50-step MD simulation at 300K, using the computationally cheap "EMT" calculator for speed.
- In the `sampler` section: randomly sample a total of 20 frames from all the resulting trajectories.
- In the `storage` section: save the final results to a database file named `CuAu_test.db`.

**WHEN**
The user executes the command `mlip_autopipec run-pipeline` in their terminal, pointing it to the directory containing the configuration file.

**THEN**
The command should execute to completion without printing any errors and should exit with a success status code (0).
And the user should see log messages indicating the progression through the four stages: Generation, Exploration, Sampling, and Storage.
And a file named `initial_structures.xyz` should be created in the output directory, and this file should contain exactly 10 atomic configurations.
And a set of trajectory files (e.g., `traj_0.xyz`, `traj_1.xyz`, ..., `traj_9.xyz`) should be created, one for each initial structure.
And a final database file named `CuAu_test.db` must be created.
And when the user inspects this `CuAu_test.db` database, it must contain exactly 20 rows, corresponding to the 20 requested samples.
And each row in the database, when inspected, must contain valid, non-null data for the structure's energy and atomic forces.

---

### **Scenario: UAT-C01-002 - Verify Pipeline Resumption (Checkpointing)**

This test verifies the critical robustness feature of the pipeline: its ability to be stopped and resumed without losing progress.

**GIVEN**
The user is in a clean directory.
And the user is using the same valid Hydra configuration file for the CuAu binary alloy as in the previous scenario.

**WHEN**
The user first executes the `mlip_autopipec run-pipeline` command.
And the user observes the log output and confirms that the Generation stage has completed and the file `initial_structures.xyz` has been created.
And the user then terminates the process (e.g., by pressing `Ctrl+C`) while the Exploration stage is still running.
And the user then immediately re-executes the exact same `mlip_autopipec run-pipeline` command.

**THEN**
The pipeline should restart its execution without any errors and should eventually run to completion.
And the log output from this second run must contain a clear message indicating that the "Generation" stage was skipped because its output (`initial_structures.xyz`) was found to already exist.
And the log output must show that the "Exploration" stage was started from the beginning.
And the final `CuAu_test.db` database should be created, and it must contain the correct number of 20 sampled structures, proving that all subsequent stages ran correctly after the resumption.

---

### **Scenario: UAT-C01-003 - Configuration Validation and Error Handling**

This test verifies that the system is robust against common user configuration errors and provides helpful, actionable feedback.

**GIVEN**
The user has created a Hydra configuration file for an alloy.
And this configuration file contains a deliberate and clear logical error. For example, in the `generator.params` section, the `composition` is set to `{"Cu": 0.4, "Au": 0.5}`, which is invalid because the fractions do not sum to 1.0.

**WHEN**
The user executes the command `mlip_autopipec run-pipeline` pointing to this invalid configuration file.

**THEN**
The command must fail immediately and must exit with a non-zero status code, indicating an error.
And the system must print a user-friendly error message to the console; it should not print a long, confusing Python stack trace.
And this error message must clearly state that the provided configuration is invalid.
And crucially, the message should be specific, mentioning that the `composition` field is the source of the error and explicitly explaining the reason for the failure (e.g., "sum of composition fractions must be 1.0").
And the user must be able to verify that no output files (such as `.xyz` files or a `.db` database) were created, confirming that the failure occurred early as intended.
