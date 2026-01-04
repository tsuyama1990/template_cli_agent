# User Acceptance Test Plan: MLIP-AutoPipe Cycle 1

This document outlines the User Acceptance Testing (UAT) plan for Cycle 1 of the MLIP-AutoPipe project. The goal of this UAT is to verify that the core command-line pipeline is functional, reliable, and meets the essential requirements of a user aiming to generate a basic MLIP training dataset.

## 1. Test Scenarios

This UAT is designed to be executed by a user from the command line. The primary artefact for this UAT will be a Jupyter Notebook (`UAT_Cycle1.ipynb`), which will guide the user through the process, providing code cells to execute commands and explanations of the expected outcomes. This approach allows for an interactive and educational testing experience, serving as a tutorial for new users.

| Scenario ID | Scenario Name                     | Priority | Description                                                                                                                                                                                                            |
| :---------- | :-------------------------------- | :------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| UAT-C1-01   | **Successful End-to-End Run**     | High     | This is the primary "happy path" scenario. The user will configure a simple binary alloy system, run the entire pipeline with a single command, and verify that a valid, non-empty ASE database is produced as the final output. This test validates the core functionality of the entire system. |
| UAT-C1-02   | **Configuration Validation**      | High     | This scenario tests the system's robustness against invalid user input. The user will attempt to run the pipeline with deliberately incorrect configuration files (e.g., invalid chemical symbols, negative temperature). The goal is to verify that the system fails gracefully and provides clear, informative error messages. |
| UAT-C1-03   | **Reproducibility Check**         | Medium   | This scenario verifies the scientific integrity of the pipeline. The user will run the exact same pipeline twice, ensuring that the use of a fixed random seed results in an identical final database. This is crucial for reproducible research. |
| UAT-C1-04   | **Inspecting the Output Database**| High     | This scenario focuses on the quality and correctness of the final output. The user will load the generated ASE database and inspect its contents programmatically. They will verify that the number of structures is as expected and that each structure contains the necessary metadata (energy, forces). |


## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behaviour for each test scenario. These will be implemented and documented within the `UAT_Cycle1.ipynb` notebook.

---

### **Scenario: UAT-C1-01 - Successful End-to-End Run**

(Minimum 500 words)

This scenario represents the most critical user journey for Cycle 1. It validates that a user can successfully configure a standard materials science problem, execute the automated pipeline, and receive a valid output without any errors. The goal is to simulate a real-world use case from start to finish. The user, likely a researcher or student in computational materials science, wants to generate a small training set for a Silicon-Germanium (SiGe) alloy, a common semiconductor material. They are starting with a simple 50/50 composition in a diamond lattice structure. They want to run a short Molecular Dynamics (MD) simulation at a high temperature to create a set of diverse, amorphous-like structures. This test ensures that all the individual components of the pipeline—Generation, Exploration, Sampling, and Storage—are correctly integrated and functioning as a cohesive unit. A successful outcome here provides high confidence in the overall system's capability to deliver on its core promise. The user should be able to achieve this by simply preparing a standard set of configuration files and executing a single command in their terminal. The process should feel seamless and automated, truly embodying the project's philosophy of "removing the human expert from the loop" for standard tasks. The test will be considered a pass if the command completes with a zero exit code and the specified output files, most importantly the final database, are created in the output directory.

**GIVEN** a directory containing valid configuration files for a SiGe binary alloy system.
**AND** the configuration specifies generating 2 initial structures.
**AND** the configuration specifies running a 10-step MD simulation at 1500K.
**AND** the configuration specifies randomly sampling 5 structures for the final dataset.
**AND** a valid MACE model file is available at the path specified in the configuration.

**WHEN** the user executes the `mlip-autopipec run` command from the terminal, pointing to the configuration directory.

**THEN** the command should execute without raising any errors and exit with a status code of 0.
**AND** an output directory should be created.
**AND** inside the output directory, an ASE database file named `final_structures.db` must exist.
**AND** the `final_structures.db` file must be a valid SQLite database and be larger than zero bytes.
**AND** log messages should be printed to the console, indicating the progress of each stage (Generation, Exploration, Sampling, Storage).

---

### **Scenario: UAT-C1-02 - Configuration Validation**

(Minimum 500 words)

A robust application must handle incorrect user input gracefully. This scenario is designed to test the resilience and user-friendliness of the MLIP-AutoPipe's configuration system. The Pydantic-based schema validation is a cornerstone of our design, and this test verifies that it is working as intended from the user's perspective. A user, especially one new to the system, will inevitably make mistakes when writing their YAML configuration files. They might make a typo in an element name, specify a negative simulation temperature, or provide a file path that does not exist. A poorly designed system would either crash with a cryptic traceback or, even worse, proceed with the invalid input, leading to nonsensical results. This UAT scenario ensures that our system does neither. It validates that the CLI provides immediate, clear, and actionable feedback when it detects a problem in the configuration. This is a critical aspect of the user experience. By testing multiple common error types, we ensure the validation is comprehensive. A successful outcome for this test is, paradoxically, a failed execution of the pipeline, but a successful validation of the input. The user should feel confident that the system protects them from simple mistakes, rather than feeling frustrated by unhelpful error messages. This builds trust in the tool and its results.

**GIVEN** a directory with a configuration file containing an invalid element symbol (e.g., "Xx") in the elements list.
**WHEN** the user executes the `mlip-autopipec run` command.
**THEN** the application should fail to start.
**AND** it must print a user-friendly error message to the console that clearly indicates the configuration error (e.g., "ValidationError: 'Xx' is not a valid chemical symbol").
**AND** the process must exit with a non-zero status code.

**GIVEN** a directory with a configuration file containing a negative value for the MD temperature (e.g., `temperature: -100`).
**WHEN** the user executes the `mlip-autopipec run` command.
**THEN** the application should fail to start.
**AND** it must print a user-friendly error message indicating that the temperature must be a positive value.
**AND** the process must exit with a non-zero status code.

**GIVEN** a directory with a configuration file pointing to a non-existent MLIP model file.
**WHEN** the user executes the `mlip-autopipec run` command.
**THEN** the application should fail to start.
**AND** it must print a user-friendly error message indicating that the file path for the model could not be found.
**AND** the process must exit with a non-zero status code.

---

### **Scenario: UAT-C1-03 - Reproducibility Check**

(Minimum 500 words)

Reproducibility is a fundamental pillar of computational science. This UAT scenario is designed to verify that the MLIP-AutoPipe pipeline is deterministic. For a given set of inputs and a fixed random seed, the output must be identical every single time. This is crucial for researchers who need to be able to reproduce their own results, and for others who may wish to validate or build upon their work. The stochasticity in the pipeline comes from several sources: the initial random placement of atoms in an alloy, the random velocities assigned to atoms at the start of an MD simulation, and the random selection of structures during the sampling phase. By setting a global random seed (for both Python's `random` module and `numpy.random`), we should be able to control all these sources of randomness and make the entire workflow deterministic. This test will involve running the exact same pipeline twice, back-to-back. The first run will produce a database. The second run will do the same. The test will then programmatically compare the two databases to ensure they are bit-for-bit identical. A successful test provides strong evidence that the pipeline is behaving in a scientifically rigorous manner. It gives the user confidence that their results are not a fluke of a particular random run and can be reliably reproduced in the future.

**GIVEN** a valid configuration for a SiGe alloy.
**AND** the configuration specifies a fixed integer for the global random seed.
**AND** an empty output directory named `run_1`.

**WHEN** the user executes the `mlip-autopipec run` command, directing the output to `run_1`.

**THEN** the run should complete successfully, producing a `final_structures.db` file in the `run_1` directory.

**GIVEN** the exact same configuration as the first run.
**AND** a new empty output directory named `run_2`.

**WHEN** the user executes the `mlip-autopipec run` command again, directing the output to `run_2`.

**THEN** the second run should also complete successfully, producing a `final_structures.db` file in the `run_2` directory.
**AND** a script that compares the file hashes (e.g., SHA256) of `run_1/final_structures.db` and `run_2/final_structures.db` must show that the hashes are identical.
