# User Acceptance Testing (UAT): Cycle 1

This document outlines the User Acceptance Testing scenarios for the deliverables of Cycle 1. The focus of this cycle is the core command-line functionality for generating an MLIP dataset. The UATs are designed to be executed by a user (e.g., a materials scientist) to confirm that the software meets its foundational requirements in a user-friendly and effective manner. These tests are not just about correctness, but also about the quality of the user experience, ensuring the tool is intuitive and provides clear, actionable feedback.

## 1. Test Scenarios

These scenarios represent the key user stories and success criteria for this cycle. They are designed to be run from the command line, using a text editor to create configuration files and a simple database browser (like `sqlite3` or `ase db`) to inspect the results. Each scenario is designed to build the user's confidence in the tool's reliability and functionality, moving from the ideal case to error handling and data validation.

| Scenario ID | Title                                       | Priority |
|-------------|---------------------------------------------|----------|
| UAT-C1-01   | Successful End-to-End Pipeline Run (Happy Path) | High     |
| UAT-C1-02   | Handling of Invalid Configuration File      | High     |
| UAT-C1-03   | Verifying Database Output Content         | Medium   |

### Scenario UAT-C1-01: Successful End-to-End Pipeline Run (Happy Path)

**Description:** This is the most critical UAT for this cycle, serving as the primary validation of the core value proposition. It verifies that a user can successfully and smoothly execute the entire data generation pipeline from a valid configuration file and receive a well-formed output database without any errors. This test is comprehensive, implicitly confirming that all the individual components of the pipeline—the configuration parser, the structure generator, the MD explorer, the random sampler, and the database storage module—are correctly integrated and functioning as a cohesive unit. The user experience is a key aspect of this test; it should be straightforward and transparent. The user should feel in control and informed throughout the process. Clear, concise feedback on the command line is expected, indicating the progress through each distinct stage of the pipeline. A successful run should build a strong sense of trust and reliability, encouraging the user to explore more complex systems. This test simulates the ideal user journey, where a well-defined scientific goal is translated into a configuration file, and the tool effortlessly produces the desired result, demonstrating the tool’s power and simplicity.

**User Story:** As a materials scientist, I want to define a simple binary alloy system in a YAML file and run a single command to generate a database of atomic structures, so that I can quickly create a foundational dataset for my MLIP training without needing to write complex scripts or perform manual data manipulation.

**Execution Steps:**
1.  The user creates a new YAML file named `config_valid.yaml` with plausible and valid content for a Copper-Gold alloy.
2.  The user executes the main CLI command: `mlip-autopipec run-pipeline --config config_valid.yaml`.
3.  The user observes real-time, sequential feedback on the command line, such as "Starting Generation stage...", "Generation complete.", "Starting Exploration stage...", etc. for all four stages.
4.  After the command finishes without any error messages, the user checks their working directory and confirms the presence of a new file named according to the configuration, `CuAu_dataset.db`.

### Scenario UAT-C1-02: Handling of Invalid Configuration File

**Description:** This UAT is essential for ensuring the system is robust, user-friendly, and resilient to common user errors. It verifies that the application provides clear, specific, and actionable feedback when the user provides a configuration file containing invalid or physically impossible parameters. A robust application should anticipate user mistakes and guide them towards a solution. The system must fail gracefully and informatively, rather than crashing with an obscure traceback or, even worse, proceeding with the invalid parameters and producing nonsensical results. This test is crucial for building user trust and making the tool easy to debug. The error message itself is the most important artifact of this test; it should be user-centric, avoiding internal jargon, and should pinpoint the exact parameter that is incorrect and explain why. For instance, `Error: Invalid configuration. 'exploration.temperature_k' must be a positive number.` is much more helpful than `pydantic.ValidationError`. A successful outcome of this test means the user can self-correct their configuration without needing to consult documentation or the source code, leading to a much smoother and less frustrating user experience.

**User Story:** As a materials scientist, when I inevitably make a mistake in my configuration file, such as a typo or an invalid value, I want the tool to immediately and clearly tell me exactly what is wrong, so that I can easily fix it and rerun the process without wasting time or computational resources.

**Execution Steps:**
1.  The user creates a deliberately invalid YAML file named `config_invalid.yaml`, for example, by setting the simulation temperature to a negative value.
2.  The user runs the CLI command: `mlip-autopipec run-pipeline --config config_invalid.yaml`.
3.  The user observes the terminal output. The application must not proceed with the pipeline. Instead, it must print a user-friendly error message that explicitly states that the `temperature_k` value is invalid and explains the constraint (e.g., it must be greater than zero).
4.  The user confirms that no part of the pipeline was executed by verifying that no output database file (e.g., `NiAl_dataset.db`) has been created.

### Scenario UAT-C1-03: Verifying Database Output Content

**Description:** This UAT focuses on the integrity and correctness of the final output, which is the most important product of the tool. It's not sufficient for the output file to simply exist; the user must have high confidence that the data contained within it accurately reflects the parameters they specified in their configuration. This test involves a "look inside the box" to inspect the generated database and verify its contents. The user will confirm that the number of stored structures precisely matches the number of requested samples. Furthermore, they will verify the chemical and structural properties of the data. For instance, the atomic composition of the structures in the database should match the composition defined in the system configuration. This test ensures that the data is not just present, but correct, and that the user's input parameters are being respected throughout the entire pipeline. Successfully passing this test gives the user the confidence needed to proceed with using the generated data for the computationally expensive task of training an MLIP, trusting that the input data is sound.

**User Story:** As a materials scientist, I need to be able to inspect the generated database and trust that it contains the correct number and type of structures that I requested, so that I can be completely confident in the quality and correctness of the data I use for training my machine learning models.

**Execution Steps:**
1.  The user first successfully completes the "happy path" scenario, UAT-C1-01, ensuring a valid database file exists.
2.  The user utilizes a standard, external tool like the `ase db` command-line utility to connect to and inspect the output file, for example: `ase db CuAu_dataset.db`.
3.  From the utility's output, the user verifies two key pieces of information:
    -   The total number of rows (which corresponds to the number of atomic structures) in the database is exactly what was specified in the `sampling.num_samples` parameter of the configuration file.
    -   By inspecting the details of a few rows, the user confirms that the chemical formula and the elements present in the structures match the `system.elements` and `system.composition` parameters.

## 2. Behavior Definitions

These Gherkin-style definitions provide a more formal, structured, and detailed description of the expected system behavior for each of the UAT scenarios. They are written in a way that can be understood by all stakeholders, including developers, testers, and end-users, and they serve as the unambiguous contract against which the application's functionality will be judged. Each definition follows the "Given-When-Then" structure to clearly separate the context, the action, and the expected, verifiable outcomes. This level of detail is crucial for ensuring that the implementation perfectly aligns with the user's expectations for a robust and reliable scientific tool. The definitions elaborate not just on the final state, but also on the intermediate feedback and the user experience, such as the nature of console messages and error reports, which are critical for a high-quality command-line application.

**Scenario: Successful End-to-End Pipeline Run**
```gherkin
GIVEN the user has created a valid configuration file named "config_valid.yaml"
  AND this file specifies a binary alloy of Copper ('Cu') and Gold ('Au') with a 50/50 composition
  AND the file requests the generation of 2 initial seed structures
  AND the file requests a Molecular Dynamics exploration at 500K for 100 steps
  AND the file specifies a random sampling method to select exactly 5 final structures
  AND the file defines the output database path as "CuAu_dataset.db"
WHEN the user executes the command "mlip-autopipec run-pipeline --config config_valid.yaml" in their terminal
THEN the application should start without any initial configuration errors
  AND the command should execute to completion and exit with a success code (0)
  AND during execution, the console should display clear, sequential status messages indicating the start and successful completion of the "Generation", "Exploration", "Sampling", and "Storage" pipeline stages
  AND upon completion, a new file named "CuAu_dataset.db" must exist in the user's current working directory
  AND this file should be a valid SQLite database, accessible by tools like the ASE database GUI or CLI.
```

**Scenario: Handling of Invalid Configuration File**
```gherkin
GIVEN the user has created an invalid configuration file named "config_invalid.yaml"
  AND this file contains a specific, deliberate error: the exploration temperature 'temperature_k' is set to a negative value (-150.0)
  AND the file defines the output database path as "NiAl_dataset.db"
WHEN the user executes the command "mlip-autopipec run-pipeline --config config_invalid.yaml"
THEN the application must immediately parse the configuration and detect the validation error before starting any computational work
  AND the command must fail, exiting with a non-zero status code to indicate an error
  AND the console must display a clear, user-friendly error message that is easy to understand
  AND this error message must specifically identify 'temperature_k' as the problematic parameter and explain the validation rule that it violated (e.g., "must be greater than zero")
  AND absolutely NO part of the computational pipeline (Generation, Exploration, etc.) should be executed
  AND consequently, NO output file named "NiAl_dataset.db" should be created in the user's directory.
```

**Scenario: Verifying Database Output Content**
```gherkin
GIVEN the "Successful End-to-End Pipeline Run" scenario has been completed successfully
  AND the output database file "CuAu_dataset.db" exists and is a valid SQLite file
WHEN the user inspects the contents of "CuAu_dataset.db" using an external tool like the `ase db` command
THEN the inspection tool must successfully connect to and read the database
  AND the tool must report that the database contains exactly 5 rows, which corresponds to the 5 atomic structures requested in the 'num_samples' configuration parameter
  AND for each row inspected, the chemical formula must only contain the elements 'Cu' and 'Au'
  AND the overall composition of the atoms in each structure must be approximately 50% Copper and 50% Gold, consistent with the 'composition' parameter defined in the configuration file.
```
