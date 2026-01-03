# User Acceptance Testing (UAT): Cycle 1

This document outlines the User Acceptance Testing (UAT) scenarios for the first development cycle of the MLIP-AutoPipe project. The exclusive focus of this cycle is the core command-line interface (CLI) pipeline. The primary goal of this UAT is to rigorously verify that an end-user, particularly a computational scientist, can successfully and reliably generate a complete dataset of atomic structures starting from a simple, human-readable configuration file. Successfully passing these tests will confirm that the foundational architecture is robust, the core components are correctly integrated, and the system is ready for the more advanced features planned for the next cycle. This UAT is designed not just as a verification step, but also as a practical, hands-on tutorial for a new user, demonstrating both the successful "happy path" and the application's robust error-handling capabilities. The use of Jupyter Notebooks as the testing tool is a deliberate choice to enhance the user experience, making the tests interactive, transparent, and easily reproducible.

## 1. Test Scenarios

### Scenario ID: UAT-C1-001
*   **Priority:** High
*   **Title:** Successful End-to-End Pipeline Run for a Binary Alloy
*   **Description:** This scenario represents the quintessential "happy path" and is the single most critical test for Cycle 1. Its purpose is to validate the entire pipeline's functionality from start to finish for a very common and representative use case: generating a dataset for a binary alloy. The user will begin by creating a simple, clean YAML configuration file that specifies the elements, composition, and basic simulation parameters. They will then execute the main CLI command, pointing to this configuration. The successful completion of this scenario will provide tangible proof that all the core architectural components—the Generator, the simplified Explorer, the Random Sampler, and the Storage writer—are correctly wired together and are functioning as a cohesive unit. This test is designed to be the first task a new user would perform, acting as a "hello, world" for the application. Its success will instill confidence in the tool's basic capabilities and provide a clear, working example from which users can build their own more complex configurations. It is the ultimate measure of whether the cycle's primary goal of delivering a functional core pipeline has been met.
*   **UAT Tool:** A carefully crafted Jupyter Notebook (`C1_UAT_Alloy.ipynb`) will be the primary vehicle for this test. This choice of tool is deliberate, as it allows for an interactive and transparent testing process. The notebook will be structured with clear markdown instructions and executable code cells that will:
    1.  Programmatically create the `config.yaml` file in the notebook's working directory. This ensures the test is self-contained and reproducible.
    2.  Use the `!mlip_autopipec` shell command syntax to execute the pipeline directly from within the notebook, capturing the command's output.
    3.  After the command completes, the notebook will use the `ase.db` library to connect to the resulting `final_dataset.db` file.
    4.  It will then perform a series of critical assertions: checking the number of rows in the database, verifying the chemical species of the atoms within the structures, and potentially plotting one of the structures for visual confirmation. This interactive approach is far more compelling and informative for the user than a simple pass/fail message from a script, as it allows them to directly inspect and interact with the final product, thus serving as both a test and a tutorial.

### Scenario ID: UAT-C1-002
*   **Priority:** Medium
*   **Title:** Configuration Validation and Graceful Error Handling
*   **Description:** This scenario is designed to test the robustness and user-friendliness of the system's configuration validation layer, which is handled by Pydantic. A powerful tool is not just one that works when given correct input, but one that behaves gracefully and informatively when given incorrect input. This is crucial for a positive user experience, as it prevents cryptic, hard-to-debug runtime errors and instead guides the user towards crafting a correct configuration. We will systematically test multiple common failure modes, such as providing incorrect data types (e.g., a string instead of a number for `md_steps`), missing required fields (e.g., omitting the `elements` list), and providing values that are syntactically correct but logically invalid (e.g., a negative number for `num_samples`, or composition fractions that do not sum to 1.0).
*   **UAT Tool:** A second, dedicated Jupyter Notebook (`C1_UAT_Validation.ipynb`) will be provided to demonstrate and verify the validation logic. This notebook will contain a series of tests, each targeting a specific type of invalid configuration. For each test case, the notebook will:
    1.  Clearly explain in markdown the error in the configuration.
    2.  Write the "bad" YAML content to a file.
    3.  Attempt to run the CLI command, configured to expect a failure.
    4.  The notebook will capture the standard error output stream from the failed CLI process.
    5.  It will then perform crucial assertions on the captured output, verifying not only that the command failed (i.e., exited with a non-zero status code), but also that the error message printed to the user is helpful and contains specific, expected keywords related to the Pydantic validation error (e.g., "field required", "value is not a valid integer", "Composition fractions must sum to 1.0"). This confirms that the validation is working and that the user is receiving actionable feedback.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN for UAT-C1-001

**GIVEN** a user has a clean working directory, signifying a fresh start.
**AND** the user has created a valid YAML configuration file named `config.yaml`, which is the primary input for the system. This file contains the following well-defined parameters:
```yaml
# A complete and valid configuration for a simple binary alloy.
system:
  elements: [Cu, Au]
  composition: {Cu: 0.5, Au: 0.5}
  generator_type: "alloy"
  num_initial_structures: 5

exploration:
  md_steps: 50
  temperature_k: 500.0

sampling:
  num_samples: 10

storage:
  db_path: "final_dataset.db"
```
**WHEN** the user executes the primary command-line tool from their terminal, providing the path to the configuration file, with the exact command: `mlip_autopipec --config config.yaml`.
**THEN** the pipeline process must execute and complete without any errors, ultimately finishing with a success exit code of 0, which is the standard indicator of a successful run in a command-line environment.
**AND** upon completion, a new file named `final_dataset.db` must be created in the user's current working directory. The presence of this file is the primary artifact of the run.
**AND** this `final_dataset.db` file must be a valid SQLite database file, readable by standard database tools and, more importantly, by the ASE library.
**AND** when a connection is made to this database using the `ase.db.connect` function, it must be found to contain exactly 10 rows. This number directly corresponds to the `sampling.num_samples` parameter in the configuration file, verifying that the sampling stage is respecting its configuration.
**AND** each row in the database, when read as an `ase.Atoms` object, must be a valid instance of this class.
**AND** a detailed inspection of the chemical symbols within each `Atoms` object must confirm that they contain only Copper (Cu) and Gold (Au) atoms, as specified in the `system.elements` configuration. This verifies that the entire pipeline is correctly propagating the system definition.

### GIVEN/WHEN/THEN for UAT-C1-002

**GIVEN** a user has a clean working directory.
**AND** the user has intentionally created an invalid YAML configuration file named `invalid_config.yaml`. This file contains a logical error where `num_initial_structures` is a negative number, which is physically nonsensical:
```yaml
# An invalid configuration file to test Pydantic validation.
system:
  elements: [Cu, Au]
  composition: {Cu: 0.5, Au: 0.5}
  generator_type: "alloy"
  num_initial_structures: -5 # This value is invalid as it must be a positive integer.

exploration:
  md_steps: 100
  temperature_k: 300

sampling:
  num_samples: 10

storage:
  db_path: "test.db"
```
**WHEN** the user executes the command-line tool, pointing it at this invalid configuration file with the command: `mlip_autopipec --config invalid_config.yaml`.
**THEN** the application must immediately fail to execute the pipeline and must exit with a non-zero status code, signaling that an error has occurred.
**AND** the application must print a clear and informative error message to the console (standard error).
**AND** this error message must be specific and actionable, clearly identifying the source of the problem. It must contain text explicitly stating that the `num_initial_structures` field failed its validation check because its value must be greater than zero (`gt=0`), directly guiding the user on how to fix their configuration file. This confirms the Pydantic validation layer is working correctly and providing a good user experience.
