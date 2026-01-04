# User Acceptance Testing (UAT): Cycle 1 - Core CLI Pipeline

This document outlines the User Acceptance Testing (UAT) plan for the first development cycle of the MLIP-AutoPipe project. The focus of this cycle is squarely on the core command-line interface (CLI) functionality, which forms the foundational engine of the entire application. The primary goal of this UAT is to provide a structured process for a target user (e.g., a computational materials scientist) to verify that the CLI tool can successfully generate a complete, valid, and scientifically sound dataset for MLIP training. The tests are designed to confirm that the pipeline operates reliably, handles user errors gracefully, and produces outputs that meet the project's core objectives of physical realism and automation. The UAT will be conducted by following a series of realistic test scenarios that cover the main use-cases of the Cycle 1 functionality. A key tool for this process will be a provided Jupyter Notebook (`UAT_Cycle1.ipynb`), which will contain pre-written code snippets to help the user execute the tests, analyze the outputs, and visualize the generated atomic structures. This ensures that the validation process is not only thorough but also accessible and educational for the user.

## 1. Test Scenarios

The following scenarios are designed to be executed by a user to validate the functionality and robustness of the core CLI pipeline. They are structured to cover the most critical "golden path" workflow, as well as common error conditions and important features like checkpointing. Each scenario is given a priority to indicate its importance to the overall success of the cycle.

| Scenario ID | Priority | Description |
| :--- | :--- | :--- |
| **UAT-C1-001** | High | **Successful End-to-End Run for a Binary Alloy**: This is the most critical test. It verifies that a user can successfully execute a complete pipeline run for a standard binary alloy system (e.g., Silicon-Germanium) from start to finish. The test confirms the generation of a valid ASE database containing a diverse set of structures, proving that all four pipeline stages are correctly integrated and functional. |
| **UAT-C1-002** | High | **Robust Handling of Invalid User Input**: This scenario tests the system's resilience and user-friendliness by verifying that it provides clear, understandable, and actionable error messages when supplied with deliberately invalid configuration files. This is crucial for ensuring a positive user experience and preventing wasted computational time due to simple configuration mistakes. |
| **UAT-C1-003** | Medium | **Successful Restart from a Checkpoint**: This test validates the pipeline's ability to correctly resume from an intermediate state. The user will simulate an interrupted workflow by running only the generation stage, and then restarting the full pipeline. The test will verify that the second run correctly uses the previously generated structures, saving computational resources and demonstrating robustness. |

---

### **Scenario UAT-C1-001: Successful End-to-End Run (Alloy)**

**(Min 300 words)**
This scenario represents the primary "golden path" for the application and is the most critical validation of its core functionality. It is designed to simulate a typical, real-world use case from the perspective of a computational scientist. The user will be tasked with generating a small but complete training dataset for a Silicon-Germanium (SiGe) alloy, a widely studied semiconductor material. The process will begin with a well-documented YAML configuration file, `sige_config.yml`. This file will specify a 50/50 Si/Ge composition on an FCC lattice, and will configure the pipeline to generate 4 initial seed structures. The exploration stage will be configured for a short, 100-step Molecular Dynamics simulation running in the NVT ensemble at a temperature of 300K. For speed and simplicity in this test, a fast classical potential like ASE's built-in EMT will be used instead of a full MLIP. The sampling stage will be set to use the `random` method to select 10 final structures for the dataset.

The user will execute the pipeline by running a single command in their terminal: `mlip-autopipec run --config sige_config.yml`. The primary expected outcome is a successful execution with no errors, which will be indicated by a clean exit to the command line (exit code 0) and a final "Pipeline completed successfully" message printed to the console. Following the run, the user will perform a detailed verification of the outputs. They will inspect the generated output directory and confirm the presence of the final ASE database file (`sige_dataset.db`), as well as intermediate files like `initial_structures.xyz` and the trajectory files. Using the provided `UAT_Cycle1.ipynb` Jupyter Notebook, the user will connect to the database and programmatically verify that it contains exactly 10 atomic structures, as requested in the configuration. They will then extract several of these structures at random and use the notebook's visualization tools to render them in 3D, visually confirming that they are physically plausible (e.g., atoms have reasonable bond lengths and are not overlapping). As a final scientific sanity check, they will plot a histogram of the potential energies of all structures in the database to confirm that the MD simulation has indeed explored a range of energies, not just a single static configuration. Successful completion of all these verification steps will confirm that every stage of the pipeline—Generation, Exploration, Sampling, and Storage—has worked together correctly to produce a valid and useful result.

### **Scenario UAT-C1-002: Input Validation Failure**

**(Min 300 words)**
A critical feature of user-friendly scientific software is its ability to handle user error gracefully and provide clear, actionable feedback. This scenario is designed to rigorously test the system's input validation capabilities, which are primarily implemented using Pydantic models. The goal is to ensure that the application fails fast when given bad input, preventing the user from launching a long, computationally expensive run that is doomed to fail due to a simple typo. The user will be provided with a series of deliberately broken configuration files and will attempt to run the pipeline with each.

The first test case will be `invalid_composition.yml`, a configuration file where the chemical `composition` fractions do not sum to 1.0 (e.g., `Si: 0.5`, `Ge: 0.4`). When the user runs the CLI with this file, the application is expected to exit immediately with a non-zero status code and print a clear, human-readable error message to the console, explicitly stating that the composition fractions must sum to 1.0. The second test case, `bad_type.yml`, will involve a type error, such as providing a string value for a parameter that should be an integer (e.g., `md_steps: "one_hundred"`). The expected outcome is a Pydantic `ValidationError` that is caught and presented to the user in a clean format, indicating the exact field (`md_steps`), the expected type (`integer`), and the incorrect value that was provided. A third test case, `unknown_field.yml`, will include a parameter that does not exist in the schema (e.g., `temperature_celcius: 100`). Because the Pydantic models are configured with `extra="forbid"`, the application should fail immediately and report that `temperature_celcius` is not a valid configuration key. This UAT is considered passed if, for every invalid configuration file tested, the system refuses to run the pipeline and instead outputs a precise and helpful error message that directly guides the user on how to fix their configuration.

## 2. Behavior Definitions

The following Gherkin-style definitions provide a structured, formal description of the expected behavior for each test scenario. This format ensures that the requirements are unambiguous and can be used as a clear reference for verification.

---

**Scenario: UAT-C1-001 - Successful End-to-End Run (Alloy)**

```gherkin
GIVEN a valid YAML configuration file named "sige_config.yml" is present in the user's working directory.
AND this file specifies a system configuration for a Si-Ge alloy with a 50/50 composition.
AND the file configures the generation of 4 initial structures.
AND the file configures a 100-step MD exploration in the NVT ensemble at 300K.
AND the file configures the sampling stage to use the 'random' method to select 10 final structures.
AND the file specifies the output database name as "sige_dataset.db".

WHEN the user executes the following command in their terminal:
  `mlip-autopipec run --config sige_config.yml`

THEN the command should execute without any errors being printed to the console.
AND the command should terminate with a final success message and a status code of 0.
AND an ASE database file named "sige_dataset.db" must be created in the output directory.
AND this database file must contain exactly 10 atomic structures.
AND a log file must be created, and it should contain log entries indicating the successful start and completion of all four pipeline stages (Generate, Explore, Sample, Store).
AND upon inspection with the provided Jupyter Notebook, the structures within the database must be physically plausible and show a distribution of energies.
```

---

**Scenario: UAT-C1-002 - Input Validation Failure**

```gherkin
GIVEN a YAML configuration file named "invalid_composition.yml" is present.
AND in this file, the composition fractions for the elements are defined as `Si: 0.5` and `Ge: 0.4`, which sum to 0.9.

WHEN the user executes the command:
  `mlip-autopipec run --config invalid_composition.yml`

THEN the command must fail immediately and exit with a non-zero status code.
AND the application must print a clear error message to the console.
AND this error message must contain text explicitly stating that the composition values are invalid because they do not sum to 1.0.

# ---

GIVEN a YAML configuration file named "bad_type.yml" is present.
AND in this file, the `md_steps` parameter is incorrectly specified as a string: `md_steps: "500"`.

WHEN the user executes the command:
  `mlip-autopipec run --config bad_type.yml`

THEN the command must fail immediately and exit with a non-zero status code.
AND the application must print a clear validation error message.
AND this error message must indicate that the field `md_steps` expects an integer but received a string.
```

---

**Scenario: UAT-C1-003 - Restart from Checkpoint**

```gherkin
GIVEN a valid YAML configuration file named "checkpoint_config.yml" is present.
AND the user has previously initiated a run that was interrupted after the Generation stage.
AND as a result, an intermediate file named `initial_structures.xyz` exists in the output directory.

WHEN the user executes the command:
  `mlip-autopipec run --config checkpoint_config.yml`

THEN the application should detect the presence of the existing `initial_structures.xyz` file.
AND the application's log output should contain a clear message indicating that it is skipping the Generation stage and is using the pre-existing structures.
AND the Exploration stage should commence immediately, using the structures loaded from the `initial_structures.xyz` file.
AND the rest of the pipeline (Explore, Sample, Store) should complete successfully, producing the final database as expected.
```
