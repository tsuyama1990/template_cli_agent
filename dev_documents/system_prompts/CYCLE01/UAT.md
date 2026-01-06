# User Acceptance Testing (UAT): MLIP-AutoPipe - Cycle 1

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for the foundational first cycle of the MLIP-AutoPipe project. The primary objective of this UAT is to empower the end-user—typically a computational materials scientist or a researcher in a related field—to independently verify that the core functionalities delivered in this cycle are working correctly, are robust to common errors, and are user-friendly. The focus is squarely on validating the command-line interface (CLI), the robustness of the configuration system, and the integrity of the basic, end-to-end data pipeline. Success in this UAT is defined by the user's ability to confidently run the tool with a valid configuration, understand the error messages from an invalid one, and trust the contents of the generated output database. This builds the necessary trust and confidence in the tool's foundational architecture before more complex features are introduced in subsequent cycles.

To facilitate a smooth and informative testing experience, the UAT will be structured around a Jupyter Notebook (`UAT_CYCLE_01.ipynb`). This interactive format is deliberately chosen for its pedagogical value in a scientific context. It allows the user to not just execute commands, but to see the inputs, the outputs, and the verification logic side-by-side. The notebook will serve as both a tutorial and a test script, combining explanatory markdown text with executable code cells. This approach demystifies the testing process, allowing the user to actively participate in verifying the software's correctness in a transparent and reproducible manner. Each scenario is designed to test a specific aspect of the Cycle 1 deliverables, ensuring comprehensive coverage of the core features from the user's perspective.

| Scenario ID | Test Scenario | Priority |
| :--- | :--- | :--- |
| UAT-C1-01 | **Successful End-to-End Run with Valid Configuration** | **High** |
| UAT-C1-02 | **Graceful Failure with a Comprehensive Set of Invalid Configurations** | **High** |
| UAT-C1-03 | **Deep Verification of Database Output Integrity and Correctness** | **High** |

---

### **Scenario UAT-C1-01: Successful End-to-End Run with Valid Configuration** (Min 300 words)

**Description:** This scenario represents the primary "happy path" and is the most fundamental test of the application's basic functionality. Its purpose is to confirm that a user, provided with a correctly formatted configuration file, can successfully execute the entire (simplified) pipeline from start to finish without encountering any errors. A successful run validates several critical aspects of the system at once: that the command-line entry point is correctly installed and accessible in the user's environment, that the YAML configuration file is parsed and validated as expected by the Pydantic models, and that all the core components of the pipeline—the `AlloyGenerator`, the mock `LabelingEngine`, and the `AseDBWrapper`—are correctly instantiated and orchestrated by the `WorkflowOrchestrator`. It is, in essence, the first end-to-end test from a user's perspective, and its success is a prerequisite for any further testing. The user should feel that the tool is reliable and performs its basic advertised function as expected.

**Jupyter Notebook Steps:**
1.  **Environment Setup:** The notebook will begin with a setup cell. This cell will use shell commands (`!pip freeze`, `!which mlip-autopipec`) to confirm for the user that the `mlip_autopipec` package is installed in the current Python environment and that the command-line tool is accessible from the system's PATH. This step prevents confusion and ensures the user's environment is correctly configured before proceeding.
2.  **Configuration File Creation:** To ensure a reproducible test, the notebook will contain a Python code cell that programmatically generates a valid `config.yaml` file in the current directory. The user will be shown the Python code that defines the configuration dictionary and writes it to a file using `yaml.dump()`. The configuration will specify a simple, easily verifiable system: a Copper-Gold (Cu-Au) alloy with a 50/50 composition, and it will request the generation of exactly 10 structures. The output database path will be explicitly set to `cycle_01_test.db`. The content of the generated YAML file will then be printed to the screen, so the user can see exactly what is being used as input.
3.  **Pipeline Execution:** The notebook will have a cell containing the command to execute the main application: `!mlip-autopipec run --config config.yaml`. The `!` prefix allows the shell command to be run directly from the notebook. The standard output and standard error streams from the command will be captured and displayed in real-time within the notebook's output cell. This provides the user with direct, immediate feedback from the application.
4.  **Initial Verification:** After the command finishes, the user will be prompted to perform the first layer of verification. They should see that the application printed logging messages indicating a successful run and, most importantly, that the command completed with a zero exit code (which is the standard for success in command-line tools). The notebook will then contain a Python cell that uses `os.path.exists('cycle_01_test.db')` to programmatically check for the existence of the output database file, printing a clear, color-coded "PASS" or "FAIL" message. This confirms the most basic success criterion: the tool ran and produced the expected output file.

---

### **Scenario UAT-C1-02: Graceful Failure with a Comprehensive Set of Invalid Configurations** (Min 300 words)

**Description:** This scenario is arguably more important for building user trust than the happy path test. It rigorously evaluates the system's robustness and user-friendliness in the face of common configuration errors. A tool that crashes with an unhelpful traceback when given slightly incorrect input is frustrating to use. This test ensures that our Pydantic-based validation is working as intended, catching errors early and providing the user with clear, actionable feedback. The goal is for the user to feel that the tool is helpful even when they make a mistake. By testing several distinct types of invalid configurations—missing fields, incorrect data types, and out-of-bounds values—we demonstrate a comprehensive approach to input validation. The user should come away from this scenario confident that the tool will protect them from wasting time and computational resources on misconfigured runs.

**Jupyter Notebook Steps:**
1.  **Test 1: Missing Required Field:** The notebook will first create a `config_missing.yaml` file. In this file, a required field, such as `num_structures` in the `generation` section, will be deliberately omitted. The notebook will then execute the CLI with this file. The user will be directed to examine the output. The key verification is that the application exits with a non-zero status code and that the error message from Pydantic is displayed, which should clearly state: `Field required` for the missing key.
2.  **Test 2: Incorrect Data Type:** Next, a `config_wrong_type.yaml` file will be created. Here, `num_structures` will be set to a string value, for example, `"ten"`. When the CLI is run with this file, the user must verify that the error message clearly indicates a data type mismatch, for example, `Input should be a valid integer`. This confirms that the type hints in our Pydantic models are being correctly enforced.
3.  **Test 3: Out-of-Bounds Value:** A third file, `config_invalid_value.yaml`, will be created where `num_structures` is set to `0`. This value is of the correct type (integer) but violates the `gt=0` (greater than 0) constraint defined in the schema. The user will run the CLI and verify that the resulting error message is specific and helpful, stating something like `Input should be greater than 0`.
4.  **Test 4: Custom Validator Failure:** A final invalid file, `config_bad_composition.yaml`, will be created to test our custom validator. The `composition` dictionary will be set to `{'Cu': 0.6, 'Au': 0.5}`, which sums to 1.1. The user will run the CLI and verify that the error message is the one we defined in our custom validator: `Composition probabilities must sum to 1.0`.
5.  **Side-Effect Verification:** For each of the above tests, a final, crucial step will be included. The notebook will contain a cell that checks for the existence of the output database file and asserts that it does *not* exist. This confirms that the application is failing fast and is not creating partial or corrupt output files when it detects an invalid configuration.

---

### **Scenario UAT-C1-03: Deep Verification of Database Output Integrity and Correctness** (Min 300 words)

**Description:** This scenario completes the user's journey by focusing on the integrity and correctness of the final product: the output database. It is not sufficient for the tool to simply run without errors and create a file; the user must be able to trust that the data *within* that file is correct and structured as expected. This test provides the user with the tools and guidance to perform a deep inspection of the database contents, thereby verifying that the data has been correctly generated, processed through the (mock) pipeline, and written to storage without corruption or misinterpretation. This builds confidence in the entire data pipeline, from the `AlloyGenerator`'s logic to the `AseDBWrapper`'s serialization.

**Jupyter Notebook Steps:**
1.  **Prerequisite:** This scenario directly follows the successful run from UAT-C1-01. It will assume that the `cycle_01_test.db` file has been successfully created. The first cell will connect to this database using `ase.db.connect`.
2.  **Database Inspection and Verification:** The notebook will then guide the user through a series of programmatic checks using the ASE DB API. Each check will be in its own cell, with a clear explanation of what is being tested and a final `print("PASS")` or `print("FAIL")` statement.
    *   **Row Count Verification:** The first check will be on the total number of entries. The code `db.count()` will be executed, and the result will be compared against the `num_structures` (which was 10) from the `config.yaml`. The user will verify that the numbers match exactly.
    *   **Composition Verification:** To verify the generator's logic, the notebook will select the first row from the database using `db.get(id=1)`. This returns an `AtomsRow` object, which can be converted back into a full `ase.Atoms` object. The notebook will then get the list of chemical symbols from this object and calculate the composition. It will assert that the ratio of Cu to Au atoms is 50/50, as specified in the configuration. This provides direct proof that the structure generation logic is working correctly.
    *   **Data Presence and Correctness:** The final check is to ensure that the data from the mock labeling engine was correctly written. The test will access the `data` dictionary associated with the row (`row.data`). It will verify that the key `energy` exists and that its value is `0.0`, which was the dummy value supplied by the mock engine. It will also retrieve the `forces` numpy array and check that its shape is correct (N x 3, where N is the number of atoms in the structure). This confirms that the `DFTResult` Pydantic model was correctly serialized, passed through the orchestrator, and written to the database by the wrapper. The successful completion of this scenario gives the user high confidence in the integrity of the entire data pipeline.

## 2. Behavior Definitions

**Feature: Core CLI Pipeline Execution**

**Scenario: A user successfully runs the pipeline with a valid configuration file for a binary alloy.**
*   **GIVEN** a well-formed YAML configuration file named `config.yaml` is present in the current directory.
*   **AND** this file specifies the creation of `15` structures for a Silicon-Germanium ('Si', 'Ge') system with a `50/50` composition.
*   **AND** the output database path is explicitly defined in the configuration as `sige_test.db`.
*   **WHEN** the user executes the command `mlip-autopipec run --config config.yaml` in their terminal.
*   **THEN** the application should start execution and display informative log messages to the console.
*   **AND** the application should terminate gracefully with a clear success message and a zero exit code.
*   **AND** a new file named `sige_test.db` must be created in the current directory.
*   **AND** upon inspection, the `sige_test.db` file must be a valid SQLite database.
*   **AND** this database must contain a table with exactly `15` rows, corresponding to the requested number of structures.
*   **AND** each row in the database must represent an atomic structure containing both Silicon and Germanium atoms in approximately the correct ratio.

---

**Feature: Robust Configuration Validation**

**Scenario: A user attempts to run the pipeline with a configuration file that contains a clear validation error.**
*   **GIVEN** a YAML configuration file named `invalid_config.yaml` is present.
*   **AND** this file contains a logically invalid value, specifically `0`, for the `num_structures` parameter, which is constrained to be a positive integer.
*   **WHEN** the user executes the command `mlip-autopipec run --config invalid_config.yaml`.
*   **THEN** the application must not start the main processing pipeline.
*   **AND** it must terminate immediately with a non-zero exit code, indicating failure.
*   **AND** before exiting, it must display a clear, human-readable error message to the standard error stream.
*   **AND** this error message must specifically mention that the value for the `num_structures` field is invalid and explain the constraint it violates (e.g., "Input should be greater than 0").
*   **AND** no output database file should be created or modified on the filesystem.
