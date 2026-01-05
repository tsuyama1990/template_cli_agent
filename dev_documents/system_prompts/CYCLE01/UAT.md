# User Acceptance Testing (UAT): CYCLE01

This document outlines the User Acceptance Testing (UAT) plan for the first development cycle of the MLIP-AutoPipe project. The focus of this UAT is to verify that the foundational features of the application—the command-line interface, the configuration system, and the basic data generation pipeline—work correctly from an end-user's perspective and provide a satisfactory, intuitive, and robust user experience. The tests are designed not just to check for correctness, but to ensure the tool is practical and trustworthy for its intended audience of computational scientists. The scenarios are presented in a way that mimics a real user's workflow, starting from a simple command and progressively verifying the integrity and usability of the output. A key aspect of this UAT is the use of Jupyter notebooks, which is a very common tool in the scientific community. By demonstrating seamless integration with this tool, we can significantly increase user confidence and adoption. This UAT plan serves as the final quality gate before the foundational components of Cycle 1 are considered complete and ready to be built upon in the next cycle.

## 1. Test Scenarios

The following scenarios are designed to be executed by a user to confirm that the core functionality delivered in Cycle 1 meets the project's requirements. They cover the "happy path" as well as common user error conditions, ensuring the application is not only functional but also resilient and user-friendly.

| Scenario ID | Description                                                                 | Priority |
| :---------- | :-------------------------------------------------------------------------- | :------- |
| UAT-C1-001  | Successful End-to-End Run with a Valid Configuration                        | High     |
| UAT-C1-002  | Application Handles a Missing Configuration File Gracefully               | High     |
| UAT-C1-003  | Application Rejects a Configuration with Invalid Schema (Bad Composition) | High     |
| UAT-C1-004  | Verify the Contents and Integrity of the Output Database                    | High     |
| UAT-C1-005  | Interactive Verification and Analysis using a Jupyter Notebook              | Medium   |

---

### **Scenario UAT-C1-001: Successful End-to-End Run with a Valid Configuration**

*   **Objective:** To verify that the user can successfully run the entire pipeline from the command line using a standard, valid configuration file and that the process completes without any errors or unexpected warnings. This is the primary "happy path" test and represents the most common use case for the application. It is designed to build the user's initial confidence in the tool's reliability and ease of use. A successful run should be a smooth, informative experience.
*   **User Story:** As a computational materials scientist, I want to be able to start a data generation run for a simple binary alloy using a single, unambiguous command. I expect the tool to run to completion and provide me with a database, so that I can quickly and easily generate a basic dataset for my initial MLIP training experiments.
*   **Preconditions:**
    *   The MLIP-AutoPipe CLI application has been correctly installed and is accessible in the shell environment (i.e., it is in the system's PATH).
    *   A valid YAML configuration file named `config_valid.yaml` has been created in the user's working directory. This file is the sole input to the application.
*   **Test Steps:**
    1.  Create a file named `config_valid.yaml` with the following content. This configuration defines a simple, representative task.
        ```yaml
        project_name: uat_test_fept
        system:
          elements: ['Fe', 'Pt']
          composition: {'Fe': 0.5, 'Pt': 0.5}
          lattice: 'fcc'
          num_structures: 10
        exploration:
          temperature: 300.0
        sampling:
          method: 'random'
          fraction: 0.8
        ```
    2.  Open a terminal or command prompt.
    3.  Navigate to the directory containing the `config_valid.yaml` file.
    4.  Execute the command: `mlip-autopipec run --config config_valid.yaml`
*   **Expected Result:**
    *   The command executes immediately and finishes without raising any unhandled Python exceptions or printing any stack traces. The entire process should feel clean and professional.
    *   The console displays a series of user-friendly status messages, ideally with color and formatting (from the Rich library), indicating the progress of the pipeline (e.g., "[1/4] Starting structure generation...", "[2/4] Executing exploration stage...", "[3/4] Sampling structures...", "[4/4] Writing to database...", "Pipeline complete.").
    *   A new file named `uat_test_fept.db`, corresponding to the `project_name` in the configuration, is created in the current working directory.
    *   The process exits with a status code of 0, which is the standard convention for a successful command-line application execution.

---

### **Scenario UAT-C1-003: Application Rejects a Configuration with Invalid Schema (Bad Composition)**

*   **Objective:** To verify that the application's schema validation is working correctly and provides clear, actionable feedback to the user when the configuration file contains logical errors. This test is crucial for ensuring a good user experience, as it prevents the user from wasting time running a flawed calculation and helps them to fix their input easily.
*   **User Story:** As a user, if I make a mistake in my configuration file, such as providing element compositions that do not sum to 100%, I want the program to fail fast and tell me exactly what is wrong, so that I can correct my input without having to guess the cause of the error.
*   **Preconditions:**
    *   The MLIP-AutoPipe CLI application is installed.
    *   A configuration file named `config_invalid_comp.yaml` has been created with a clear logical error.
*   **Test Steps:**
    1.  Create a file named `config_invalid_comp.yaml` with a composition that sums to 1.1 instead of 1.0.
        ```yaml
        project_name: uat_test_invalid
        system:
          elements: ['Fe', 'Pt']
          composition: {'Fe': 0.6, 'Pt': 0.5} # Invalid: sums to 1.1
          lattice: 'fcc'
          num_structures: 10
        ...
        ```
    2.  In a terminal, execute the command: `mlip-autopipec run --config config_invalid_comp.yaml`
*   **Expected Result:**
    *   The application should exit almost immediately, within a second or two. It should not attempt to start any of the computationally intensive pipeline stages.
    *   A clear, user-friendly error message should be printed to the console. The message must explicitly state that the validation failed because the compositions do not sum to 1.0. It should point out the specific fields that are incorrect. For example: "Error: Configuration validation failed. The sum of compositions in `system.composition` must be 1.0, but it is 1.1."
    *   The application should exit with a non-zero status code (e.g., 1), indicating that an error occurred.
    *   No output database file (`uat_test_invalid.db`) should be created.

---

### **Scenario UAT-C1-005: Interactive Verification and Analysis using a Jupyter Notebook**

*   **Objective:** To provide a highly user-friendly and practical way to interact with and verify the generated data, using a tool that is ubiquitous in the scientific Python ecosystem. This test demonstrates the usability and interoperability of the output format with standard scientific libraries, which is a key requirement for user adoption. It serves as a tutorial for how a user would typically begin their own analysis.
*   **User Story:** As a scientist, I want to be able to easily load and analyse the generated structures in a Jupyter Notebook, which is my primary environment for research and analysis. I need to be able to inspect the structures, verify their properties, and integrate the data into my existing analysis and visualization workflows seamlessly.
*   **Preconditions:**
    *   Scenario UAT-C1-001 has been executed successfully.
    *   The output database file `uat_test_fept.db` exists.
    *   Jupyter Notebook or JupyterLab is installed in the user's Python environment, along with the ASE library.
*   **Test Steps:**
    1.  Create a new Jupyter Notebook named `verify_data.ipynb`.
    2.  In the first cell, add and execute the Python code to import the necessary libraries (ASE and NumPy) and to connect to the database. This step verifies that the database file is a standard, accessible SQLite file.
        ```python
        from ase.db import connect
        import numpy as np

        # Connect to the database generated by the pipeline
        db = connect('uat_test_fept.db')
        print(f"Successfully connected to the database.")
        print(f"Database contains {len(db)} structures.")
        ```
    3.  In a second cell, add and execute code to load the first atomic structure from the database and perform a series of critical checks on its physical properties. This verifies the scientific integrity of the generated data.
        ```python
        # Retrieve the first structure as an ASE Atoms object
        atoms = db.get_atoms(id=1)
        print("Successfully loaded structure with ID 1.")
        print("Chemical symbols:", atoms.get_chemical_symbols())

        # Check for a reasonable, non-zero cell volume
        assert atoms.get_volume() > 1.0, "Cell volume is zero or negative!"

        # Check for the minimum distance between any two atoms to ensure no overlaps
        # This is a critical check for physical validity
        distances = atoms.get_all_distances(mic=True) # mic=True for periodic boundary conditions
        min_dist = np.min(distances[np.nonzero(distances)])
        print(f"Minimum inter-atomic distance: {min_dist:.2f} Angstrom")
        assert min_dist > 1.0, "Atoms are overlapping, structure is not physically valid!"

        print("\nVerification successful: The structure is valid and has the correct properties.")
        ```
*   **Expected Result:**
    *   The first cell executes without errors and prints the message: "Successfully connected to the database. Database contains 8 structures."
    *   The second cell also executes without any errors or failed assertions. The chemical symbols will be a list containing the correct ratio of 'Fe' and 'Pt'.
    *   The printed minimum inter-atomic distance is a positive number greater than 1.0 Angstrom, which provides strong confirmation of the physical validity of the generated structure. The test demonstrates that the user can immediately start working with the data using standard tools.

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behaviour of the system in response to user actions in a formal, structured way. This provides an unambiguous specification of the system's requirements from a user's point of view, which can be used to guide both development and testing.

**Feature: Core Data Generation Pipeline**

As a user of the MLIP-AutoPipe tool, I want to be able to run the data generation pipeline via a simple command-line interface so that I can efficiently and reliably create a dataset for training my machine learning interatomic potential models. The process should be transparent, providing clear feedback, and the output should be in a standard, usable format.

**Scenario: Running the pipeline successfully with a valid configuration file**

*   **GIVEN** I have created a valid configuration file named "config_valid.yaml".
*   **AND** this file specifies the generation of 10 initial structures for an FePt alloy with a 50/50 composition.
*   **AND** the configuration specifies a random sampling method with a fraction of 0.8.
*   **WHEN** I execute the command `mlip-autopipec run --config config_valid.yaml` in my terminal.
*   **THEN** the application should start and run to completion without any errors or stack traces being displayed.
*   **AND** it should produce a final database file named "uat_test_fept.db", as specified by the `project_name` in the configuration.
*   **AND** this database should contain exactly 8 atomic structures (10 initial structures * 0.8 sampling fraction).
*   **AND** upon inspection, each structure in the database should be composed of Iron and Platinum atoms in a 50/50 ratio.
*   **AND** each structure should be physically plausible, with no overlapping atoms.

---

**Feature: Robust User Input Validation**

As a user, I need the application to provide immediate and clear feedback when my configuration is incorrect, so that I can easily identify and fix the errors without needing to debug the application's internal state. This is essential for a good user experience and for preventing wasted time and computational resources.

**Scenario: Attempting to run the pipeline with a path to a non-existent configuration file**

*   **GIVEN** I am in a directory where there is no file named "non_existent_config.yaml".
*   **WHEN** I execute the command `mlip-autopipec run --config non_existent_config.yaml`.
*   **THEN** the application should exit gracefully within a few seconds.
*   **AND** it must display a clear, human-readable error message in the console, for example: "Error: The specified configuration file 'non_existent_config.yaml' was not found."
*   **AND** it must not create any new files, such as an empty database file.
*   **AND** the command's exit code should be non-zero to indicate failure.

**Scenario: Attempting to run the pipeline with a configuration file containing a logical validation error**

*   **GIVEN** I have a configuration file where the chemical compositions are logically flawed (e.g., they sum to 1.1 instead of 1.0).
*   **WHEN** I execute the `run` command pointing to this invalid file.
*   **THEN** the application must immediately reject the configuration before starting any computational steps.
*   **AND** it must display a clear error message that precisely identifies the source of the validation error, for instance: "Configuration Error: The sum of compositions in `system.composition` must be 1.0."
*   **AND** it must not proceed to the structure generation stage or any subsequent stages of the pipeline.
*   **AND** it must not create any database file.
*   **AND** the command's exit code should be non-zero.
