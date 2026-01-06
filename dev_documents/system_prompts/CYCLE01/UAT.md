# User Acceptance Testing (UAT): MLIP-AutoPipe Cycle 1

This document outlines the User Acceptance Testing (UAT) scenarios for the first cycle of the MLIP-AutoPipe project. The goal of this UAT is to rigorously verify that the core command-line functionality is robust, user-friendly, and correctly generates a physically plausible dataset for a simple, representative use case. This UAT is designed to be performed using an interactive Jupyter Notebook, which not only serves as a transparent and reproducible testing artifact but also doubles as a tutorial for new users, providing a clear, step-by-step demonstration of the tool's capabilities and the quality of its output.

## 1. Test Scenarios

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-C1-01   | Successful End-to-End Pipeline Execution    | High     |
| UAT-C1-02   | Handling of Invalid Configuration File      | High     |
| UAT-C1-03   | Verification of Output Database Integrity   | Medium   |

### UAT-C1-01: Successful End-to-End Pipeline Execution

**(Min 300 words)**
This scenario represents the primary "happy path" and is the most critical validation of the core functionality for Cycle 1. The user persona is a computational materials science researcher who needs to quickly generate a preliminary training dataset for a copper-gold (CuAu) alloy, a classic benchmark system. The user begins by creating a simple, clean YAML configuration file. This file will specify the fundamental parameters: elements as `['Cu', 'Au']`, composition as `{'Cu': 0.5, 'Au': 0.5}`, a modest 3x3x3 supercell to keep the atom count manageable, and parameters for a brief but meaningful molecular dynamics run (e.g., 500 steps at 500K). The core expectation is that the user can execute a single, straightforward command, `mlip-autopipec run --config-file config.yml`, from their terminal and have the entire data generation process run to completion without any errors or need for further intervention. The tool must provide clear, real-time feedback to the command line, indicating the start, progress, and successful completion of each of the four distinct pipeline stages: Generation, Exploration, Sampling, and Storage. This verbose feedback is crucial for user trust, as it makes the process transparent rather than a "black box." Upon completion, the user anticipates finding a new SQLite database file, precisely named as specified in the config (e.g., `CuAu_dataset.db`), in their working directory. This test is of paramount importance as it validates the entire core workflow, confirming that all the individual service components (Generator, Explorer, Sampler, Storage) are correctly integrated and functioning as a cohesive unit. A successful run will be a powerful demonstration that the tool can reliably translate a simple, high-level declarative request into a concrete, useful, and structured data artifact, thereby fulfilling the central promise of the application. The user experience should be smooth, predictable, and confidence-inspiring, encouraging the user to trust the tool for larger and more complex tasks in the future. Furthermore, the output should be deterministic and reproducible; running the pipeline again with the same configuration and a fixed random seed should produce an identical database, a cornerstone of scientific computing.

### UAT-C1-02: Handling of Invalid Configuration File

**(Min 300 words)**
This scenario is designed to test the tool's robustness, intelligence, and user-friendliness when confronted with common and expected user errors in the input. The user persona is a scientist who is new to the tool and is likely to make mistakes while authoring their first configuration file. We will test several distinct error cases. In one case, the user specifies a negative temperature (`temperature_k: -100`) for the molecular dynamics simulation, a physically nonsensical value. In another case, they define a composition for their alloy that does not sum to 1.0 (e.g., `{'Cu': 0.6, 'Au': 0.5}`). A third case might involve a simple typo, using an invalid chemical symbol (e.g., `elements: ['Cu', 'Az']`). The non-negotiable expected behavior is that the tool should fail gracefully and, most importantly, immediately upon parsing the invalid configuration. It must not proceed to any computationally expensive tasks like structure generation or launching parallel MD simulations. Instead, the application must provide a clear, concise, and genuinely helpful error message to the user directly on the command line. The message must be specific, stating exactly which configuration parameter is invalid and precisely why (e.g., "Validation Error in ExplorationConfig: `temperature_k` must be a positive number," or "Validation Error in SystemConfig: Composition values must sum to 1.0," or "Validation Error: 'Az' is not a valid chemical symbol."). This immediate and informative feedback is a critical feature for a positive user experience, as it empowers the user to quickly identify and rectify the mistake in their configuration file without the frustrating process of digging through log files or consulting extensive documentation. This test rigorously verifies that the Pydantic-based validation layer is working correctly and is properly integrated into the CLI entry point. A successful outcome will demonstrate that the tool is robust against common user errors and actively guides the user towards providing valid input, thereby preventing wasted time, computational resources, and user frustration.

### UAT-C1-03: Verification of Output Database Integrity

**(Min 300 words)**
This scenario focuses on the scientific correctness and integrity of the final output, which is the most important measure of the tool's value. Following a successful pipeline execution from scenario UAT-C1-01, the user's goal is to meticulously inspect the generated database to ensure the data is physically valid and strictly conforms to all input specifications. The primary tool for this verification will be a Jupyter Notebook, provided as part of the UAT materials. This notebook will contain boilerplate Python code using the ASE and pandas libraries to connect to the output SQLite database file and perform a series of analytical checks. The user will execute these checks step-by-step. First, they will query the total number of entries (rows) in the database and verify that this number exactly matches the expected count based on the configuration (`num_initial_structures` * `num_samples`). Second, they will extract several random `Atoms` objects from the database and perform detailed assertions on their properties. They will confirm that the chemical species present in the structures are exclusively Cu and Au, and that the overall composition of each structure is approximately the 50-50 mix specified in the input configuration. Third, and most critically, the user will perform a physical validity check by calculating all interatomic distances within a sample structure, ensuring that no two atoms are closer than a reasonable physical bond length (e.g., ~1.5 Angstroms). This confirms that the overlap-checking logic in the generation phase is working correctly. The notebook will also include simple 3D visualizations of a few of the structures using `nglview`, providing a qualitative "eyeball" check that they look like reasonable, non-pathological crystal structures. This scenario is absolutely vital for building user trust in the scientific validity of the tool's output.

## 2. Behavior Definitions

**(Min 500 words)**

### Gherkin-style Definitions

**Scenario: UAT-C1-01 - Successful End-to-End Pipeline Execution**

*   **GIVEN** a user has created a valid YAML configuration file named `config.yml`.
*   **AND** the file specifies a system with `elements: ['Cu', 'Au']` and `composition: {'Cu': 0.5, 'Au': 0.5}`.
*   **AND** the file specifies the creation of 5 initial structures and the sampling of 10 frames from each.
*   **AND** the file specifies an output database path of `CuAu_dataset.db`.
*   **AND** the user is in a terminal session in the same directory as the `config.yml` file.
*   **WHEN** the user executes the command `mlip-autopipec run --config-file config.yml`.
*   **THEN** the application should print a clear message to the console indicating the start of the "Generation" stage.
*   **AND** the application should subsequently print messages indicating the successful completion of "Generation" and the start of "Exploration".
*   **AND** this pattern of start and completion messages should continue for the "Sampling" and "Storage" stages.
*   **AND** the application should finally print a message like "Pipeline executed successfully."
*   **AND** the application should exit with a success code (0).
*   **AND** a file named `CuAu_dataset.db` must exist in the current directory.
*   **AND** a connection to `CuAu_dataset.db` should reveal exactly 50 rows of data.

**Scenario: UAT-C1-02 - Handling of Invalid Configuration File**

*   **GIVEN** a user has a configuration file named `config_invalid.yml` where the `temperature_k` key under `exploration` is set to `-100`.
*   **AND** the user is in a terminal in the same directory as the `config_invalid.yml` file.
*   **WHEN** the user executes the command `mlip-autopipec run --config-file config_invalid.yml`.
*   **THEN** the application should immediately exit with a non-zero error code.
*   **AND** the application must print a user-friendly error message to the standard error stream.
*   **AND** the error message must contain the name of the invalid field (`temperature_k`) and a human-readable explanation of the validation rule it violated (e.g., "must be positive").
*   **AND** no output database file should be created.
*   **AND** no messages indicating the start of any pipeline stage (like "Generation" or "Exploration") should be printed.
*   **GIVEN** a user has another configuration file `config_invalid_comp.yml` where the `composition` values sum to 1.1.
*   **WHEN** the user executes the command `mlip-autopipec run --config-file config_invalid_comp.yml`.
*   **THEN** the application should immediately exit with a non-zero error code.
*   **AND** the error message must clearly indicate that the composition values must sum to 1.0.

**Scenario: UAT-C1-03 - Verification of Output Database Integrity**

*   **GIVEN** a user has successfully executed the pipeline, which has generated a `CuAu_dataset.db` file from a valid configuration.
*   **AND** the configuration specified 10 initial structures and to sample 5 frames from each trajectory.
*   **AND** the user has launched a Jupyter Notebook environment with the ASE library available.
*   **WHEN** the user executes a notebook cell containing the Python code `from ase.db import connect; db = connect('CuAu_dataset.db')`.
*   **THEN** the database connection object `db` should be created successfully without raising an exception.
*   **AND** when the user executes the cell `len(db)`, the output must be exactly 50.
*   **AND** when the user executes a cell to retrieve a random row (`row = db.get(3)`) and inspect its atomic numbers (`set(row.numbers)`), the resulting set must be `{29, 79}`.
*   **AND** when the user executes a cell to calculate the composition of that row, the result should be approximately 0.5 for both elements.
*   **AND** when the user executes a cell to get all interatomic distances (`row.get_all_distances()`) and finds the minimum of the non-zero distances, the result must be greater than a physically realistic threshold of 1.5 Angstroms.
*   **AND** when the user visualizes the structure, it should appear as a plausible, non-fragmented crystal lattice.
