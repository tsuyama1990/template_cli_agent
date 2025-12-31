# UAT.md: Cycle 02 - Initial Structure Generation & Configuration

## 1. Test Scenarios

User Acceptance Testing for Cycle 02 is critically focused on the new user-facing entry point of the system and the automation it promises. The key goal of this UAT is to rigorously verify that a non-expert user can successfully and confidently initiate a complex, multi-stage workflow by providing only a minimal, human-readable configuration file. This cycle's UAT is designed to test the system's ability to intelligently and robustly automate the entire setup and initial data generation phases, providing tangible confirmation that the core project philosophy of "removing the human expert from the loop" is being successfully implemented in practice. The user must be assured that the system is making sensible, physically-grounded decisions on their behalf.

| Scenario ID | Scenario Description                                       | Priority |
| :---------- | :--------------------------------------------------------- | :------- |
| UAT-C02-01  | **Successful and Verifiable Workflow Initiation from a Minimal `input.yaml`** | High     |
|             | This scenario represents the primary "happy path" and the most common user interaction for this cycle. The user's goal is to confirm that they can start the entire, complex pipeline with the simplest possible input. They will create a very basic `input.yaml` file specifying only the chemical elements of their system. They will then execute the main command-line tool. The user's expectation is that the system will automatically perform two major tasks: first, generate a full, detailed, and explicit configuration file for inspection and provenance; and second, produce a set of physically plausible initial atomic structures in the database, ready for the subsequent DFT labelling phase. Successfully passing this test validates the entire initialization chain from user input to database population. |          |
| UAT-C02-02  | **Correctness and Physical Sensibility of Heuristic-Based Parameter Generation**             | High     |
|             | This scenario is designed to build the user's trust in the system's "intelligence" by verifying the quality of the decisions made by the Config Expander. The user will provide a minimal config for a well-known and non-trivial material, such as the magnetic alloy FePt. They will then meticulously inspect the system-generated `exec_config_dump.yaml`. The user expects to see not just valid, but physically sensible and well-chosen parameters that an expert would select. This includes the correct identification of the material's structural type ('alloy'), the selection of appropriate high-precision DFT settings from a standard protocol (like SSSP), and the enabling of relevant physics (like magnetism). This test confirms that the embedded heuristics are genuinely useful and reliable. |          |
| UAT-C02-03  | **Robust Validation of User Input and Graceful, Informative Error Handling**        | Medium   |
|             | This scenario tests the system's robustness and user-friendliness when confronted with incorrect input. The user will create a deliberately invalid `input.yaml` file, for example, with a typo in a required key name (`elemants` instead of `elements`), or by specifying a chemically nonsensical system. The user expects the system to "fail fast" and provide a clear, user-friendly error message that helps them diagnose and fix the problem. The system should not proceed with incorrect data, nor should it crash with an unreadable technical stack trace. This UAT is crucial for ensuring a positive user experience, especially for beginners, and for preventing the waste of valuable computational resources on misconfigured runs. |          |

## 2. Behavior Definitions

These behaviors will be tested by a user interacting with the system solely through the command line, as intended for the final product. The verification will involve executing the main `mlip-pipe` command and then inspecting the generated files (`exec_config_dump.yaml`, `mlip.db`) and console log output to confirm the expected outcomes.

### Scenario: UAT-C02-01 - Successful and Verifiable Workflow Initiation from a Minimal `input.yaml`

*   **GIVEN** a clean, newly created project directory that contains no existing database, configuration files, or other artifacts.
*   **AND** the user has used a simple text editor to create a new file named `input.yaml` containing the absolute minimal required information for a run:
    ```yaml
    system:
      elements: ["Si"]
      composition: "Si"
    ```
*   **AND** the `mlip-pipe` command-line tool has been correctly installed and is accessible in the user's shell environment.
*   **WHEN** the user executes the primary command `uv run mlip-pipe input.yaml` from their terminal within the project directory.
*   **THEN** the command must execute to completion without raising any unhandled exceptions and must terminate with a clear success message printed to the log.
*   **AND** the user must be able to verify that a new file, named `exec_config_dump.yaml`, has been created in the current directory. This file should be human-readable and contain a comprehensive list of all execution parameters.
*   **AND** the user must also verify that a new SQLite database file, named `mlip.db`, has been created.
*   **AND** upon connecting to this new database with a SQLite client, the user must find that it is populated with multiple, distinct atomic structures of Silicon. Each of these entries must be correctly marked with an initial state, such as 'initial_generated', to clearly indicate their origin and readiness for the next workflow stage.

### Scenario: UAT-C02-02 - Correctness and Physical Sensibility of Heuristic-Based Parameter Generation

*   **GIVEN** a clean project directory.
*   **AND** the user has created a slightly more complex `input.yaml` file for a common binary alloy that has important physical properties like magnetism:
    ```yaml
    system:
      elements: ["Fe", "Pt"]
      composition: "FePt"
    ```
*   **WHEN** the user executes the command `uv run mlip-pipe input.yaml`.
*   **THEN** a detailed `exec_config_dump.yaml` file must be generated.
*   **AND** the user will open this `exec_config_dump.yaml` file and meticulously verify the correctness and physical sensibility of the following key auto-generated settings:
    *   Under the `system` section, the `structure_type` key must have been correctly inferred and set to the string `alloy`.
    *   Under the `dft_compute` section, the `pseudopotentials` configuration must list the specific, recommended SSSP (Standard Solid State Pseudopotentials) precision library file names for both Iron (Fe) and Platinum (Pt).
    *   Still under `dft_compute`, the `ecutwfc` (wavefunction cutoff energy) value must be a sensible floating-point number that is appropriate for the selected high-precision pseudopotentials (e.g., a value greater than or equal to 80 Rydberg).
    *   Crucially, under `dft_compute`, the `magnetism` setting should be automatically set to a default starting guess, such as `ferromagnetic`, because the heuristic engine correctly detected the presence of a known magnetic element (Fe) in the input.

### Scenario: UAT-C02-03 - Robust Validation of User Input and Graceful, Informative Error Handling

*   **GIVEN** a clean project directory.
*   **AND** the user deliberately creates an `input.yaml` file containing a common typographical error in a required key:
    ```yaml
    system:
      elemants: ["Si"]  # Intentional typo: should be "elements"
      composition: "Si"
    ```
*   **WHEN** the user attempts to execute the command `uv run mlip-pipe input.yaml`.
*   **THEN** the program must terminate immediately with a non-zero exit code, indicating failure.
*   **AND** a clear, concise, and helpful error message must be printed to the console. The message should ideally diagnose the specific problem, for example: "Validation Error: a required field 'system.elements' is missing" or "Validation Error: 'system.elemants' is not a recognised field. Did you mean 'elements'?".
*   **AND** the user must verify that no partial or empty `exec_config_dump.yaml` or `mlip.db` files have been created. The program should have failed before producing any persistent artifacts.
