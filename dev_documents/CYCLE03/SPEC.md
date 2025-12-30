# CYCLE 03: SPEC.md

## 1. Summary

This document provides the detailed technical specifications for Cycle 03 of the Autonomous Development Environment (AC-CDD) project. The main focus of this cycle is to implement the Coder agent and integrate the E2B sandboxed execution environment. This is a critical phase of the project, as it introduces the core functionality of autonomous code generation and testing. The Coder agent will be responsible for reading the specifications generated in Cycle 02 and writing the corresponding Python code and unit tests.

A key feature of this cycle is the integration of the E2B sandbox, which will provide a secure and isolated environment for all coding and testing activities. This ensures that the development process is reproducible and does not interfere with the user's local machine. The `run-cycle` command will be introduced as the user-facing entry point for this functionality, allowing users to initiate a development cycle for a specific feature.

This cycle builds upon the work of the previous two cycles. The CLI and configuration system from Cycle 01 will be used to manage the `run-cycle` command and the E2B API key, while the documents generated in Cycle 02 will provide the input for the Coder agent. The successful completion of this cycle will result in a system that can autonomously write and test code in a secure, sandboxed environment, representing a major milestone in the project's development.

## 2. System Architecture

The system architecture for Cycle 03 introduces several new components to support the Coder agent and the sandboxed execution environment.

The new and updated components are:

*   **`run-cycle` Command:** A new command will be added to the CLI to initiate a development cycle. It will take a cycle ID as an argument to specify which feature to implement.

*   **Coder Service:** This new service will encapsulate the business logic for the `run-cycle` command. It will be responsible for setting up the E2B sandbox, invoking the Coder agent, and running the tests within the sandbox.

*   **Coder Agent:** The Coder agent will be responsible for writing the source code and unit tests. It will be powered by an LLM, such as Jules, and will be provided with the `SPEC.md` for the current cycle as its primary input.

*   **E2B Sandbox Service:** A dedicated service will be created to manage all interactions with the E2B sandbox. This will include methods for creating a new sandbox, executing commands within it, and syncing files between the local machine and the sandbox.

*   **Jules API Client:** The existing Jules API client will be used to power the Coder agent.

The workflow for this cycle is as follows: The user runs the `run-cycle` command with a specific cycle ID. The command invokes the Coder Service, which then uses the E2B Sandbox Service to create a new sandbox. The Coder Service then invokes the Coder Agent, providing it with the relevant `SPEC.md`. The Coder Agent writes the code and unit tests and saves them to the local filesystem. The Coder Service then uses the E2B Sandbox Service to sync the new files to the sandbox and run the tests.

## 3. Design Architecture

The design architecture for Cycle 03 focuses on integrating the new components in a modular and maintainable way.

**File Structure:**

The following new and updated files will be created within `dev_src/ac_cdd_core/`:

*   `cli.py`: Updated to include the new `run-cycle` command.
*   `agents/coder.py`: A new file for the `CoderAgent` class.
*   `services/coder.py`: A new file for the `CoderService` class.
*   `services/sandbox.py`: A new file for the `E2bSandboxService` class.

**Class and Function Definitions:**

*   **`cli.py`:**
    *   `@app.command() def run_cycle(cycle_id: str):`: The function that defines the `run-cycle` command. It will instantiate and run the `CoderService`.

*   **`agents/coder.py`:**
    *   `class CoderAgent:`:
        *   `def __init__(self, api_client):`: The constructor will take an instance of the `JulesApiClient`.
        *   `def write_code(self, spec):`: The main method, which takes the cycle specification as a string and returns the generated code and unit tests.

*   **`services/coder.py`:**
    *   `class CoderService:`:
        *   `def run(self, cycle_id):`: The main method that orchestrates the `run-cycle` process, including setting up the sandbox, invoking the agent, and running the tests.

*   **`services/sandbox.py`:**
    *   `class E2bSandboxService:`:
        *   `def __init__(self, api_key):`: The constructor takes the E2B API key.
        *   `def start(self):`: Starts a new sandbox session.
        *   `def execute(self, command):`: Executes a command in the sandbox.
        *   `def sync_files(self, local_path, remote_path):`: Syncs files to or from the sandbox.
        *   `def stop(self):`: Stops the sandbox session.

This design ensures a clear separation of concerns, with the agent responsible for code generation, the sandbox service for environment management, and the coder service for orchestration.

## 4. Implementation Approach

The implementation of Cycle 03 will be carried out in the following steps:

1.  **Add `run-cycle` command:** The `cli.py` file will be updated to include the new `run-cycle` command.

2.  **Implement `E2bSandboxService`:** The `E2bSandboxService` class will be implemented in `services/sandbox.py`. This will involve using the E2B Python SDK to interact with the sandbox API.

3.  **Implement `CoderAgent`:** The `CoderAgent` class will be implemented in `agents/coder.py`. It will use the `JulesApiClient` to generate the code and unit tests based on the provided specification.

4.  **Implement `CoderService`:** The `CoderService` class will be created in `services/coder.py`. This class will orchestrate the entire `run-cycle` process, integrating the `E2bSandboxService` and `CoderAgent`.

5.  **Integrate service with CLI:** The `run-cycle` command in `cli.py` will be updated to instantiate and run the `CoderService`.

6.  **Write unit tests:** Comprehensive unit tests will be written for each of the new classes. The E2B API and Jules API will be mocked.

7.  **Write integration tests:** Integration tests will be written for the `run-cycle` command using the `CliRunner`. These tests will use a mocked E2B sandbox to verify the end-to-end functionality.

## 5. Test Strategy

The test strategy for Cycle 03 will focus on ensuring the correctness of the Coder agent and the sandboxed execution environment.

**Unit Testing Approach:**

*   **`CoderAgent`:** The unit tests will verify that the agent correctly formats the prompt for the Jules API and that it correctly parses the response to extract the code and unit tests. The `JulesApiClient` will be mocked.

*   **`CoderService`:** The tests will verify that the service correctly orchestrates the `run-cycle` process. The `CoderAgent` and `E2bSandboxService` will be mocked.

*   **`E2bSandboxService`:** The unit tests will verify that the `E2bSandboxService` correctly interacts with the E2B SDK. The SDK will be mocked to avoid making actual API calls.

**Integration Testing Approach:**

The integration tests will use the `CliRunner` to test the `run-cycle` command.

*   **Happy Path:** The test will run the `run-cycle` command for a specific cycle. The Jules API will be mocked to return a sample set of code and unit tests. The E2B sandbox will also be mocked to simulate the successful execution of the tests. The test will verify that the command completes successfully and that the generated files are saved to the correct location.

*   **Failing Tests:** The test will simulate a scenario where the unit tests generated by the Coder agent fail in the sandbox. It will verify that the `run-cycle` command exits with an appropriate error message.
