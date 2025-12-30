# CYCLE 02: SPEC.md

## 1. Summary

This document provides the detailed technical specifications for Cycle 02 of the Autonomous Development Environment (AC-CDD) project. The primary focus of this cycle is to implement the Architect agent and the associated `gen-cycles` command. This represents a significant step forward in the project, as it introduces the first AI-powered component of the system. The Architect agent will be responsible for taking the high-level, raw requirements defined in `ALL_SPEC.md` and transforming them into a comprehensive set of architectural documents, including the `SYSTEM_ARCHITECTURE.md` file, and the `SPEC.md` and `UAT.md` files for each development cycle.

The `gen-cycles` command will serve as the user-facing entry point for this functionality, orchestrating the process of invoking the Architect agent, creating a dedicated Git branch for the design phase, and committing the generated documents to the repository. The implementation will involve integrating with the Jules API, which will power the Architect agent's reasoning and document generation capabilities. The successful completion of this cycle will result in a fully functional architectural design phase, capable of automatically generating a detailed and actionable plan from a set of high-level requirements.

This cycle builds upon the foundational CLI and configuration management system created in Cycle 01. The API keys, which are now securely stored in the `.env` file, will be used to authenticate with the Jules API. The development process will be guided by a robust testing strategy, with a focus on mocking the external API to ensure that the tests are fast, reliable, and do not incur any cost. By the end of this cycle, the AC-CDD system will be able to autonomously generate its own development plans, a core feature of its self-driving capabilities.

## 2. System Architecture

The system architecture for Cycle 02 expands upon the foundation laid in Cycle 01, introducing the first AI agent and the services required to support it.

The new and updated components for this cycle are:

*   **`gen-cycles` Command:** A new command will be added to the Typer-based CLI. This command will be responsible for initiating the architectural design phase. It will delegate the core logic to the Architect Service.

*   **Architect Service:** This new service will encapsulate the business logic for the `gen-cycles` command. It will be responsible for creating a new Git branch, invoking the Architect agent, and committing the generated files to the repository.

*   **Architect Agent:** The core of this cycle's functionality, the Architect agent will be responsible for generating the architectural documents. It will be powered by the Jules API and will be provided with the content of the `ALL_SPEC.md` file as its primary input.

*   **Jules API Client:** A new client class will be created to handle all interactions with the Jules API. This will provide a clean and reusable interface for making API calls, handling authentication, and processing the responses.

*   **Git Service:** A dedicated service will be created to manage all interactions with the Git repository. This will include functionality for creating branches, checking the status of the repository, and committing files. This service-oriented approach ensures that the Git-related logic is centralised and reusable.

The workflow for this cycle is as follows: The user runs the `gen-cycles` command. The command invokes the Architect Service, which in turn calls the Git Service to create a new design branch. The Architect Service then invokes the Architect Agent, which uses the Jules API Client to communicate with the Jules API and generate the architectural documents. Finally, the Architect Service uses the Git Service to commit the newly created documents to the design branch.

## 3. Design Architecture

The design architecture for Cycle 02 focuses on integrating the new components into the existing structure in a clean and modular way.

**File Structure:**

The following new and updated files will be created within `dev_src/ac_cdd_core/`:

*   `cli.py`: Updated to include the new `gen-cycles` command.
*   `agents/architect.py`: A new file to house the `ArchitectAgent` class.
*   `services/architect.py`: A new file for the `ArchitectService` class.
*   `services/git.py`: A new file for the `GitService` class.
*   `clients/jules.py`: A new file for the `JulesApiClient` class.

**Class and Function Definitions:**

*   **`cli.py`:**
    *   `@app.command() def gen_cycles():`: The function that defines the `gen-cycles` command. It will instantiate and run the `ArchitectService`.

*   **`agents/architect.py`:**
    *   `class ArchitectAgent:`:
        *   `def __init__(self, api_client):`: The constructor will take an instance of the `JulesApiClient`.
        *   `def generate_documents(self, requirements):`: The main method, which takes the raw requirements as a string and returns the generated documents.

*   **`services/architect.py`:**
    *   `class ArchitectService:`:
        *   `def run(self):`: The main method that orchestrates the `gen-cycles` process, including creating the Git branch, invoking the agent, and committing the files.

*   **`services/git.py`:**
    *   `class GitService:`:
        *   `def create_branch(self, branch_name):`: Creates a new Git branch.
        *   `def commit_files(self, files, message):`: Commits a list of files with a given message.

*   **`clients/jules.py`:**
    *   `class JulesApiClient:`:
        *   `def __init__(self, api_key):`: The constructor takes the API key.
        *   `def generate(self, prompt, context):`: A method for making a generation request to the Jules API.

This design ensures that each component has a single, well-defined responsibility, which makes the system easier to understand, test, and maintain.

## 4. Implementation Approach

The implementation of Cycle 02 will be carried out in a series of well-defined steps:

1.  **Add `gen-cycles` command:** The `cli.py` file will be updated to include the new `gen-cycles` command, with a placeholder for the `ArchitectService`.

2.  **Implement `GitService`:** The `GitService` class will be implemented in `services/git.py`. This will involve using the `gitpython` library to interact with the Git repository.

3.  **Implement `JulesApiClient`:** The `JulesApiClient` class will be created in `clients/jules.py`. This will use the `httpx` library to make asynchronous requests to the Jules API.

4.  **Implement `ArchitectAgent`:** The `ArchitectAgent` class will be implemented in `agents/architect.py`. It will use the `JulesApiClient` to generate the architectural documents based on the provided requirements.

5.  **Implement `ArchitectService`:** The `ArchitectService` class will be created in `services/architect.py`. This class will orchestrate the entire `gen-cycles` process, integrating the `GitService` and `ArchitectAgent`.

6.  **Integrate service with CLI:** The `gen-cycles` command in `cli.py` will be updated to instantiate and run the `ArchitectService`.

7.  **Write unit tests:** Comprehensive unit tests will be written for each of the new classes. The `JulesApiClient` and `GitService` will be mocked to ensure that the tests are fast and isolated.

8.  **Write integration tests:** Integration tests will be written for the `gen-cycles` command using the `CliRunner`. These tests will mock the external Jules API but will use a real Git repository in a temporary directory to verify the end-to-end functionality.

## 5. Test Strategy

The test strategy for Cycle 02 will focus on ensuring the correctness of the Architect agent and its integration into the CLI.

**Unit Testing Approach:**

*   **`ArchitectAgent`:** The unit tests will verify that the agent correctly formats the prompt and context for the Jules API and that it correctly parses the response. The `JulesApiClient` will be mocked to return a predefined response, allowing us to test the agent's logic in isolation.

*   **`ArchitectService`:** The tests will verify that the service correctly orchestrates the `gen-cycles` process. The `ArchitectAgent` and `GitService` will be mocked to ensure that the tests are focused on the service's own logic.

*   **`GitService`:** The unit tests will verify that the `GitService` correctly interacts with the underlying Git repository. The `gitpython` library will be mocked to avoid making actual changes to the filesystem.

*   **`JulesApiClient`:** The tests will verify that the client correctly formats the API requests and handles different types of responses (e.g., success, error). The `httpx` library will be mocked to avoid making actual network calls.

**Integration Testing Approach:**

The integration tests will use the `CliRunner` to test the `gen-cycles` command.

*   **Happy Path:** The test will run the `gen-cycles` command in a clean Git repository. The Jules API will be mocked to return a set of sample documents. The test will then verify that a new Git branch is created and that the generated documents are correctly committed to it.

*   **Dirty Repository:** The test will run the `gen-cycles` command in a repository with uncommitted changes. It will verify that the command exits with an error message, instructing the user to commit or stash their changes first.

By combining these testing approaches, we can ensure that the Architect agent and the `gen-cycles` command are implemented correctly and are robust to various real-world scenarios.
