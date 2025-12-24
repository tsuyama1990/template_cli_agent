# CYCLE05 Specification: Finalisation, UI, and Documentation

## 1. Summary

CYCLE05 is the final phase of the MLIP-AutoPipe project, focusing on transitioning the powerful backend system developed in previous cycles into a polished, user-friendly, and robust tool. The primary objectives of this cycle are to enhance the user experience (UX), harden the system against real-world usage, and produce comprehensive documentation. While the core scientific functionality is now complete, this cycle is critical for making the tool accessible, reliable, and maintainable.

Key initiatives for this cycle include a significant overhaul of the Command-Line Interface (CLI). The current CLI is functional but lacks the polish expected of a finished product. We will enhance it with features like rich progress bars for long-running tasks (e.g., DFT calculations, MD simulations), more detailed and structured logging, and clearer error reporting. This will provide users with better visibility into the system's status and help them diagnose issues more effectively.

Another major focus is hardening the system. This involves comprehensive error handling for all foreseeable failure modes, such as issues with external dependencies (Quantum Espresso, LAMMPS), file system problems (permissions, disk space), and invalid user inputs. The goal is to replace cryptic stack traces with clear, actionable error messages.

Finally, this cycle will deliver a complete set of documentation. This will not only include technical, API-level documentation for developers but also user-focused guides and tutorials. A "Getting Started" guide will walk new users through the installation process and a basic example. In-depth tutorials will cover more advanced use-cases, such as configuring the active learning loop or choosing the right structure generation strategy. This documentation is essential for the project's long-term success and adoption by the wider scientific community.

## 2. System Architecture

The architectural changes in CYCLE05 are not focused on adding new modules but on improving the interfaces and interactions between existing components and the user. The core logic remains the same, but it will be wrapped in a more robust and communicative shell.

The main architectural enhancements are:
1.  **CLI Refinement**: The `cli.py` module will be refactored to use a more advanced CLI library like `Rich` in conjunction with `Typer`. This will allow for the implementation of progress bars, spinners for waiting periods, and color-coded, formatted output.
2.  **Structured Logging**: A centralized logging system will be implemented. Instead of simple `print` statements, all modules will use Python's `logging` library. Log messages will be configured to output to both the console (with user-friendly formatting) and a detailed log file (for debugging). Log levels (INFO, WARNING, ERROR) will be used to control verbosity.
3.  **Centralized Exception Handling**: A global exception handler will be added to the main CLI entry point. This will catch any unhandled exceptions that occur during the workflow, log the full stack trace to the debug log file, and present a clean, user-friendly error message to the console, often with a suggestion on how to fix the problem.
4.  **Documentation Generation**: A documentation generator like `MkDocs` or `Sphinx` will be integrated into the project. The codebase will be annotated with docstrings (e.g., in Google or NumPy format), which the generator will use to automatically create an API reference. The user guides and tutorials will be written as separate Markdown files within a `docs/` directory.

The overall system diagram does not change, but the user's perception of and interaction with the system will be substantially improved.

## 3. Design Architecture

The focus of this cycle is on improving existing code and adding new user-facing and documentation-related infrastructure.

**Key Files and Classes to be Modified/Added:**

*   `src/mlip_autoprope/cli.py`:
    *   This file will be heavily refactored. `Typer` commands will be decorated to provide better help messages.
    *   The `rich` library will be used to create a `Console` object. `print` statements will be replaced with `console.log` and `console.print`.
    *   `rich.progress.Progress` will be used to wrap long loops, such as the iteration over structures in the `LabelingEngine` or the steps in the `SimulationEngine`.
*   `src/mlip_autoprope/core/logging.py`:
    *   A new file to configure the application's logger. It will set up handlers for both console and file output, with different formatting and levels for each.
    *   All other modules will import and use the logger configured here (`import logging; logger = logging.getLogger(__name__)`).
*   `src/mlip_autoprope/core/exceptions.py`:
    *   A new file to define custom exception classes, e.g., `DFTCalculationError`, `ConfigurationError`, `LAMMPSExecutionError`.
    *   Modules will raise these specific exceptions, allowing the global handler in `cli.py` to catch them and provide targeted error messages.
*   `docs/`:
    *   A new top-level directory.
    *   `docs/index.md`: The main landing page for the documentation.
    *   `docs/getting_started.md`: Installation and first-run tutorial.
    *   `docs/tutorials/`: A subdirectory for more detailed guides.
    *   `mkdocs.yml`: The configuration file for the `MkDocs` site generator.
*   `pyproject.toml`: New development dependencies will be added, including `rich`, `mkdocs`, and any necessary `MkDocs` plugins.

## 4. Implementation Approach

1.  **Dependency Setup**: Add `rich` and `mkdocs` to the `[tool.uv.dev-dependencies]` section of `pyproject.toml`.
2.  **Logging Implementation**:
    *   Create the `core/logging.py` module and implement the setup logic.
    *   Go through all existing modules (`a_...`, `b_...`, etc.) and replace all `print()` statements with `logger.info()`, `logger.warning()`, etc.
3.  **CLI Overhaul**:
    *   Refactor `cli.py` to use a `rich.console.Console` object for all output.
    *   Identify the main loop in the `LabelingEngine` and wrap it in a `rich.progress.Progress` context manager to show a bar tracking the DFT calculations.
    *   Do the same for the main MD/kMC loop in the `SimulationEngine`. Add a spinner for the training step, which is long but not iterative.
4.  **Exception Handling**:
    *   Define the custom exception classes in `core/exceptions.py`.
    *   Go through the modules and replace generic `raise Exception(...)` calls with the new custom exceptions. For example, in the `LabelingEngine`, if QE fails, `raise DFTCalculationError(...)`.
    *   Wrap the main logic in `cli.py` in a `try...except` block that catches these custom exceptions and prints friendly messages.
5.  **Documentation Writing**:
    *   Set up the `MkDocs` site structure with `mkdocs.yml` and the initial Markdown files.
    *   Write the content for the `getting_started.md` guide. This should be a complete, self-contained tutorial.
    *   Write at least one advanced tutorial, for example, on "Understanding and Configuring the Active Learning Loop".
    *   Go through the codebase and add comprehensive docstrings to all public classes and methods, explaining their purpose, arguments, and return values.
6.  **Final Review and Testing**: Perform a full, end-to-end run of the system, paying close attention to the user experience. Check for clarity in log messages, usefulness of progress bars, and correctness of documentation. Fix any typos or confusing language.

## 5. Test Strategy

Testing for CYCLE05 is primarily focused on the user interface and experience, as well as the robustness of the error handling.

### Unit Testing Approach

*   **Custom Exceptions**: While hard to unit test the raising itself, we can test the error handling logic. Write tests that call functions that are mocked to raise specific custom exceptions. The test should assert that the CLI's global exception handler catches these exceptions and returns the expected user-friendly error message. For example, mock the `LabelingEngine` to raise a `DFTCalculationError` and check that the CLI prints "A DFT calculation failed. Check the log file for details." and exits gracefully.
*   **CLI Output**: Use a library like `pytest-capture` to capture the standard output of CLI commands. Write tests that run commands and then assert that the captured output contains the expected text, such as the headers from `rich.table` or the completion message.

### Integration Testing Approach

*   **"Bad Weather" Testing**: This is the core of the integration testing for this cycle. Create a series of tests that intentionally break the system's dependencies to verify the error handling.
    *   **Missing QE**: Run the workflow with the path to `pw.x` pointing to a non-existent file. Assert that the system fails quickly with a clear error message like "`pw.x` not found at the specified path."
    *   **Missing LAMMPS**: Do the same for the `lammps` executable in an active learning run.
    *   **Invalid Config**: Run the workflow with an `input.yaml` that is syntactically correct but semantically wrong (e.g., specifying an element that doesn't exist). Assert that the Pydantic validation catches this and provides a clear error.
    *   **No Disk Space**: (Hard to test automatically, may be a manual test). Run a large simulation in an environment with a very small disk quota. Verify the system fails gracefully when it can no longer write to the trajectory or database files.
*   **Documentation Build**: Add a step to the continuous integration (CI) pipeline that runs `mkdocs build`. This test will fail if the documentation has syntax errors or broken links, ensuring the docs are always in a deployable state.
*   **Full "Happy Path" UAT**: Perform one final end-to-end run of the entire system for a simple case (e.g., generating a potential for Si). This UAT will be evaluated based on the user experience. Are the progress bars accurate? Is the logging clear? Is the final output easy to understand? This serves as the final sign-off on the project's usability.
