# Cycle 5 Specification: User Interface & Finalisation

## 1. Summary

Cycle 5 is the final phase of the MLIP-AutoPipe project, focusing on usability, robustness, and packaging. While the previous cycles built the powerful and autonomous core engine, this cycle ensures that the system is easy to use, provides clear feedback, and is delivered as a professional and polished tool. The main deliverables of this cycle are a clean and intuitive Command-Line Interface (CLI), comprehensive logging, clear error reporting, and final user documentation.

The primary goal is to refine the user experience. The CLI will be the main entry point for users, and it must be well-documented and user-friendly. We will use a modern Python library like Typer or Click to create a structured and self-documenting interface. Logging will be improved to provide clear, actionable information about the pipeline's progress and to help with debugging when things go wrong. Error handling will be made more robust to catch common issues (e.g., missing dependencies, invalid inputs) and provide helpful guidance to the user. Finally, this cycle includes the creation of user-facing documentation, such as a README file and tutorials, to ensure that new users can quickly get started with the system. This final polish is crucial for the project's success, as it makes the powerful technology developed in earlier cycles accessible and reliable for the target audience of materials scientists.

## 2. System Architecture

The architectural focus of this cycle is on the outermost layer of the application—the user interface.

**Command-Line Interface (CLI):**
The CLI will be the single entry point for the entire application. It will be designed to be simple and declarative.
-   **Main Command:** A single main command, `mlip-pipe`, will be the entry point.
-   **Argument:** The primary argument will be the path to the user's `input.yaml` file.
-   **Example Usage:** `uv run mlip-pipe input.yaml`
-   **Options/Flags:** The CLI will include options for controlling the application's behaviour, for example:
    *   `--verbose` or `-v`: To increase the logging verbosity for debugging.
    *   `--log-file PATH`: To specify a file for saving the logs.
    *   `--resume PATH`: An advanced feature to resume a previously failed run from a specific state.
-   **Help Text:** The CLI will have auto-generated, helpful text (e.g., `mlip-pipe --help`) that clearly explains the commands and options.

**Logging System:**
A structured logging system will be implemented to provide insight into the pipeline's execution.
-   **Log Levels:** It will use standard log levels (INFO, DEBUG, WARNING, ERROR).
    *   `INFO`: Default level, showing major progress steps (e.g., "Starting Cycle 2: Structure Generation", "DFT calculation for 10 structures complete").
    *   `DEBUG`: More verbose output for developers, showing details of internal operations.
-   **Console Output:** By default, logs will be printed to the console with clear formatting and timestamps. Progress bars for long-running tasks (like MD simulations or model training) will be used to improve user experience.
-   **File Output:** Logs will also be saved to a file (`mlip-autopyke.log`) in the output directory for each run, creating a permanent record for debugging and traceability.

**Error Handling:**
The application will have a global exception handler to ensure that unexpected errors are caught and presented to the user in a friendly way, avoiding scary stack traces for common issues. For example, if a required external program like Quantum Espresso is not found in the system's PATH, the application will exit gracefully with a clear message: "Error: Quantum Espresso executable 'pw.x' not found. Please ensure it is installed and in your PATH."

## 3. Design Architecture

The implementation will focus on a new `cli.py` module and enhancements to the main application orchestrator to integrate logging and error handling.

**`src/mlip_autopyke/cli.py`:**
This file will be the main entry point defined in `pyproject.toml`.
-   **Library:** We will use `Typer` to create the CLI application. Typer is modern, based on Python type hints, and generates help messages automatically.
-   **`main` function:**
    ```python
    import typer
    from .workflow import WorkflowOrchestrator

    app = typer.Typer()

    @app.command()
    def run(
        input_file: Path = typer.Argument(..., help="Path to the minimal input.yaml file."),
        verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable debug logging."),
        log_file: Path = typer.Option(None, help="Path to save the log file.")
    ):
        # 1. Setup logging configuration based on verbose flag and log_file path.
        # 2. Instantiate WorkflowOrchestrator.
        # 3. Wrap the main execution in a try...except block for graceful error handling.
        try:
            orchestrator.execute_full_pipeline(input_file)
        except Exception as e:
            # Log the error and exit with a user-friendly message.
            typer.echo(f"An unexpected error occurred: {e}")
            raise typer.Exit(code=1)
    ```

**`src/mlip_autopyke/workflow/orchestrator.py`:**
The `WorkflowOrchestrator` will be modified to accept a logger instance.
-   **`__init__(self, logger)`**: The constructor will store the logger.
-   **Logging Calls:** Throughout the execution methods, calls to the logger will be added to report progress, e.g., `self.logger.info("Starting structure generation...")`.
-   **Progress Bars:** For long-running loops, like iterating through DFT calculations, we will use a library like `rich` or `tqdm` to display a progress bar in the console.

**Documentation:**
-   **`README.md`:** A high-quality README file will be created at the root of the project. It will include:
    *   A brief project description.
    *   Installation instructions.
    *   A "Quick Start" guide showing a basic usage example.
    *   A link to more detailed documentation.
-   **`docs/` directory:** A directory for more detailed user documentation will be created. This may include tutorials for different material types and a more detailed explanation of the configuration options.

## 4. Implementation Approach

1.  **Integrate Typer:** Add `typer` as a dependency. Refactor the main entry point of the application into the `cli.py` file using the `@app.command()` decorator.
2.  **Setup Logging:** Integrate Python's built-in `logging` module. Write a setup function that configures the log format, level, and handlers (for console and file output) based on the CLI options.
3.  **Refactor Orchestrator:** Pass the configured logger into the `WorkflowOrchestrator` and sprinkle informative log messages throughout the code.
4.  **Add Progress Bars:** Identify the key long-running loops in the pipeline. Integrate the `rich.progress` library to wrap these iterators and provide visual feedback to the user.
5.  **Implement Global Error Handling:** Add the main `try...except` block in `cli.py` to catch exceptions from the pipeline. Write specific checks for common pre-flight errors, such as missing files or external dependencies.
6.  **Write `README.md`:** Draft a clear and concise `README.md` file following the structure outlined above.
7.  **Create User Tutorials:** Write at least one tutorial (e.g., in `docs/tutorial_alloy.md`) that walks a user through the entire process of generating an MLIP for a simple alloy system, from creating the `input.yaml` to interpreting the final output.
8.  **Final Review and Packaging:** Review all code and documentation for clarity and correctness. Ensure the `pyproject.toml` is correctly configured for packaging and distribution if the project were to be published.

## 5. Test Strategy

Testing for Cycle 5 focuses on the user interface and the application's robustness.

**Unit Testing Approach (Min 300 words):**
We will write unit tests for the CLI itself. `Typer` provides a `CliRunner` utility that allows you to invoke CLI commands within a test script and assert on their output and exit codes. We will write tests to verify that:
-   Running with `--help` prints the help message and exits with code 0.
-   Running with a non-existent input file prints a "File not found" error and exits with a non-zero code.
-   Running with the `--verbose` flag correctly sets the logging level to DEBUG.
We will also write unit tests for any specific error handling logic. For instance, we can mock a scenario where a required program is missing and assert that our pre-flight check function correctly raises a custom `DependencyNotFoundError`.

**Integration Testing Approach (Min 300 words):**
The integration tests will focus on the end-to-end user experience. We will have a set of small, complete example `input.yaml` files in the test suite. An integration test will use the `CliRunner` to execute the full `mlip-pipe run` command on one of these examples. The test will run a simplified, fast version of the entire pipeline. The assertions will focus on the user-facing output. The test will capture the console output and assert that the expected progress messages are logged. It will check that a log file is created if the `--log-file` option is used. Finally, it will assert that the command finishes with an exit code of 0, indicating success. Another integration test will be designed to fail intentionally (e.g., by providing an invalid parameter in the `input.yaml`). This test will assert that the application exits with a non-zero code and that a user-friendly error message is printed to the console, rather than a full Python stack trace.
