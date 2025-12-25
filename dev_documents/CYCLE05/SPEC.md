# Cycle 5 Specification: User Interface and Finalisation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

This document provides the detailed technical specification for Cycle 5, the final development phase of the MLIP-AutoPipe project. Having successfully implemented the core scientific engine (Cycle 1), automated data generation (Cycle 2), efficient exploration (Cycle 3), and the autonomous active learning loop (Cycle 4), the project is now a technologically complete and powerful scientific tool. However, power without usability is of limited value. The objective of this final cycle is to address all aspects of the user experience, transforming the collection of sophisticated backend modules into a polished, professional, and accessible software package. The focus shifts from the development of new scientific capabilities to the crucial tasks of **human-computer interaction, documentation, and software distribution**.

The primary deliverable of this cycle is a robust, intuitive, and self-documenting Command Line Interface (CLI). This CLI will serve as the single, authoritative entry point for all user interactions with the pipeline. It will move beyond a simple script-like execution model to a well-structured application with clear subcommands, arguments, and options, built using a modern framework like Typer. A significant effort will be dedicated to providing rich, real-time feedback to the user. Long-running computational tasks will be accompanied by progress bars, and the logging output will be structured and colorized for readability, ensuring that the user is always aware of the system's status and actions. Beyond the CLI itself, this cycle encompasses the creation of a comprehensive documentation suite. This is not merely an afterthought but a core feature, designed to lower the barrier to entry for new users. The documentation will include a quick-start guide, a detailed step-by-step tutorial, and a complete reference for all configuration parameters. Finally, this cycle will address the practicalities of software distribution. The project will be properly packaged according to modern Python standards (PEP 621), enabling straightforward installation with a single `pip` command. At the end of Cycle 5, MLIP-AutoPipe will not just be a functional proof-of-concept; it will be a complete, well-documented, and easily distributable product, ready for its initial release to the broader scientific community.

## 2. System Architecture

The architecture of Cycle 5 is unique in that it introduces no new core scientific modules. Instead, it focuses on building a sophisticated, user-facing layer that sits "on top" of the existing pipeline. The main architectural change is the formalization of the `cli.py` module as the exclusive entry point for users. This design choice enforces a clean separation between the user interface and the core application logic (the `WorkflowManager` and the scientific modules). This separation is crucial for maintainability; changes to the CLI's appearance or command structure can be made with minimal risk of breaking the underlying scientific code. The introduction of a dedicated `analyze` subcommand also hints at a future architectural direction where the CLI can be used not just to run the pipeline, but also to interact with its results in a more sophisticated way. The use of the `rich` library for output is another key architectural decision, centralizing the responsibility for user-facing presentation in the CLI layer and keeping the core modules focused on their computational tasks.

```mermaid
graph TD
    A[User at Command Line] --> B{CLI Application (`mlip-pipe`)};

    subgraph "CLI Subcommands"
        B -- `run --input ...` --> C[Run Command Logic];
        B -- `analyze --db ...` --> D[Analysis Command Logic];
        B -- `--help` --> E[Help Message Display];
    end

    C -- Initiates & Configures --> F[Workflow Orchestrator (from previous cycles)];
    F -- Executes --> G[Core Pipeline (Modules A-E)];
    G -- Produces --> H[ASE Database & MLIP Model];

    D -- Reads from --> H;
    D -- Generates & Displays --> I[Results Summary / Plots];

    subgraph "New Architectural Layer in Cycle 5"
        direction LR
        B; D; E;
    end
```

**Detailed Workflow Description:**

-   **CLI as the Sole Entry Point:** All user interaction with the pipeline is now channeled through a single, well-defined CLI application, which will be installed on the user's system as `mlip-pipe`. This replaces any previous script-based execution methods.
-   **Structured Subcommands:** The CLI's functionality will be organized using a subcommand structure for clarity and extensibility.
    -   The primary function of running the entire pipeline will be handled by a `run` subcommand. A typical user invocation will be `mlip-pipe run --input input.yaml --output-dir /path/to/results`. This command will be responsible for parsing these arguments, performing initial validation (e.g., ensuring the input file exists), and then instantiating and launching the main `WorkflowManager`.
    -   A new `analyze` subcommand will be introduced for post-processing. A user could, for example, run `mlip-pipe analyze --db /path/to/results/results.db` to get a summary of a completed run, such as the number of structures, the distribution of energies, or the final model's performance on the training set.
-   **Rich Logging and Reporting:** The CLI layer will take full responsibility for the user's view of the process. It will initialize and configure the logging system (e.g., setting the log level based on a `--verbose` flag). It will use a library like `rich` to render live progress bars for long-running, iterative tasks, such as the surrogate MD in Module B or the retraining cycles in Module E. This provides crucial real-time feedback, transforming the user experience from a "black box" execution to a transparent process. Upon successful completion of a run, the CLI will call a reporting function to display a clean, well-formatted summary of the key results.
-   **Packaging and Distribution:** The `pyproject.toml` file will be finalized to include an `[project.scripts]` entry point. This tells the packaging tools to create the `mlip-pipe` executable and place it on the user's PATH during installation, making the CLI accessible system-wide. This completes the transition from a development project to a distributable software package.

## 3. Design Architecture

This cycle's design is focused on the top-level application layer, user documentation, and packaging standards.

**Key Classes and Modules:**

-   **`src/mlip_autopipe/cli.py`:** This file becomes the application's main entry point and the heart of the user interface.
    -   The implementation will be based on the **Typer** library, which is chosen for its modern, type-hint-driven approach that automatically generates help messages and performs command-line argument validation.
    -   **`@app.command() run(input_file: Path, ...)` function:** This will be the decorated function for the `run` subcommand. It will take arguments for all user-configurable aspects of a run, such as the input file path, the output directory, and the logging level. Its logic will be lean, focused on validating inputs and then delegating the actual work to the `WorkflowManager`.
    -   **`@app.command() analyze(...)` function:** The decorated function for the `analyze` subcommand. It will contain the logic for connecting to an existing ASE database and generating summary statistics or plots.
    -   **Global Exception Handler:** A top-level `try...except` block will be wrapped around the main application logic. This is a critical feature for user-friendliness. If any exception propagates up from the core pipeline, this handler will catch it, log the full stack trace to a debug file, and present the user with a simple, clean error message indicating what went wrong at a high level (e.g., "Error: DFT calculation failed. Check the log file for details.") rather than the intimidating default Python stack trace.

-   **Project Root Directory:**
    -   **`README.md`:** This file will be significantly overhauled to become a comprehensive, user-focused landing page for the project. It will include a project summary, a list of key features, clear and concise installation instructions, and a "Quick Start" example that allows a new user to get their first successful run in minutes.
    -   **`docs/` directory (New):** This new directory will house the source files for the full user documentation. The documentation will be built using a modern static site generator like **MkDocs** with the `mkdocs-material` theme for a professional and searchable website.
        -   `docs/index.md`: The landing page of the documentation site.
        -   `docs/installation.md`: Detailed installation instructions, including prerequisites like Quantum Espresso and setting up the Python environment with `uv`.
        -   `docs/tutorial.md`: A detailed, narrative-style, step-by-step tutorial that walks the user through a complete, realistic example, from creating the `input.yaml` to analyzing the final results.
        -   `docs/configuration.md`: A crucial reference document that provides a detailed explanation for every single parameter available in the `input.yaml` file, including its type, default value, and effect on the pipeline.
        -   `docs/architecture.md`: A high-level overview of the system's architecture for advanced users who want to understand how the pipeline works under the hood.

## 4. Implementation Approach

The implementation of Cycle 5 is a series of well-defined steps to add the final layer of polish and professionalize the project.

1.  **CLI Framework Integration:**
    a.  The first step is to add `typer` and `rich` to the project's dependencies in `pyproject.toml`.
    b.  The existing, likely simple, script-based entry point will be completely refactored. A new Typer application object will be created in `src/mlip_autopipe/cli.py`.
    c.  The `run` command will be implemented as a Typer-decorated function. The existing workflow orchestration logic will be moved into a `WorkflowManager` class, and the `run` function's body will be simplified to just: parse arguments, instantiate `WorkflowManager`, and call `workflow_manager.execute()`. This enforces the separation of UI and logic. Options like `--input-file`, `--output-dir`, and `--log-level` will be added using Typer's declarative syntax.
2.  **Rich Logging and Reporting:**
    a.  The standard Python `logging` module will be configured to use `rich.logging.RichHandler`. This will immediately make all existing log messages more readable with color-coding and better formatting.
    b.  Long-running loops within the core modules (e.g., the surrogate MD in `ExplorerSampler`, the active learning loop in `SimulationEngine`) will be identified and refactored. The `rich.progress.Progress` context manager will be wrapped around these loops. This will automatically display a dynamic progress bar in the terminal, showing percentage completion, iteration rate, and estimated time remaining, which vastly improves the user experience for long jobs.
    c.  A new reporting function will be created. At the end of a successful run, the `WorkflowManager` will call this function, which will use `rich.table.Table` to print a clean, formatted summary of the key results of the run to the console.
3.  **Documentation Writing:**
    a.  The `README.md` file will be rewritten from scratch to be a user-centric guide.
    b.  A MkDocs project will be initialized in the `docs/` directory.
    c.  The content for the documentation will be written section by section. This is a significant effort that requires clear, concise, and accurate technical writing. The tutorial, in particular, must be carefully crafted and tested to ensure it works perfectly for a new user.
4.  **Packaging for Distribution:**
    a.  The `pyproject.toml` file will be reviewed and finalized to ensure all required metadata (author, license, project URLs, etc.) is present.
    b.  The `[project.scripts]` section will be added to define the `mlip-pipe` entry point, linking it to the main function in `src/mlip_autopipe/cli.py`.
    c.  The build process will be tested by running `python -m build` to create the source distribution (`.tar.gz`) and the wheel (`.whl`) files.
    d.  The most critical step, the installation test, will be performed. This involves creating a completely new, empty virtual environment, installing the newly built wheel file into it using `pip`, and then verifying that the `mlip-pipe` command is available on the PATH and that running `mlip-pipe --help` executes correctly. This test will be automated as part of the CI pipeline.

## 5. Test Strategy

The testing for Cycle 5 is fundamentally different from previous cycles. It focuses less on the scientific correctness of the algorithms and more on the robustness and usability of the user-facing application layer.

**Unit Testing Approach:**
-   **CLI Commands:** The CLI is the primary new component, and it will be tested thoroughly using `typer.testing.CliRunner`. This tool allows us to invoke the CLI application from within a pytest script and assert on its output and exit code, without needing to use `subprocess`. We will create a suite of tests that simulate various user interactions. For example, a test for invalid input will run the equivalent of `mlip-pipe run --input non_existent_file.yaml` and assert that the runner's exit code is non-zero and that its standard output contains the expected "File not found" error message. Another test will run `mlip-pipe run --help` and assert that the output contains the names of all the available options, ensuring the help message is being generated correctly. We will also test the successful parsing of arguments. For instance, we will run with `--log-level DEBUG` and then use mocking to inspect the logger configuration and assert that its level was correctly set to DEBUG.

**Integration Testing Approach:**
-   **Installation and Execution Test:** The single most important test for Cycle 5 is the end-to-end installation and execution test. This test provides the ultimate validation that the final product, as a user would experience it, is functional. This test will be a critical part of the CI/CD pipeline and will be configured to run in a completely sterile, clean environment on the CI runner. The test script will perform the following sequence of operations:
    1.  First, it will invoke the build process (`python -m build`) to create the project's distributable wheel file.
    2.  Next, it will create a new, empty virtual environment using `venv` or `virtualenv`.
    3.  It will then use the `pip` from this new environment to install the wheel file created in the first step. This tests that the package's declared dependencies are correct and can be resolved by `pip`.
    4.  Once installed, the test will execute the installed CLI command, `mlip-pipe run ...`, on a simple, fast-running test case (likely using a mock DFT engine to ensure it completes in seconds).
    5.  Finally, the test will assert that the command completes with an exit code of 0 and that the expected output files (a model and a database) are created in the correct location.
    The successful completion of this entire sequence provides a very high degree of confidence that the packaging is configured correctly, all necessary dependencies are properly specified in `pyproject.toml`, and the final installed product is functional "out of the box."
