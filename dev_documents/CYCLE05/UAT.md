# Cycle 5 User Acceptance Tests: User Interface and Finalisation

**Version:** 1.0.0
**Status:** Final

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 5. After four cycles of intensive development on the core scientific capabilities of MLIP-AutoPipe, this final set of tests focuses exclusively on the **end-user experience**. The goal is to validate that the final, packaged software is not just powerful, but also **usable, accessible, and robust** from the perspective of a typical user—a computational materials scientist who is an expert in their domain but not necessarily in software engineering. These scenarios will test the intuitiveness of the Command Line Interface (CLI), the clarity and usefulness of the documentation, the quality of the real-time feedback during a run, and the overall process of installation and distribution. Passing these UATs is the final quality gate, ensuring that the project is ready to be released and can be successfully adopted by its target audience.

| Scenario ID | Description | Priority |
| :--- | :--- | :--- |
| **UAT-C5-001** | **Successful Installation via Package Manager:** This test validates the most common first step in any user's journey: installation. It ensures that the project is correctly packaged and that the installation process is smooth, automatic, and error-free within a standard Python environment. A seamless installation is critical for a positive first impression and for lowering the barrier to entry for new users. | **High** |
| **UAT-C5-002** | **Clear and Usable CLI:** This test focuses on the primary human-computer interface of the application. The CLI must be easy to navigate and understand. This scenario verifies that the command structure is logical, that the built-in help messages are genuinely helpful, and that a user can discover how to perform a standard task simply by interacting with the CLI itself. This is a key test of the application's overall usability. | **High** |
| **UAT-C5-003** | **Informative Progress and Logging:** The pipeline performs long, complex calculations that can take hours or days. A "black box" experience is unacceptable. This test ensures that the user is kept informed of the system's progress and current status. By validating the presence and clarity of real-time feedback mechanisms like progress bars and structured logging, it confirms that the user can effectively monitor their jobs, estimate completion times, and understand what the autonomous system is doing at any given moment. | **High** |
| **UAT-C5-004**| **Comprehensive Documentation:** Good documentation is a critical feature of any successful software project. This test validates that the documentation is accurate, complete, and effective. It simulates the journey of a brand-new user, from installation to running their first successful calculation, based solely on the provided documentation. This ensures that the documentation is a reliable and sufficient resource for user self-service and learning. | **Medium** |
| **UAT-C5-005** | **Graceful Handling of User Errors:** Users will inevitably make mistakes, such as providing a wrong file path or a malformed configuration file. This test ensures the system responds to these common errors in a helpful and non-intimidating way. It validates that the CLI can catch user errors early, fail fast, and provide clear, actionable error messages, rather than crashing with an obscure internal stack trace. This is crucial for building a positive and frustration-free user experience. | **Medium** |

---

## 2. Behavior Definitions

### **UAT-C5-001: Successful Installation via Package Manager**

**Scenario:** A new PhD student has been told to use MLIP-AutoPipe for their project. They have a standard university-provided computing environment with Python and pip installed. They want to install the software in a clean virtual environment to avoid conflicts.
> **GIVEN** a standard Linux environment with Python 3.10+ and `pip` available.
> **AND** the user has created and activated a new, empty Python virtual environment using `python -m venv my_env`.
> **AND** this environment has no previous version of MLIP-AutoPipe or its specific dependencies installed.
>
> **WHEN** the user runs the standard pip command from the terminal: `pip install mlip-autopipe` (assuming the package is hosted on PyPI, or `pip install /path/to/wheel.whl` for a local file).
>
> **THEN** the `pip install` command should execute and complete without any dependency resolution errors or other installation failures, and it should exit with a success status.
> **AND** all of the project's required dependencies (such as `ase`, `torch`, `numba`, `typer`, `rich`) should be automatically downloaded and installed alongside the main package.
> **AND** after the installation is complete, running the command `mlip-pipe --version` from the terminal should immediately print the correct, latest version number of the software.
> **AND** running `which mlip-pipe` should show that the executable has been placed correctly in the virtual environment's `bin` directory.

### **UAT-C5-002: Clear and Usable CLI**

**Scenario:** A user has successfully installed the software but is unsure of the exact command or options needed to start a run. They want to explore the CLI's functionality to figure out how to use it.
> **GIVEN** the MLIP-AutoPipe software is successfully installed and the `mlip-pipe` command is available on the PATH.
>
> **WHEN** the user, seeking an overview of the tool, runs `mlip-pipe --help` in their terminal.
>
> **THEN** the output must display a clear, well-formatted, and comprehensive help message.
> **AND** this top-level message must list all available subcommands, clearly indicating that `run` is the primary command for executing a workflow.
> **AND** when the user, seeking more detail, runs the command `mlip-pipe run --help`.
> **THEN** the output must display a new, specific help message for the `run` command.
> **AND** this message must list and provide a concise, one-line explanation for all the specific options available for this command, such as `--input` (describing the path to the YAML file) and `--output-dir` (describing the path to the results directory).

### **UAT-C5-003: Informative Progress and Logging**

**Scenario:** A user has launched a long active learning simulation that is expected to run for several hours. They want to be able to check on its progress without having to manually inspect log files.
> **GIVEN** a pipeline run has been initiated that involves a long surrogate MD simulation (e.g., millions of steps) and a subsequent active learning loop.
>
> **WHEN** the user observes the terminal output where the pipeline is running.
>
> **THEN** during the surrogate MD simulation, the console must display a dynamic, continuously updating progress bar from the `rich` library.
> **AND** this progress bar must show the percentage of MD steps completed, the number of steps processed per second (iteration rate), and a reasonably accurate estimated time remaining.
> **AND** during the active learning phase, key milestones must be clearly logged to the console with color-coding and timestamps. For example, messages like "INFO: High uncertainty detected. Pausing simulation." or "INFO: Retraining MLIP, cycle 3/10..." should be clearly visible.
> **AND** upon the successful completion of the entire run, a final summary `rich.table.Table` must be printed to the console, showing key results like the total number of DFT calculations performed, the number of active learning cycles completed, and the file path of the final trained model.

### **UAT-C5-004: Comprehensive Documentation**

**Scenario:** A new postdoctoral researcher, who is an expert in materials science but new to this specific software, needs to get up to speed quickly. They need to install the software, run a tutorial example, and then understand a specific configuration parameter for their own research.
> **GIVEN** the user has been given the URL to the project's documentation website (e.g., hosted on GitHub Pages).
>
> **WHEN** the user navigates to the "Installation" section of the documentation.
> **THEN** they should find clear, step-by-step instructions that allow them to successfully install the software and its prerequisites in their environment.
>
> **AND WHEN** the user proceeds to the "Tutorial" section, which provides an example `input.yaml` file.
> **THEN** by following the instructions and using the provided example file, they should be able to successfully run their first calculation and see the expected output files being created.
>
> **AND WHEN** the user, now preparing their own input file, wants to understand the `uncertainty_threshold` parameter, they navigate to the "Configuration Reference" page.
> **THEN** they should be able to easily locate the `uncertainty_threshold` parameter in the reference list, and find a clear explanation of what the parameter does, its data type (float), and its default value.

### **UAT-C5-005: Graceful Handling of User Errors**

**Scenario:** An experienced user is working quickly and makes a common mistake: they try to launch a run but provide the path to an input file that does not exist.
> **GIVEN** the software is installed correctly.
> **AND** the file `my_non_existent_input.yaml` does not exist in the current directory.
>
> **WHEN** the user runs the command `mlip-pipe run --input my_non_existent_input.yaml`.
>
> **THEN** the program must not proceed to the complex core logic. It must fail fast, exiting within a second or two.
> **AND** it must not print a long, intimidating Python stack trace to the console.
> **AND** instead, it must print a single, clear, user-friendly error message, formatted by the CLI framework. For example: "Error: Invalid value for '--input': File 'my_non_existent_input.yaml' not found." This allows the user to immediately diagnose and correct their own mistake without needing to understand the internal workings of the program.
