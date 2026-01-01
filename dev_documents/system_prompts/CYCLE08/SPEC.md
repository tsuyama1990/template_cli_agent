# Specification: CYCLE08 - Advanced Simulation (kMC) & UI/UX

## 1. Summary

CYCLE08 is the final cycle of this project phase, focusing on expanding the capabilities of `Module E: Simulation Engine` to include more advanced simulation techniques and on refining the overall user experience (UX). The primary scientific goal is to implement a Kinetic Monte Carlo (kMC) simulation engine. While MD is excellent for exploring phenomena on the picosecond-to-nanosecond timescale, kMC is essential for reaching the microsecond-to-second timescales required to study rare events like diffusion, phase nucleation, and catalysis. This cycle will deliver a basic kMC engine that uses the trained MLIP to evaluate the energies of different atomic states and predict long-timescale system evolution.

The second major goal of this cycle is to improve the usability and robustness of the command-line interface (CLI) and the system's logging output. This involves providing clearer feedback to the user, better progress indicators for long-running tasks, more informative error messages, and finalising the documentation. This UX focus is critical for turning the powerful backend developed in the previous cycles into a tool that is pleasant and intuitive for the end-user. The successful completion of this cycle will not only add a significant new scientific capability (kMC) to the pipeline but will also polish the entire application into a mature, user-friendly package.

## 2. System Architecture

This cycle enhances the `SimulationEngine` with a new kMC mode and refines the top-level `cli.py` and workflow orchestration to improve user feedback.

**File Structure (CYCLE08 Focus):**

Files to be created or modified are in **bold**.

```
mlip-autopipec/
├── pyproject.toml
├── uv.lock
├── src/
│   └── mlip_autopipec/
│       ├── __init__.py
│       ├── **cli.py**              # Major UX improvements (progress bars, better help)
│       ├── **workflow.py**         # Add logging and orchestrate kMC runs
│       ├── **config.py**           # Add config for kMC simulations
│       ├── database.py
│       └── modules/
│           ├── __init__.py
│           ├── structure_generator.py
│           ├── explorer_sampler.py
│           ├── labeling_engine.py
│           ├── training_engine.py
│           └── **simulation_engine.py** # Add kMC method
└── tests/
    ├── conftest.py
    ├── unit/
    │   └── **test_kmc_engine.py** # New test file
    └── integration/
        └── **test_cli_ux.py**     # New test for CLI experience
```

**Component Breakdown:**

*   **`config.py`**: A new Pydantic model, `KMCSimulationConfig`, will be added.
    *   `KMCSimulationConfig`: Will contain fields like `temperature: float`, `num_events: int`, `event_catalog: list[str]` (e.g., ["vacancy_hop", "adatom_diffusion"]).
*   **`modules/simulation_engine.py`**: The `SimulationEngine` class will be expanded.
    *   A new public method, `run_kmc()`, will be added.
    *   This method will take the `KMCSimulationConfig`. It will first need to identify all possible events and their corresponding transition states from the initial atomic configuration.
    *   It will use the trained MLIP to calculate the energy barrier for each possible event (`E_a`).
    *   It will then use the Arrhenius equation (`k = v * exp(-E_a / k_B * T)`) to calculate the rate of each event. Initially, the pre-exponential factor `v` will be a fixed constant.
    *   The main kMC loop will then be executed: advance time based on event rates, select an event to perform stochastically, update the atomic configuration, and repeat.
*   **`cli.py`**: The `click`-based CLI will be significantly overhauled for better UX.
    *   **Progress Bars:** Long-running processes like exploration, labelling, and training will be wrapped with progress bars (e.g., using the `rich` library). For example, the user will see a bar showing the progress of the MD simulation or the number of structures labelled.
    *   **Better Logging:** Logging will be made more structured. Use a library like `loguru` or `rich` to provide color-coded, timestamped logs with clear levels (INFO, WARNING, ERROR).
    *   **Help Text:** All commands and options will have clear, well-written help text.
*   **`workflow.py`**: The `WorkflowOrchestrator` will be updated to include better logging at the start and end of each major step. It will also be responsible for calling the new `simulation_engine.run_kmc()` method when requested by the CLI.
*   **`tests/unit/test_kmc_engine.py`**: A new test file to unit-test the kMC algorithm's logic.
*   **`tests/integration/test_cli_ux.py`**: A new test file that uses `click.testing.CliRunner` to capture the CLI's output and assert that it contains the expected progress bar and logging formats.

## 3. Design Architecture

The design for this cycle focuses on a simplified but correct kMC implementation and a professional CLI experience.

*   **kMC Engine Design:**
    *   The initial kMC implementation will be a "lattice" kMC, assuming a known crystal structure and a predefined catalog of possible events (e.g., an atom hopping to a vacant lattice site).
    *   **Event Finding:** The engine will first build a neighbor list for the current `ase.Atoms` structure. It will then iterate through the atoms to find all possible "events" based on the `event_catalog`. For a "vacancy_hop", it would find a vacancy and see which of its neighbors could hop into it.
    *   **Barrier Calculation:** For each potential event, a transition state search (e.g., Nudged Elastic Band) would ideally be needed. However, to simplify for this cycle, we can use a heuristic or a pre-computed barrier for a given event type. The MLIP is used to calculate the energy of the initial and final states.
    *   **Rate Calculation:** The rates for all possible events are calculated and stored.
    *   **Event Selection:** The core kMC algorithm uses the calculated rates to determine both *when* the next event will happen and *which* event it will be. This is a standard stochastic algorithm (BKL or rejection-free).
    *   **State Update:** The `ase.Atoms` object is updated by moving the atom that was selected to perform the event. The process repeats.

*   **CLI UX Design (`rich` library):**
    *   The `rich` library will be the standard for all user-facing output.
    *   **`rich.progress`**: Will be used for loops with a known number of iterations (labelling structures, running kMC events). A `Progress` object will be created and updated within the loop.
    *   **`rich.spinner`**: Will be used for tasks of indeterminate length (e.g., SCF convergence in Quantum Espresso). A spinner will be displayed with a status message like "[bold green]Running SCF calculation...[/]"
    *   **`rich.console`**: A global `Console` object will be used for all logging. `console.log()` will be used for standard messages, and `console.print()` with styling will be used for important summaries or warnings (e.g., `console.print("[bold yellow]Warning: SCF failed to converge...[/]")`).
    *   This consistent use of a single, powerful library ensures a professional and uniform user experience across the entire application.

## 4. Implementation Approach

The kMC engine and CLI improvements can be developed semi-independently.

1.  **Dependencies and Config:** Add `rich` (and potentially `loguru`) to `pyproject.toml`. Update `config.py` to include the `KMCSimulationConfig`.
2.  **TDD for kMC:** Create `test_kmc_engine.py`. Write a test for a very simple kMC scenario.
    *   **Test Case:** A 2D square lattice with one vacancy. There are 4 possible hops into the vacancy. Assume two have a low energy barrier and two have a high barrier.
    *   **Mocking:** Mock the MLIP calculator to return the pre-defined energies for the initial, final, and transition states.
    *   **Assertions:** Run the kMC for many steps. Assert that the vacancy has moved. More importantly, assert that the number of times the low-barrier hops occurred is significantly greater than the number of times the high-barrier hops occurred. This validates the stochastic selection algorithm is correctly biased by the energy barriers.
3.  **Implement kMC Engine:** Implement the `run_kmc` method in the `SimulationEngine`. This includes the event-finding logic, rate calculation, and the main BKL/rejection-free kMC loop algorithm.
4.  **Refactor CLI Output:** Go through `cli.py` and `workflow.py`.
    *   Replace all `print()` statements with `console.log()` or `console.print()`.
    *   Identify the main loops (e.g., the loop over structures to be labelled in the orchestrator). Wrap these loops with `rich.progress.track` or a manual `Progress` context manager.
    *   Add spinners for operations like the QE subprocess call.
5.  **TDD for CLI UX:** Create `test_cli_ux.py`.
    *   Write a test that uses `click.testing.CliRunner` to invoke a command that has a progress bar.
    *   Capture the `stdout` of the command.
    *   Assert that the output string contains the special characters and formatting (`"━"`, `"[.../...]"` etc.) that are characteristic of a `rich` progress bar. This confirms the new UX has been successfully implemented.
6.  **Integrate kMC:** Add a new CLI command, `cdd run-kmc`, which calls a new orchestrator method that in turn calls `simulation_engine.run_kmc()`.

## 5. Test Strategy

This final cycle requires testing a new scientific simulation method and the user-facing presentation layer.

**Unit Testing Approach (Min 300 words):**
The unit tests for the kMC engine will be crucial for verifying the correctness of its statistical mechanics implementation.

*   **`test_kmc_engine.py`**:
    *   **`test_rate_calculation`**: This test will focus solely on the rate calculation. It will provide a fixed energy barrier (`E_a`), temperature (`T`), and pre-exponential factor (`v`). It will assert that the calculated rate `k` matches the expected value from the Arrhenius equation. This will be tested for a few different inputs to ensure numerical stability.
    *   **`test_event_selection`**: This is the most important unit test. It will set up a scenario with two possible events. Event A will be configured to have a rate that is exactly 9 times higher than Event B. The kMC loop will be run, for example, 10,000 times. After the loop, the test will count how many times Event A was chosen versus Event B. It will then assert that the ratio of `count(A) / count(B)` is close to 9.0 (e.g., within a 10% tolerance to account for stochastic noise). This provides strong statistical validation of the core event selection algorithm.
    *   **`test_time_advancement`**: This test will verify that the kMC simulation time is being updated correctly. It will run a simulation with known, fixed rates and assert that the total elapsed time in the simulation state matches the expected value `t = -ln(u) / R_total`, where `u` is a random number and `R_total` is the sum of all rates.

**Integration Testing Approach (Min 300 words):**
The integration tests will focus on the CLI experience and the wiring of the new kMC command.

*   **`test_cli_ux.py`**:
    *   **`test_labeling_progress_bar`**: This test will simulate a workflow where three structures need to be labelled. It will mock the `LabelingEngine`'s `label` method. The test will invoke the CLI and capture the output. It will then assert that the output contains strings like `"[1/3]"`, `"[2/3]"`, and `"[3/3]"`, confirming that the progress bar was displayed correctly during the labelling loop.
    *   **`test_error_message_styling`**: This test will configure a mock to raise an exception (e.g., a simulated SCF convergence failure). It will invoke the CLI and capture the output. The test will assert that the captured error message is styled with the expected `rich` formatting for errors (e.g., it contains the ANSI escape codes for red color or is prefixed with a specific error emoji/text).
*   **`test_kmc_workflow_wiring`**: A test in `test_active_learning_loop.py` or a similar file will verify the new `run-kmc` command. It will mock the `SimulationEngine.run_kmc` method. The test will invoke `cdd run-kmc` via the `CliRunner` and assert that the mocked `run_kmc` method was called exactly once, confirming the command is correctly wired to the engine.
