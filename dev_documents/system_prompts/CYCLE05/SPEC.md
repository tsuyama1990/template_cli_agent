# SPEC.md: Cycle 05 - Advanced Simulations & User Interface

## 1. Summary

Cycle 05 represents the final and culminating phase of the MLIP-AutoPipe project, a cycle dedicated to polishing the powerful backend engine into a robust, user-friendly, and distributable software package ready for widespread scientific research. The development efforts in this cycle are strategically focused on two primary and equally important objectives. The first is to significantly enhance the scientific capabilities of the **Simulation Engine (Module E)** by integrating advanced, long-timescale simulation techniques, specifically Adaptive Kinetic Monte Carlo (kMC). This new functionality is designed to enable the study of critical material phenomena like diffusion, phase nucleation, and chemical reactions involving rare events, which are fundamentally inaccessible to the shorter timescales of standard Molecular Dynamics (MD). The second major objective is to develop a clean, professional, and intuitive **Command-Line Interface (CLI)** that makes the full, sophisticated power of the automated pipeline accessible and manageable for the end-user.

The integration of an Adaptive kMC algorithm will address a key limitation of MD simulations, which are highly effective at exploring local energy basins but are often trapped in them, unable to observe the rare activated processes that govern many important material transformations. By implementing a kMC engine that uses the trained MLIP to discover saddle points on the potential energy surface, the system will be able to efficiently discover transition pathways and predict their corresponding reaction rates. To make this computationally tractable, a "Tiered Rate Calculation" approach will be implemented. This smart, multi-stage process will ensure that the most expensive parts of the calculation are reserved for only the most probable events, dramatically expanding the scope of scientific questions the generated MLIPs can answer, pushing simulation timescales into the realm of microseconds or even longer.

On the user-experience front, the development of a polished, intuitive, and self-documenting CLI using the `click` library is paramount. This will replace any temporary, ad-hoc developer scripts with a professional, stable, and versioned interface. The user will be able to initiate and manage the entire complex workflow, from initial generation to active learning and advanced kMC simulations, with a small set of simple, intuitive commands. This cycle also explicitly includes the critical "last mile" tasks of software engineering: the creation of comprehensive user documentation, the writing of clear tutorials and examples, and the final packaging of the software for public distribution. Upon the successful completion of Cycle 05, the MLIP-AutoPipe will have transitioned from a powerful technology into a finished product: a fully-documented, easily-installable, and accessible software tool ready to be deployed by the broader materials science community.

## 2. System Architecture

The architectural changes in Cycle 05 are highly focused and represent the final layers of functionality and user interaction. The modifications are concentrated in two main areas: a significant functional enhancement of the `Simulation Engine (Module E)` to support the new kMC simulation mode, and the introduction of a formal `main.py` CLI entry point that will serve as the single, unified interface for the entire application.

The new workflow for an advanced kMC simulation will be as follows:
1.  **Initiation via CLI**: The user, wishing to study a long-timescale process, will initiate a kMC simulation via a new, dedicated CLI command, for example, `mlip-pipe run-kmc --config input.yaml`.
2.  **Loading of Trained MLIP**: The `Orchestrator`, triggered by the CLI, will first load the most mature, fully-trained MLIP from the project's database. This assumes that a robust potential has already been generated using the active learning workflow.
3.  **kMC Simulation Start**: The `Orchestrator` will then call the `Simulation Engine`, but this time with a specific mode or method that invokes the kMC logic instead of the standard MD loop.
4.  **Saddle-Point Search (Event Discovery)**: The kMC algorithm, running within the Simulation Engine, will begin by searching for possible transition events from the system's current atomic state. This involves using specialized algorithms, such as the Dimer method, which leverage the MLIP's energy and force predictions to efficiently find the saddle points on the potential energy surface that connect the current energy minimum to adjacent ones.
5.  **Tiered Rate Calculation**: For each discovered saddle point, the system will calculate the transition rate using Transition State Theory. To do this efficiently, it will employ a two-tiered approach:
    *   **Tier 1 (Screening)**: A quick, approximate calculation of the rate for all potential events is performed using a fixed, physically reasonable frequency prefactor (e.g., 10^12 Hz). This cheap calculation allows the system to rapidly identify a small number of events that are most likely to occur.
    *   **Tier 2 (Refinement)**: For only this small subset of high-probability events, a more expensive Hessian calculation is performed (again, using the MLIP as the calculator) to determine an accurate vibrational frequency prefactor for each specific event, thus significantly refining the rate calculation.
6.  **Event Selection and System Evolution**: Based on these refined rates, a single event is stochastically selected to be the "next" event in the simulation. The atomic configuration is updated to the final state of the chosen transition, and the simulation time is advanced by an amount inversely proportional to the sum of all rates.
7.  **Looping for Long Timescales**: The process then repeats from step 4, continuing to discover and execute events, allowing the system to simulate timescales that are many orders of magnitude longer than what is possible with direct MD. Uncertainty quantification can still be integrated into this process, particularly during the saddle-point searches, to trigger retraining if the model is found to be unreliable in these critical transition-state regions.

The CLI architecture will be designed for clarity and extensibility:
*   The file `src/mlip_autopipec/main.py` will be established as the single, canonical entry point for all user interactions with the application.
*   It will leverage the `click` library to define a main command group (`mlip-pipe`) and a series of clear, verb-oriented subcommands (e.g., `run` for the main workflow, `run-kmc` for specialised simulations, and potentially `analyze` or `status` in the future).
*   Each subcommand will be responsible for parsing its own specific command-line arguments and options. Its sole action will be to instantiate the `Orchestrator` and call the appropriate high-level method from its public API. This cleanly and effectively decouples the user interface and argument parsing logic (which resides in `main.py`) from the core business logic of the application (which resides in `orchestrator.py`).

This final architecture completes the system by providing both the advanced scientific functionality required by expert users and the professional, maintainable, and user-friendly interface required for broad adoption.

## 3. Design Architecture

The design for the final cycle, Cycle 05, involves significant functional additions to the `SimulationEngine` to incorporate the kMC logic, and the creation of the final, polished command-line interface. This will also involve updates to the `Orchestrator` to expose the new functionality to the CLI.

**New and Updated Classes/APIs:**

1.  **`mlip_autopipec.modules.simulation_engine.SimulationEngine`** (Updated)
    *   **New Method**: `run_akmc(mlip_model: Any)`: This new public method will be the primary entry point for running an Adaptive Kinetic Monte Carlo simulation. It will contain the main kMC loop that evolves the system through time.
    *   **Key Internal Methods for kMC**:
        *   `_find_saddle_points(atoms: ase.Atoms) -> List[SaddlePointInfo]`: This private method will implement a sophisticated algorithm (such as the Dimer method or a similar saddle-point search technique) to discover the transition states accessible from the current energy minimum. It will return a list of data structures, each containing the atomic configuration of the saddle point and its energy.
        *   `_calculate_rates(saddle_points: List[SaddlePointInfo]) -> Dict[str, float]`: This method will implement the "Tiered Rate Calculation" strategy. It will first perform a cheap screening of all saddle points with a fixed prefactor. Then, for the most promising candidates, it will call the `_calculate_hessian` method to get an accurate prefactor and refine the final rate.
        *   `_calculate_hessian(atoms: ase.Atoms) -> np.ndarray`: This method will calculate the Hessian matrix (the matrix of second derivatives of the energy) at a specific atomic configuration (either a minimum or a saddle point) by using numerical finite differences of the forces provided by the MLIP. This is needed to compute vibrational frequencies.
        *   `_select_and_execute_event(...)`: This method will contain the logic for the stochastic kinetic Monte Carlo algorithm itself. It will select the next event to execute based on the calculated rates and will update the system's atomic state and advance the simulation time accordingly.

2.  **`mlip_autopipec.main.py`** (New File)
    *   **Purpose**: To define the complete, user-facing command-line interface using the `click` library.
    *   **Structure (using `click` decorators)**:
        ```python
        import click
        from .orchestrator import Orchestrator

        @click.group()
        @click.version_option()
        def cli():
            """MLIP-AutoPipe: The Autonomous ML Potential Generator."""
            pass

        @cli.command()
        @click.argument('config_file', type=click.Path(exists=True))
        def run(config_file):
            """Run the full, end-to-end active learning pipeline to generate a new MLIP."""
            orchestrator = Orchestrator(config_file)
            orchestrator.run_active_learning_loop()

        @cli.command()
        @click.argument('config_file', type=click.Path(exists=True))
        def run_kmc(config_file):
            """Run a long-timescale Kinetic Monte Carlo simulation with a pre-trained MLIP."""
            orchestrator = Orchestrator(config_file)
            orchestrator.run_kmc_simulation()

        # Other commands like 'init' or 'status' could be added in the future.
        ```

3.  **`mlip_autopipec.orchestrator.Orchestrator`** (Updated)
    *   **`__init__` Modification**: The constructor will be updated to take the `config_file` path as its primary argument. It will be responsible for loading this file, immediately running the `ConfigExpander` to generate the full configuration, and storing this `FullConfig` object as an instance variable for use by all other methods. This ensures that the orchestrator is always in a valid, fully configured state.
    *   **New Method**: `run_kmc_simulation()`: This new public method will be added to the orchestrator's API. Its responsibility is to load the best available MLIP from the database and then call the `SimulationEngine.run_akmc` method to execute the simulation.

**Documentation and Packaging Strategy**:
*   **`pyproject.toml`**: The `pyproject.toml` file will be finalised. A crucial addition will be the definition of the CLI entry point under the `[project.scripts]` table. For example: `mlip-pipe = "mlip_autopipec.main:cli"`. This is the standard mechanism that allows the `mlip-pipe` command to be automatically available in the user's PATH when the package is installed via `pip` or `uv`.
*   **`docs/` directory**: A new top-level directory, `docs/`, will be created to hold all user-facing documentation. This will be written in a clear, accessible format like Markdown and will include a "Getting Started" guide, a comprehensive explanation of all the options in the `input.yaml` file, and several tutorials for common scientific use cases.
*   **`README.md`**: The main project `README.md` file will be significantly updated and rewritten to be user-facing rather than developer-facing. It will provide a concise, high-level overview of the project's purpose and capabilities, and will prominently link to the full documentation in the `docs/` directory.

This comprehensive design completes the system by adding the final layer of user-facing abstraction (the CLI) and the advanced scientific methods (kMC) that build upon the solid foundation of all the previously developed components.

## 4. Implementation Approach

The implementation for Cycle 05 will be strategically divided into two parallel tracks that can be developed largely independently before a final integration step. This allows for concurrent progress on the complex scientific algorithm development and the more standard software engineering task of building a user interface.

1.  **Track 1: kMC Algorithm Development (in Module E)**:
    *   The first and most challenging step in this track is to research and select a suitable, robust Python library for performing saddle-point searches (like the Dimer method or Nudged Elastic Band). If a suitable library that can interface with ASE calculators is not available, a simplified version of the Dimer method will be implemented from scratch. This is the most complex part of the scientific implementation.
    *   With the saddle-point search method in place, the `_find_saddle_points` method will be implemented.
    *   Next, the `_calculate_hessian` method will be implemented. This will likely involve writing a numerical finite-difference routine that repeatedly calls the MLIP's `get_forces` method with small atomic displacements.
    *   The `_calculate_rates` method will then be implemented to orchestrate the tiered calculation logic: a fast screening followed by a refined calculation for promising candidates.
    *   Finally, the main `run_akmc` loop will be implemented. This will be a stateful method that calls these helper functions in the correct sequence to simulate the system's kinetic evolution over time.

2.  **Track 2: CLI and UX Polish**:
    *   The `click` library will be added to the project's dependencies in `pyproject.toml`.
    *   The `main.py` file will be created, and the `click` command group (`@click.group()`) and the primary `run` subcommand will be defined as specified in the design.
    *   The `Orchestrator`'s `__init__` method will be updated to be initialised with the config file path, immediately loading and expanding the configuration.
    *   The `run` command's implementation will be connected to the `Orchestrator.run_active_learning_loop` method, ensuring the existing functionality is now accessible through the new CLI.

3.  **Integration of kMC into the CLI**:
    *   The two tracks will be merged. The `run-kmc` subcommand will be added to `main.py`.
    *   The corresponding `run_kmc_simulation` method will be implemented in the `Orchestrator`, which will simply call the new `run_akmc` method on the `SimulationEngine`.
    *   The `FullConfig` Pydantic model will be extended to include a new section for kMC-specific parameters, ensuring that these simulations can be configured by the user.

4.  **Final Packaging and Documentation**:
    *   The `pyproject.toml` file will be updated to formally register the `mlip-pipe` script entry point.
    *   The project will then be built and installed locally using `uv pip install .`. A crucial manual test will be to open a new terminal and confirm that the `mlip-pipe --help` command works as expected.
    *   The `docs/` directory will be created, and the user guide and tutorials will be written. A static site generator like MkDocs may be used to render these Markdown files into a professional-looking HTML documentation website.
    *   Finally, the main `README.md` will be rewritten to be clear, concise, and aimed squarely at new users, providing them with the motivation and initial steps needed to get started with the software.

This dual-track approach is an efficient strategy that allows the two distinct types of work in this cycle to proceed in parallel, with the two tracks converging seamlessly during the final integration and packaging step.

## 5. Test Strategy

The testing strategy for the final cycle, Cycle 05, must be rigorous and multi-faceted, covering the complex scientific logic of the new kMC algorithm, the user-facing functionality and robustness of the CLI, and the successful integration of all components in final end-to-end tests.

**Unit Testing Approach (Min 600 words):**

The unit tests for this cycle will be highly targeted, using mock objects and simplified models to validate the logic of the new components in isolation.

*   **`SimulationEngine` (kMC methods)**: Testing the complex kMC logic with a high-dimensional, black-box MLIP is nearly impossible. Therefore, these tests will rely on using a simple, known, 2D analytical potential (like a Mueller-Brown potential) as a "mock MLIP". This toy potential has well-characterised minima and saddle points, allowing for precise and deterministic testing.
    *   The `_find_saddle_points` method will be tested by programmatically placing an "atom" at one of the known minima of the 2D potential. The test will then execute the saddle-point search algorithm and assert that the coordinates it finds for the nearby saddle point are numerically very close to the known analytical coordinates of that saddle point.
    *   The `_calculate_hessian` method will be tested by comparing its numerical finite-difference result with the exact, analytical second derivative of the 2D toy potential at a given point. This will validate the correctness of the finite-difference implementation.
    *   The `_calculate_rates` method will be tested to ensure that the Arrhenius equation is correctly applied using the energies from the toy potential and that the tiered logic correctly prioritises the lower-energy saddle points.
    *   These tests, while seemingly abstract, provide a very high degree of confidence in the fundamental correctness of the complex kMC implementation before it is used with a real, high-dimensional MLIP.

*   **CLI (`main.py`)**: The `click.testing.CliRunner` is the standard tool for testing `click`-based command-line interfaces, and it will be used extensively here. The tests will be designed to run without actually executing the full, expensive pipeline.
    *   The `Orchestrator` class and all its methods will be completely mocked using `unittest.mock.patch`.
    *   Tests will then invoke the CLI runner with a variety of command-line argument combinations, for example, `runner.invoke(cli, ['run', 'test.yaml'])`.
    *   The tests will then assert that the correct method on the mock `Orchestrator` was called (e.g., `run_active_learning_loop` was called) and that it was called with the correctly parsed arguments (e.g., that the `config_file` path was correctly passed).
    *   A separate suite of tests will be created to validate the CLI's robustness and user-friendliness. These tests will cover user errors, such as providing a path to a non-existent file, providing a malformed config file, or using a typo in a subcommand. The tests will assert that `click`'s error handling works as expected and that the application exits gracefully with a clear, informative, and user-friendly error message in each case.

**Integration and End-to-End Testing (Min 300 words):**

The final testing phase for the entire project will involve running a small number of complete, end-to-end tests. These tests are designed to verify that all the components from all five cycles have been successfully integrated and that the final, user-facing product works as intended.

*   **E2E Test for the `run` command**:
    *   This test will represent the most common user workflow and will be a full system test. It will start with a minimal `input.yaml` file for a very simple and computationally fast-to-calculate system (e.g., a Nitrogen dimer, N2).
    *   The test script will then execute the final `mlip-pipe run input.yaml` command as a subprocess.
    *   This test will use the *real* DFT code (Quantum Espresso) and the *real* MLIP training library, with no mocking of major components.
    *   To be feasible to run in an automated CI/CD environment, the test will be configured to run for only one or two active learning iterations and with very loose convergence parameters.
    *   The final assertion will be simple: that the entire process completes successfully with an exit code of 0 and that a trained MLIP file has been created in the output directory. This single test provides a powerful validation that all cycles (01-05) have been integrated correctly and the main use case is functional.

*   **E2E Test for the `run-kmc` command**:
    *   This test will first need a pre-trained MLIP. It can be designed to use the MLIP generated by the successful completion of the `run` command E2E test.
    *   The test will then invoke the `mlip-pipe run-kmc input.yaml` command.
    *   The test will be configured to run the kMC simulation for a very small number of steps (e.g., 5-10 events).
    *   The primary assertions will be that the process completes successfully and that it has produced a trajectory output file. A more detailed assertion can be made by reading the trajectory and confirming that the system's state has evolved (i.e., that at least one atomic hop or event has occurred and the atoms have moved to a new minimum). This confirms that the kMC machinery is correctly integrated with the rest of the system.
