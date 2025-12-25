# Specification: Cycle 5 - Advanced Simulation & Finalisation

**Version:** 1.0.0
**Status:** Final

## 1. Summary

Cycle 5 is the final development phase of the MLIP-AutoPipe project, culminating in a feature-complete, robust, and user-friendly system. Having established the core automated pipeline and the transformative active learning loop in the preceding cycles, this final cycle focuses on expanding the system's scientific capabilities and ensuring it is ready for real-world use. The primary objective is to implement advanced simulation techniques, specifically those for exploring rare events, which are often the most scientifically interesting phenomena. A secondary but equally important goal is to finalize the user experience, including a polished command-line interface (CLI), comprehensive documentation, and a suite of example use cases.

The major technical scope of this cycle is the extension of `Module E: Simulation Engine` to support not just molecular dynamics but also more advanced sampling methods like Adaptive Kinetic Monte Carlo (akMC). While MD is excellent for studying processes that occur on nanosecond timescales, many critical material processes, such as diffusion, phase nucleation, and chemical reactions, are "rare events" that may only occur once per microsecond or longer. To address this, we will implement an akMC framework. This involves integrating saddle point search algorithms (like Nudged Elastic Band or the Dimer method) to find transition states, a rate calculator to determine event probabilities, and a KMC engine to simulate the long-timescale evolution of the system. This entire process will be accelerated by JIT-compiling the performance-critical KMC loop with Numba.

The second part of the scope is finalization and polish. We will refine the CLI to be more intuitive and informative, providing clear feedback to the user about the progress of the workflow. We will write comprehensive user documentation, explaining not just how to run the software but also the underlying concepts and best practices. Finally, we will prepare a set of tutorials and example workflows for different material systems to help new users get started quickly. By the end of Cycle 5, MLIP-AutoPipe will not just be a functional proof-of-concept; it will be a powerful, production-ready tool capable of tackling a wide range of challenging problems in computational materials science. This will fulfill the project's ultimate goal of providing a truly "expert-in-a-box" solution for MLIP generation. The final deliverable will be a polished, well-documented, and scientifically powerful tool that is ready for distribution to the wider research community.

## 2. System Architecture

The architectural changes in Cycle 5 are primarily extensions to `Module E: Simulation Engine`, with a focus on enhancing its capabilities and performance. The overall cyclical architecture of the active learning loop remains the same, but the simulation module itself will become more sophisticated.

**Architectural Extensions to `Module E`:**
`Module E` will be refactored to support multiple simulation modes. The user will be able to specify in the `input.yaml` whether they want to run a standard MD simulation or an advanced akMC simulation.
*   **Mode Dispatcher:** The `SimulationEngine`'s main `run_otf_simulation` method will act as a dispatcher. Based on the configuration, it will delegate the execution to either the existing MD runner or the new akMC runner.
*   **akMC Runner Component:** This new component will encapsulate the entire logic for the Adaptive Kinetic Monte Carlo simulation. Its architecture will be a state machine that cycles through several steps:
    1.  **Saddle Point Search:** From the current state (a local minimum on the potential energy surface), it will launch a series of saddle point searches to find nearby transition states. This will use algorithms like the Dimer method. Crucially, the uncertainty quantification (UQ) mechanism developed in Cycle 4 will remain active during these searches. If the system enters an uncertain region while climbing towards a saddle point, it will trigger the re-training loop, ensuring the energy barriers are calculated accurately.
    2.  **Rate Calculation:** Once a set of escape pathways (saddle points) is found, it will calculate the transition rates for each event using Harmonic Transition State Theory.
    3.  **KMC Step:** It will then use the calculated rates to perform a standard Kinetic Monte Carlo step, selecting one event to execute and advancing the simulation time by the appropriate amount. This step is computationally intensive and will be heavily optimized.
    4.  **State Update:** The system moves to the new minimum, and the cycle repeats.

**Performance Optimization with Numba:**
A key architectural feature of this cycle is the systematic application of Just-In-Time (JIT) compilation using Numba. The KMC event loop, which involves searching and updating event lists and calculating rates, can be a significant bottleneck in pure Python. The core logic of this loop will be factored out into dedicated functions in a `utils` file and decorated with `@numba.jit(nopython=True)`. This will compile the Python code down to highly efficient machine code, enabling simulations to reach timescales that would be impossible with an interpreted implementation. This represents a commitment to high performance as a core architectural principle of the system.

## 3. Design Architecture

The design for Cycle 5 involves refactoring `Module E` to accommodate the new akMC functionality and creating a more polished user-facing CLI.

**Updated Project Structure:**
```
src/mlip_autoflow/
├── __init__.py
├── main.py                # CLI logic will be enhanced
├── config/
│   └── models.py
└── modules/
    ├── ...
    └── e_simulation_engine/   # Converted to a sub-package
        ├── __init__.py
        ├── main_engine.py     # The top-level dispatcher
        ├── md_runner.py       # Existing MD logic refactored here
        └── akmc_runner.py     # New file for akMC logic
    └── utils/
        └── kmc_kernels.py     # New file for Numba-optimized KMC loop
docs/
├── user_guide.md
├── tutorials/
│   └── si_bulk_example.md
```

**Class and Method Definitions:**

*   **`config.models.py`**: The `SimulationConfig` Pydantic model will be updated.
    *   A `simulation_mode: Literal['md', 'akmc']` field will be added.
    *   New sub-models like `AKMCConfig` will be added to hold parameters for saddle point searches (e.g., `dimer_method_params`) and the KMC simulation itself.

*   **`modules.e_simulation_engine/`**: This module is now a sub-package.
    *   **`main_engine.py`**: The `SimulationEngine` class will be refactored.
        *   The `run_otf_simulation` method will now read the `simulation_mode` and instantiate either `MDRunner` or `AKMCRunner`, passing the relevant configuration.
    *   **`md_runner.py`**: The `MDRunner` class will contain the logic previously in the `SimulationEngine`, cleaned up and focused solely on running molecular dynamics.
    *   **`akmc_runner.py`**: The new `AKMCRunner` class will be the core of this cycle's work.
        *   `run(self)`: The main loop that orchestrates the akMC workflow (saddle search -> rate calculation -> KMC step).
        *   `_find_saddle_points(self, current_state)`: A method that uses ASE's optimizers (e.g., `ASE.optimize.Dimer`) to find transition states. It will incorporate the UQ check.
        *   `_calculate_rates(self, saddle_points)`: A method to compute rates using Harmonic TST.
        *   `_perform_kmc_step(self, rates)`: A method that calls the high-performance, Numba-compiled KMC kernel to choose an event and advance time.

*   **`modules.utils.kmc_kernels.py`**: This new file will contain the core KMC algorithm.
    *   A function `execute_kmc_step(rates_array)` decorated with `@numba.jit` that implements the "residence time" algorithm efficiently using NumPy arrays.

*   **`main.py`**: The Typer-based CLI will be enhanced.
    *   The output will be made more user-friendly, using a library like `rich` to display progress bars for long-running processes like DFT calculations and simulations.
    *   Help messages and command documentation will be improved.
    *   A new command like `mlip-pipe generate-docs` could be added to create sphinx-based documentation.

*   **`docs/`**: A new top-level directory will be created to house the user documentation and tutorials.

This refactored design makes `Module E` much more modular and extensible for future simulation methods, while the focus on the CLI and documentation in `main.py` and `docs/` directly addresses the usability and finalization goals of the cycle.

## 4. Implementation Approach

The implementation of Cycle 5 will be split into two parallel efforts: the development of the akMC functionality and the finalization of the user-facing components.

**Step 1: Refactor `Module E`**
The first step is to perform the architectural refactoring. We will convert the `e_simulation_engine.py` file into a sub-package. The existing MD logic will be moved into the `md_runner.py` file and encapsulated in an `MDRunner` class. The `SimulationEngine` class in `main_engine.py` will be updated to act as the dispatcher. This refactoring will be done first to ensure that all existing functionality from Cycle 4 continues to work before we add new features.

**Step 2: Implement Saddle Point Search**
We will begin the implementation of the `AKMCRunner`. The first component will be the `_find_saddle_points` method. We will use the saddle point search algorithms available in ASE, such as the Dimer method. A key task here is to integrate the uncertainty quantification (UQ) check from Cycle 4 into the optimization loop of the saddle point search. As the optimizer takes steps "uphill" towards a transition state, we will perform a UQ check at each step. If the uncertainty is high, we will pause the search, trigger the re-training of the MLIP, and then resume the search with the improved potential.

**Step 3: Implement Rate Calculation**
Once a saddle point is found, we need to calculate the transition rate. We will implement the `_calculate_rates` method. This will involve calculating the vibrational frequencies at both the initial minimum and the saddle point to get the prefactors for Harmonic Transition State Theory. This is a standard procedure and can be implemented using ASE's `Vibrations` class.

**Step 4: Implement and Optimize the KMC Loop**
This is a performance-critical step. We will first write a pure Python/NumPy implementation of the residence time KMC algorithm in `modules/utils/kmc_kernels.py`. This algorithm involves calculating a cumulative sum of the rates and drawing a random number to select an event. Once the Python version is working and tested, we will apply the `@numba.jit(nopython=True)` decorator to it. We will benchmark the performance before and after to quantify the speedup. The `_perform_kmc_step` method in the `AKMCRunner` will then simply be a wrapper that calls this highly optimized kernel.

**Step 5: Finalize the CLI and User Experience**
In parallel with the akMC development, another developer can focus on the user-facing aspects. The `main.py` script will be updated to use the `rich` library to provide better real-time feedback. This includes:
*   Using `rich.progress` to show progress bars for loops over structures or simulation steps.
*   Using `rich.panel` to display a summary of the configuration and results.
*   Improving the clarity and formatting of log messages.

**Step 6: Write Documentation and Tutorials**
The final step of the project is to write the documentation. We will create the `docs` directory. The `user_guide.md` will be a comprehensive manual explaining the philosophy of the software, the structure of the `input.yaml` file, and the meaning of the outputs. The `tutorials` directory will contain step-by-step guides for running a complete workflow on 1-2 different material systems (e.g., creating a potential for bulk Silicon, and another for an Fe-Pt alloy). These tutorials will serve as both examples and a final end-to-end test suite for the entire system.

## 5. Test Strategy

The testing for Cycle 5 needs to validate the complex logic of the akMC algorithm and ensure the overall system is robust and user-friendly.

**Unit Testing Approach:**
*   **`AKMCRunner`**: Testing the full akMC loop is difficult in a unit test. We will test its components in isolation.
    *   **Saddle Point Search:** We will test the `_find_saddle_points` method on a simple analytical potential (like the Müller-Brown potential) where the location of the saddle points is known. We will mock the MLIP calculator and UQ check, and assert that the method correctly drives the system towards the known saddle point.
    *   **Rate Calculation:** We will provide fixed frequencies for a minimum and a saddle point and assert that the `_calculate_rates` method returns the correct rate according to the TST formula.
    *   **KMC Kernel:** We will write a dedicated test for the Numba-compiled `execute_kmc_step` function. We will provide a fixed set of rates and assert that over many runs, the function selects the different events with the correct probabilities.

**Integration Testing Approach:**
The final integration test will be a full, end-to-end run of the akMC active learning workflow on a simple, well-known system.
*   **Test Scenario:** Simulating vacancy diffusion in a simple crystal (e.g., Aluminum).
    1.  **Setup:**
        *   Create an `input.yaml` configured for an `akmc` simulation.
        *   Start with a training set that is known to be incomplete (e.g., it describes the perfect crystal well but has no information about atomic configurations near a vacancy).
    2.  **Execution:**
        *   Run the entire MLIP-AutoPipe workflow.
        *   The initial training will produce a potential that is poor at describing the saddle point for vacancy hopping.
        *   The `AKMCRunner` will start, and during its first saddle point search for the hop, the UQ mechanism should be triggered.
        *   The system should automatically label the saddle point structure, re-train the potential, and then resume the akMC simulation.
        *   The simulation should then be able to correctly identify the diffusion event and its rate.
    3.  **Verification:**
        *   Check the logs to confirm that the active learning loop was triggered during the saddle point search.
        *   Verify that the final calculated diffusion barrier is close to the known DFT value for the system.
        *   Confirm that the KMC simulation successfully simulates a hop and updates the system state.

This test validates the most complex feature of the entire project. Successfully passing this test, along with completing the documentation and CLI polishing, will signify the successful completion of the MLIP-AutoPipe project.
