# CYCLE05/SPEC.md

## 1. Summary

This document provides the detailed technical specification for Cycle 5, the final development cycle for the initial version of the MLIP-AutoPipe project. The preceding cycles have successfully built a fully autonomous, active-learning-driven pipeline for generating high-quality MLIPs. The objective of this final cycle is to harden the system, enhance its physical realism by incorporating advanced simulation methods and more sophisticated physics, and improve its overall usability and polish. This cycle focuses on transitioning the system from a functional prototype to a robust, feature-complete tool ready for practical use by materials scientists.

The scope of this cycle is twofold. First, we will enhance the simulation capabilities of **Module E (Simulation Engine)** by implementing more advanced techniques, such as Kinetic Monte Carlo (kMC), to allow for the study of long-timescale phenomena like diffusion and phase transitions. Second, we will improve the physical fidelity of our DFT calculations in **Module C (Labeling Engine)** by adding automated handling of complex but critical phenomena like magnetism and spin-orbit coupling (SOC). This is crucial for accurately modeling a wide range of important functional materials. Finally, this cycle will also address overall system usability, including improved documentation, tutorials, and a more polished user interface and output. By the end of this cycle, MLIP-AutoPipe will not only be a powerful, automated tool but also a reliable and accessible platform for the broader materials science community, fully realizing the project's vision of democratizing the creation of high-fidelity machine learning potentials.

## 2. System Architecture

The architecture in Cycle 5 involves targeted enhancements to existing modules rather than the introduction of new ones. The focus is on increasing the sophistication of the simulation and labeling engines and ensuring the overall system is robust and user-friendly.

**File Structure:**

The following ASCII tree highlights the files that will be significantly modified or added in this cycle. The changes are concentrated in Modules C and E and their corresponding utility and configuration files.

```
.
├── docs/
│   ├── **index.md**
│   └── **tutorial.ipynb**  # New user tutorial
├── src/
│   └── mlip_autopipec/
│       ├── data/
│       │   └── models.py     # Modified to include advanced physics/sim config
│       ├── modules/
│       │   ├── c_labeling_engine.py # Modified for magnetism/SOC
│       │   └── e_simulation_engine.py # Modified for kMC
│       ├── orchestrator.py # Minor modifications for new options
│       └── utils/
│           ├── dft_utils.py  # Modified to generate QE inputs with mag/SOC
│           └── **kmc_utils.py**  # New utilities for kMC simulations
├── tests/
│   ├── unit/
│   │   ├── modules/
│   │   │   └── **test_e_simulation_engine_kmc.py**
│   │   │   └── **test_c_labeling_engine_magnetism.py**
│   │   └── utils/
│   │       └── **test_kmc_utils.py**
│   └── e2e/
│       └── **test_cycle05_workflow.py**
└── ...
```

**Component Blueprint:**

*   **`modules/e_simulation_engine.py` (Modified)**: The `SimulationEngine` will be extended to support different simulation modes. Alongside the existing MD simulation for active learning, it will gain the capability to run Kinetic Monte Carlo (kMC) simulations. This will likely involve adding a new method, e.g., `run_kmc_simulation()`, which would use the trained MLIP to calculate energy barriers for atomic events (like vacancy hops) and simulate the system's evolution over much longer timescales than MD.
*   **`utils/kmc_utils.py`**: This new utility file will contain the core logic for the kMC algorithm. This includes functions for finding transition states (e.g., using a Nudged Elastic Band method with the MLIP), calculating event rates using Eyring's transition state theory, and running the main kMC loop (often using the BKL or residence-time algorithm).
*   **`modules/c_labeling_engine.py` (Modified)**: The `LabelingEngine` will be upgraded to handle more complex physics. The logic for generating Quantum Espresso inputs will be enhanced to include parameters for magnetism (e.g., setting initial magnetic moments on atoms) and spin-orbit coupling (SOC). It will also include a "two-stage" learning strategy, where initial training is done without SOC for speed, and a final refinement stage is performed with SOC enabled for maximum accuracy.
*   **`utils/dft_utils.py` (Modified)**: The helper functions in this module will be updated to support the new features in the `LabelingEngine`. For example, `create_qe_input_from_atoms` will now accept additional arguments to control magnetic and SOC settings in the generated input file.
*   **`data/models.py` (Modified)**: The Pydantic configuration models will be extended with new sections to allow the user to enable and control these advanced features (e.g., a boolean `soc_enabled` flag, a string for `magnetism_initialization`).
*   **`docs/`**: A new top-level directory will be created to house user-facing documentation, including a comprehensive tutorial in the form of a Jupyter Notebook that walks a user through a complete example of generating a potential.

## 3. Design Architecture

The design for Cycle 5 focuses on providing user-configurable "switches" to enable advanced, computationally expensive features. This allows the user to make a trade-off between speed and physical accuracy, which is a critical aspect of practical computational science.

**Pydantic Schema Design:**

*   **`MagnetismConfig` (Model)**: A new model to group magnetism settings.
    *   `initialization_strategy`: A string literal, e.g., `"ferromagnetic"`, `"antiferromagnetic"`, `"random"`.
    *   `soc_enabled`: A boolean flag to enable/disable spin-orbit coupling calculations.
    *   `refinement_stage_only`: A boolean to control the "two-stage" learning. If true, SOC will only be enabled in the final active learning generations.

*   **`SimulationModeConfig` (Model)**: A model to control the behavior of Module E.
    *   `mode`: A string literal, e.g., `"md_active_learning"`, `"kmc_analysis"`.
    *   `kmc_parameters`: An optional nested model containing settings specific to kMC (e.g., `event_types`, `temperature`).

*   **`FullConfig` (Top-level Model, modified)**:
    *   `dft_compute`: This section will be updated to include an instance of `MagnetismConfig`.
    *   `simulation`: This section will be updated to include an instance of `SimulationModeConfig`.

**Data Flow and Consumers/Producers:**

The core data flow remains the active learning loop from Cycle 4. The new features augment this flow:

*   **Consumer (`LabelingEngine`)**: The engine now consumes the `MagnetismConfig` from the main configuration. When generating QE inputs, it uses these settings to add the appropriate flags and parameters, thus **producing** more physically complex DFT calculations. The results (labeled structures) are still produced in the same format.
*   **Consumer (`SimulationEngine`)**: The engine now consumes the `SimulationModeConfig`.
    *   If `mode` is `"md_active_learning"`, the data flow is the same as in Cycle 4.
    *   If `mode` is `"kmc_analysis"`, the engine still **consumes** the trained MLIP. However, instead of producing uncertain structures for retraining, it **produces** analysis data, such as a list of diffusion events, mean squared displacement plots, or other long-timescale properties. This represents a terminal, "production science" step in the workflow.

This design allows the same core pipeline to be used for both creating a potential (the active learning mode) and using it for scientific discovery (the analysis mode), all controlled via the central configuration file.

## 4. Implementation Approach

The implementation will be split between enhancing Module C for better physics and Module E for advanced simulation capabilities.

1.  **Configuration and Heuristics**:
    *   Add the new `MagnetismConfig` and `SimulationModeConfig` models to `data/models.py`.
    *   Update the `ConfigExpander` service. For example, it should be able to heuristically enable the `magnetism` settings if it detects known magnetic elements like "Fe", "Co", or "Ni" in the user's input.

2.  **Advanced DFT Features (Module C)**:
    *   Modify `utils/dft_utils.py`. The `create_qe_input_from_atoms` function will be refactored to accept the `MagnetismConfig` and add the necessary lines to the QE input string (e.g., `nspin = 2`, `starting_magnetization`, `lspinorb = .true.`).
    *   Modify the `LabelingEngine` class. It will need to read the new configuration and pass it down to the utility functions. For the "two-stage" learning, it will need to be aware of the current active learning generation number to decide whether to enable SOC.
    *   Write unit tests in `test_c_labeling_engine_magnetism.py` that create an `Atoms` object for a magnetic material like iron and assert that the generated QE input string contains the correct magnetism-related keywords.

3.  **Advanced Simulation Features (Module E)**:
    *   Create `utils/kmc_utils.py`. This is a significant undertaking. The first implementation might focus on a specific type of event, like vacancy-mediated diffusion. It will need functions to find nearest neighbors, use the MLIP to estimate energy barriers (e.g., using a simplified drag method), and the core BKL algorithm for the kMC loop.
    *   Refactor the `SimulationEngine` class. The `execute` method will now have a conditional block based on the `simulation.mode` config setting. If it's "md", it runs the existing logic. If it's "kmc", it calls the new kMC-related methods.
    *   Unit tests in `test_e_simulation_engine_kmc.py` will be challenging. They will likely involve using a mock potential with a very simple, known energy landscape (e.g., a 1D sine wave). The test would then run the kMC logic and assert that the simulated event rates match the analytically calculated rates for that simple potential.

4.  **Documentation and Usability**:
    *   Create the `docs/` directory.
    *   Write a `tutorial.ipynb` Jupyter Notebook. This notebook will guide the user through a complete, simple example (e.g., generating a potential for bulk copper). It will include cells for setting up the `input.yaml`, running the `mlip-pipe` command, and then loading the final potential to calculate a simple property like the lattice constant, thereby validating the result.
    *   Improve the CLI output with better progress indicators and more informative logging.

5.  **End-to-End Test**:
    *   Create `test_cycle05_workflow.py`. This E2E test will specifically validate the new configuration flags. One test case will run with magnetism enabled and assert that the (mocked) call to Quantum Espresso includes the correct command-line flags or input file parameters. Another test case will set the simulation mode to "kmc" and assert that the (mocked) kMC methods in `SimulationEngine` are called instead of the MD methods.

## 5. Test Strategy

The testing for Cycle 5 focuses on verifying the correctness of the new, complex physical models and simulation methods.

**Unit Testing Approach (Min 300 words):**
*   **Magnetism and SOC**: The unit tests for the new DFT features will be centered in `test_c_labeling_engine_magnetism.py`. We will not test Quantum Espresso itself, but rather our interface to it. A key test will involve creating an ASE `Atoms` object for bcc Iron. We will then call our updated `create_qe_input_from_atoms` function with a `MagnetismConfig` that specifies a ferromagnetic initialization. The test will then assert that the resulting output string contains the exact QE keywords required, such as `nspin = 2` and the correct `starting_magnetization(1) = ...` line. Another test will do the same but with `soc_enabled = True`, and it will assert that `lspinorb = .true.` is present in the output. This confirms that our configuration is being correctly translated into the format the DFT engine expects.
*   **Kinetic Monte Carlo**: Testing the kMC logic in `test_kmc_utils.py` requires a controllable potential. We will create a mock MLIP calculator that implements a simple, 2D potential surface with a known minimum and a known energy barrier for an "atom" to hop from one site to another. We will then run our kMC simulation on this mock potential. According to transition state theory, the rate of hopping should be proportional to `exp(-Barrier/kT)`. Our test will run the kMC simulation for a set number of steps and count the number of hops that occurred. We will then assert that this observed number of hops is statistically consistent with the rate predicted by the analytical formula. This provides a strong validation that our implementation of the kMC algorithm and its use of the MLIP's energy predictions are physically correct.

**Integration Testing Approach (Min 300 words):**
*   **Configuration-driven Logic Branching**: An integration test will be designed to verify that the `Orchestrator` and the `SimulationEngine` correctly change their behavior based on the configuration file. We will run the pipeline with a config file that sets `simulation.mode = "kmc_analysis"`. We will use mock versions of all the modules. The test will assert that the `execute` method of the `SimulationEngine`'s MD part is *not* called, while the `run_kmc_simulation` method *is* called. Conversely, we will run with `mode = "md_active_learning"` and assert the opposite. A similar test will be done for the `LabelingEngine`, where we will check that the `soc_enabled` flag correctly propagates from the config file down to the `dft_utils` that generate the input files. This confirms that our high-level configuration switches are correctly controlling the program's execution flow.
*   **End-to-End Test for Advanced Features**: The `test_cycle05_workflow.py` will simulate a user invoking the CLI with a configuration that uses the new features. It will be a "dry run" test with extensive mocking. One test scenario will involve running with magnetism enabled. We will mock the `subprocess.run` call for Quantum Espresso. The assertion will be to inspect the mocked call's arguments to ensure that the input file passed to it was indeed generated with the correct magnetic parameters. Another scenario will test the "two-stage" learning for SOC. We will mock the `Orchestrator`'s loop and the `LabelingEngine`. The test will assert that during the first few mocked "generations," the labeling engine is called with `soc_enabled=False`, and in the later generations, it is called with `soc_enabled=True`. This provides end-to-end validation that these complex, multi-stage workflow configurations are being respected by the pipeline.
