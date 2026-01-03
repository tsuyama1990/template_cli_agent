# Specification: Cycle 2 - Advanced Exploration and Web UI

## 1. Summary

This document outlines the technical specification for the second and final development cycle of the MLIP-AutoPipe project. Building upon the stable foundational architecture established in Cycle 1, this cycle is dedicated to implementing the advanced, scientifically critical components of the pipeline and enhancing its usability. The primary focus is on developing the sophisticated **Exploration Engine**, which represents the core intellectual property of the application. This engine will use a hybrid Molecular Dynamics (MD) and Monte Carlo (MC) approach to generate diverse and physically relevant atomic configurations. The goal is not simply to run a simulation but to intelligently search the vast potential energy surface of a given chemical system. This involves implementing complex logic for managing simulations, handling high-energy events, and automatically adapting the simulation protocol to the type of system being studied. This component is what elevates the tool from a simple structure generator to a powerful engine for scientific discovery, capable of uncovering the rare and high-energy atomic configurations that are essential for training robust and accurate MLIPs.

Furthermore, this cycle will replace the simple random sampling method from Cycle 1 with a more intelligent **Farthest Point Sampling (FPS)** algorithm. While random sampling is adequate, it is not optimal. FPS is a state-of-the-art technique that ensures the final dataset is not just a random collection of points from the simulation, but a structurally diverse and minimally redundant set of configurations. This is achieved by representing each structure in a high-dimensional feature space and selecting points that are maximally distant from one another. This technique is computationally more demanding but yields datasets that are significantly more efficient for training machine learning models, leading to better models with less training data.

Finally, to enhance usability and accessibility for a broader scientific audience, this cycle will introduce a **Web-based User Interface (Web UI)**. While a CLI is powerful for experts and for integration into automated workflows, a graphical interface is invaluable for new users, for teaching, and for interactive exploration. This UI will provide an intuitive, form-based way for users to define their physical system, configure all the complex parameters of the pipeline, launch simulation runs, and visualise the resulting atomic structures directly in their browser. The completion of this cycle will transform the tool from a functional CLI utility into a comprehensive, user-friendly, and scientifically powerful platform for automated MLIP dataset generation, fulfilling the project's ultimate vision.

## 2. System Architecture

Cycle 2 significantly expands the project's capabilities, evolving the stubs from Cycle 1 into fully-featured components and introducing a new web application layer. The file structure will be modified to reflect these new additions, ensuring the project remains clean and maintainable.

**File Structure for Cycle 2:**

```
src/mlip_autopipec/
├── __init__.py
├── cli.py
├── config.py              # To be modified. Add ExplorationConfig and SamplingConfig.
├── database.py
├── factories.py           # To be modified. Add create_explorer and create_sampler.
├── interfaces.py          # To be modified. Define IExplorer and ISampler.
├── pipeline/
│   ├── __init__.py
│   ├── orchestrator.py
│   ├── generation.py
│   ├── exploration.py     # To be modified from stub. Implement the full MDExplorer.
│   └── sampling.py        # To be modified from stub. Implement FarthestPointSampler.
└── web/
    ├── __init__.py
    ├── app.py             # To be created. Main web application logic (Streamlit).
    └── templates/         # To be created (if needed by web framework).
        └── index.html
tests/
├── __init__.py
├── conftest.py
├── test_config.py
├── test_database.py
├── test_generation.py
├── test_exploration.py    # To be created. Unit tests for MDExplorer.
└── test_sampling.py       # To be created. Unit tests for FarthestPointSampler.
```

**Key Component Blueprints:**

*   **`config.py`**: The `FullConfig` model will be extended to include new Pydantic models: `ExplorationConfig` and `SamplingConfig`. These new models will be nested within `FullConfig` to maintain a clean, hierarchical configuration structure.
    *   `ExplorationConfig` will define a rich set of parameters for controlling the MD/MC simulation, such as `temperature_k` (with a list-of-floats option for simulated annealing), `pressure_gpa`, `timestep_fs`, `mc_swap_probability`, and the name of the MLIP model to be used. It will include validators to ensure that temperature and timestep are positive.
    *   `SamplingConfig` will define the sampling method (using a `Literal['random', 'fps']` for type safety) and its corresponding parameters, like the `num_samples_to_select`, ensuring it's a positive integer.

*   **`pipeline/exploration.py`**: This module will be completely overhauled from its Cycle 1 stub. It will house the `MDExplorer` class, which will implement the `IExplorer` interface. This class is the scientific core of the application and will be responsible for:
    *   Taking a list of `ase.Atoms` objects as input.
    *   Setting up and running MD simulations in parallel using Python's `concurrent.futures.ProcessPoolExecutor`.
    *   Implementing the "late-binding" of the MLIP calculator. This is a critical performance optimization where the large PyTorch model is loaded from disk *inside* the worker process, avoiding the extremely slow and memory-intensive process of pickling it from the main process.
    *   Implementing a `detect_vacuum` helper function that analyzes the cell geometry to determine if it is a bulk or surface system.
    *   Dynamically selecting the correct thermodynamic ensemble (NPT for bulk, NVT for surfaces) based on the output of `detect_vacuum`.
    *   Integrating a ZBL (Ziegler-Biersack-Littmark) potential for robust handling of short-range interactions, which prevents atoms from getting unrealistically close and causing the simulation to crash, especially at high temperatures.
    *   Executing hybrid MD/MC moves by periodically attempting to swap atoms of different species, with the acceptance probability determined by the Metropolis criterion based on the change in potential energy.

*   **`pipeline/sampling.py`**: The basic `RandomSampler` will be supplemented by a new `FarthestPointSampler` class, which will implement the `ISampler` interface. This class will perform the following steps:
    *   Take the raw simulation trajectories (a large list of `ase.Atoms` objects) as input.
    *   Use a library like `dscribe` to calculate high-dimensional SOAP (Smooth Overlap of Atomic Positions) descriptors for each structure. The SOAP vector serves as a "fingerprint" of the local atomic environment.
    *   Implement the iterative FPS algorithm. This starts by selecting a random structure, then iteratively finds the structure in the remaining set that is farthest away (in Euclidean distance in the SOAP feature space) from the set of already-selected structures, adding it to the selection until the desired number of samples is reached.

*   **`web/app.py`**: This new module will contain the web application logic. We will use **Streamlit** as the framework due to its rapid development cycle and simplicity for data-centric applications. The `app.py` script will:
    *   Define the UI layout, including input widgets (sliders for temperature, text boxes for elements, dropdowns for samplers) that directly map to the parameters in the `FullConfig` Pydantic model.
    *   Provide a "Generate Configuration" button that displays the corresponding YAML for the user's review.
    *   Include a "Run Pipeline" button. When clicked, this will save the configuration to a temporary file and launch the `mlip-autopipec run` CLI command as a background subprocess using `subprocess.Popen`.
    *   Display the real-time standard output and standard error from the subprocess in a scrolling text box on the page, providing live feedback.
    *   Offer a results visualisation component. After a run completes, the user can select the output database, and the app will use `ase.db.connect` to read the structures and a library like `py3Dmol` to render them interactively in the browser.

*   **`factories.py`**: The factory functions will be updated to handle the new components. New `create_explorer` and `create_sampler` functions will be added. These will read the relevant sections of the configuration and instantiate the correct class (`MDExplorer` or `FarthestPointSampler` vs. their simpler counterparts), injecting them as dependencies into the main orchestrator.

## 3. Design Architecture

The advanced components of Cycle 2 are designed for performance, robustness, and scientific accuracy. The Pydantic schemas are extended to provide fine-grained control over this new functionality, and the design patterns established in Cycle 1 are maintained to ensure consistency and maintainability.

**Pydantic Schema Design:**

*   **`ExplorationConfig(BaseModel)`**: This model provides detailed control over the simulation.
    *   `model_name: str`: The name of the MLIP model to use (e.g., 'MACE_v1'). The factory will use this to dynamically load the correct calculator.
    *   `temperature_k: float = Field(..., gt=0)`: Simulation temperature in Kelvin.
    *   `pressure_gpa: Optional[float] = None`: Simulation pressure in GPa. The presence or absence of this value is the primary trigger for NPT vs. NVT, although this can be overridden by the automatic vacuum detection.
    *   `num_steps: int = Field(..., gt=0)`: Total number of MD steps.
    *   `timestep_fs: float = Field(..., gt=0)`: Timestep for the MD integrator in femtoseconds.
    *   `use_zbl_repulsion: bool = True`: Flag to enable the ZBL short-range potential, defaulting to on for safety.
    *   `mc_swap_probability: float = Field(0.0, ge=0, le=1)`: Probability of attempting a Monte Carlo atom swap at each MD step. A value of 0.0 disables MC moves entirely.

*   **`SamplingConfig(BaseModel)`**: This model controls the data selection process.
    *   `method: Literal['random', 'fps'] = 'fps'`: The chosen sampling algorithm, defaulting to the more advanced method. The use of `Literal` provides compile-time safety and automatic validation.
    *   `num_samples: int = Field(..., gt=0)`: The target number of structures for the final dataset. A validator will be included to ensure this is not greater than the total number of structures produced by the exploration stage.

**Data Flow and Consumers/Producers:**

The core data flow remains a linear pipeline, but the `exploration` and `sampling` stages are now active and transformative processors of the data.

1.  **`MDExplorer` (Consumer/Producer)**: The explorer consumes the list of seed `ase.Atoms` objects produced by the generation stage. After running the MD/MC simulations, it acts as a producer, outputting a much larger list of `ase.Atoms` objects. This list represents the full, time-resolved simulation trajectories for all the initial seed structures. This potentially very large dataset is then passed in memory (or via temporary files if too large) to the sampling stage.

2.  **`FarthestPointSampler` (Consumer/Producer)**: The sampler consumes the large list of trajectory structures. It performs its descriptor calculation and iterative selection process, and then acts as a producer, outputting a smaller, highly curated list of `ase.Atoms` objects. This final list, which represents the most diverse and informative snapshots from the simulations, is then passed to the `AseDBWrapper` for final storage.

3.  **`web/app.py` (Producer/Consumer)**: The web application has a dual role in the data flow.
    *   **As a Producer**: The web UI acts as an interactive producer of the `FullConfig` object. Instead of the user writing YAML manually, the UI's form widgets generate the configuration data. When the user clicks "Run", the web app serializes this configuration into a temporary YAML file, which is then consumed by the CLI backend.
    *   **As a Consumer**: After a pipeline run is complete, the web UI acts as a consumer of the final output database (`.db` file). It reads the `ase.Atoms` objects from the database to present them in the interactive visualisation view, completing the feedback loop for the user.

**Key Design Considerations:**

*   **Performance and Memory Management**: The parallel execution of MD simulations is critical. The design must ensure that the `ProcessPoolExecutor` is used effectively. A key consideration is managing the data transfer between the main process and the worker processes. For this reason, workers will write their output trajectories to individual temporary files on disk. The main process will then collect and concatenate these files before passing them to the sampler, avoiding having to pass gigabytes of atomic coordinate data through inter-process queues.
*   **Robustness and Error Handling**: The MD simulation is the most likely point of failure. The `MDExplorer` must be wrapped in robust `try...except` blocks. If a single simulation worker crashes (e.g., due to a "Coulomb explosion"), the main process should log the error, save the failed input structure for debugging, and continue processing the results from the other, successful simulations. The application should not fail entirely just because one of its parallel tasks did.
*   **Decoupling the UI from the Backend**: The Web UI in `web/app.py` will be strictly decoupled from the core pipeline logic. It will trigger the pipeline via a system call to the `cli.py` entry point (`subprocess.Popen(['mlip-autopipec', 'run', ...])`). This is a simple yet powerful design choice. It ensures the UI and the backend logic remain independent; we could completely rewrite the backend in another language, and as long as it exposes the same CLI, the UI would continue to work. It also means the core team can continue to work on and test the backend without needing to run the UI, and vice versa.

## 4. Implementation Approach

The implementation of Cycle 2 will be layered, starting with the core backend logic for exploration and sampling, and then building the user-facing web interface on top of that stable foundation. A test-driven development approach will be used throughout.

1.  **Extend Configuration**: First, modify `config.py` to add the new `ExplorationConfig` and `SamplingConfig` Pydantic models. Update the `FullConfig` model to include them. Add unit tests for the new validation rules in `test_config.py`.

2.  **Implement `MDExplorer`**: This is the most complex task of the cycle. Replace the `PassThroughExplorer` in `pipeline/exploration.py` with the `MDExplorer`.
    *   Start by implementing the basic MD loop for a single structure using an ASE dynamics object (e.g., `Langevin`).
    *   Integrate the "late-binding" calculator logic, where the MLIP model is loaded only within the worker process function.
    *   Refactor the single-simulation logic into a top-level function that can be targeted by `ProcessPoolExecutor`.
    *   Implement the main `explore` method, which will use the `ProcessPoolExecutor` to map the input structures to the simulation worker function.
    *   Implement the automatic ensemble switching by adding the vacuum detection logic and using an `if/else` block within the worker function to select the appropriate ASE dynamics class (`NPT` or `NVT`).
    *   Add support for hybrid MD/MC by creating a custom ASE dynamics move that attempts an atom swap and uses the Metropolis criterion for acceptance, and attach it to the dynamics object.

3.  **Implement `FarthestPointSampler`**: Replace the `RandomSampler` in `pipeline/sampling.py` with the `FarthestPointSampler`.
    *   Add `dscribe` as a project dependency.
    *   In the `sample` method, first set up the SOAP descriptor calculator from `dscribe`.
    *   Map the input structures to the descriptor calculator to get a list of SOAP vectors.
    *   Implement the FPS algorithm itself. This will be a loop that iteratively selects the structure whose SOAP vector has the greatest minimum distance to the SOAP vectors of the already-selected structures.

4.  **Update Factories and Interfaces**: Modify `interfaces.py` to define the new `IExplorer` and `ISampler` abstract base classes. Update the factory functions in `factories.py` to instantiate `MDExplorer` and `FarthestPointSampler` based on the settings in the new configuration models.

5.  **Write Backend Tests**: Create `tests/test_exploration.py` and `tests/test_sampling.py`. Write focused unit tests for the new components, using mocks for the computationally expensive parts. For example, mock the `ase.calculators` object and the `dscribe.SOAP` object to test the setup and logic without the computational overhead.

6.  **Develop Web UI**: Create the `web/app.py` module using Streamlit.
    *   Add `streamlit` as a project dependency.
    *   Build the input form using `st.sidebar` for layout. Use widgets like `st.slider` for temperature, `st.text_input` for elements, and `st.selectbox` for samplers, mapping them to the `FullConfig` Pydantic model.
    *   Implement the "Run" button callback. This function will gather the state from the UI widgets, construct the `FullConfig` object, serialize it to a temporary YAML file, and then launch the CLI command using `subprocess.Popen`, redirecting stdout and stderr.
    *   Implement a log display area using `st.empty()` and a loop that reads the output from the subprocess in real-time.
    *   Add a results page that uses `st.file_uploader` to allow the user to select an output database. When a database is loaded, use `ase.db.connect` to read the structures and `py3Dmol` to render them in an interactive component.

7.  **Final Integration Testing**: Perform end-to-end testing of the entire, updated pipeline. This will involve creating a new integration test that uses a configuration enabling the MD explorer and FPS sampler. This test will use a fast classical potential (like EMT) to make it feasible to run quickly in a CI environment. Finally, manually run through the Web UI scenarios to ensure a smooth user experience.

## 5. Test Strategy

The testing for Cycle 2 must cover the new, complex backend logic as well as the new user interface. This requires expanding our testing pyramid to include not just unit and integration tests, but also end-to-end browser-based testing.

**Unit Testing Approach (Min 300 words):**

The unit testing for Cycle 2 will focus heavily on the new scientific components, ensuring their internal logic is correct in isolation. For the `MDExplorer`, running a full MD simulation is not feasible for a unit test. Therefore, our tests will be designed to verify the *setup* and *control flow* of the simulation logic without the computational cost. We will use `pytest-mock` to patch the `ase.calculators` object and the ASE dynamics objects (`NPT`, `NVT`). One key test will verify the automatic ensemble switching. We will create a bulk `ase.Atoms` object and pass it to the explorer, then assert that the `ase.md.NPT` class was instantiated. We will then repeat the test with a surface slab structure and assert that `ase.md.NVT` was instantiated instead. This confirms the vacuum detection and switching logic is working. We will also test the MC logic. We will mock the MLIP calculator so that we can control the energy it reports. We will run a single step of the dynamics, trigger a swap move, and mock the calculator to report a lower energy for the swapped configuration. We will then assert that the atoms' positions have indeed been swapped. We will repeat the test, mocking a higher energy, and assert that the positions have not changed, thus verifying our implementation of the Metropolis criterion.

For the `FarthestPointSampler`, we need to test the correctness of the FPS algorithm. A real SOAP calculation is too slow for a unit test, so we will mock the `dscribe.SOAP` object. The test will focus on the selection algorithm itself. We will create a small, deterministic set of input "structures" (which can be simple mock objects) and provide a fixed list of corresponding feature vectors (e.g., simple NumPy arrays) that will be returned by the mocked SOAP calculator. The test will then call the `sample` method and assert that the indices of the structures selected are exactly the ones we would expect from a manual, on-paper calculation of the FPS algorithm. For example, if we have points at (0,0), (1,0), (5,0), and (6,0) and ask for 3 samples, we expect the test to select the points at (0,0), (6,0), and then (1,0) or (5,0), demonstrating that it correctly prioritizes the most distant points. We will also test edge cases, such as requesting more samples than are available in the input, ensuring the sampler handles this gracefully by returning all input structures.

**Integration Testing Approach (Min 300 words):**

The integration tests for Cycle 2 will verify the full pipeline, now including the real exploration and sampling stages. As a full pipeline run with a real MLIP model would be far too slow and resource-intensive for an automated test suite, we will employ a crucial strategy: replacing the slow MLIP model with a very fast, classical calculator that is built into ASE, such as the Effective Medium Theory (`EMT`) potential. This allows us to run a genuine, short MD simulation that completes in seconds, exercising the entire code path of the exploration engine.

A full backend integration test will be structured as follows:
1.  **Arrange**: The test will create a `config.yaml` file in an isolated filesystem. This configuration will specify the `EMT` calculator, a short MD run (e.g., 100 steps), `fps` sampling to select 5 final structures, and a target output database.
2.  **Act**: The test will invoke the CLI runner with this configuration.
3.  **Assert**: After the run completes successfully, the test will perform several key verifications. It will connect to the output database and assert that it contains exactly 5 structures, as requested by the FPS sampler. More importantly, it will verify that the exploration stage actually worked. It will do this by reading the first structure from the database and comparing it to the initial seed structure that was generated. It will assert that their atomic positions are *not* identical and that their potential energies are different. This provides strong, direct evidence that the structures were genuinely evolved by the MD simulation and not just passed through, confirming that the entire exploration and sampling pipeline is working correctly end-to-end.

For the Web UI, we will introduce a new layer of End-to-End (E2E) testing using a browser automation framework like **Playwright**. The test script for the UI will be comprehensive:
1.  **Setup**: The test will start the Streamlit web application as a separate subprocess.
2.  **Act**: Playwright will launch a real browser, navigate to the Streamlit app's local URL, and then programmatically interact with the page. It will locate the input widgets by their labels and fill them in (e.g., `'st.text_input(label="Elements")'` will be filled with "Si"). After filling the form, it will simulate a user clicking the "Run Pipeline" button.
3.  **Assert**: The test will then wait and observe the state of the page. It will assert that the log viewer appears and that its text content is updated with the expected log messages from the pipeline. Finally, it will wait for a "Pipeline completed successfully" message to appear. This provides the ultimate confidence that the application is working correctly from the user's perspective, from browser interaction all the way to backend computation and back to the UI.
