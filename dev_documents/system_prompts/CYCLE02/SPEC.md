# CYCLE02 Specification: Advanced Features and Web UI

## 1. Summary

This document details the technical specifications for Cycle 2 of the MLIP-AutoPipe project. Building upon the stable command-line interface (CLI) and core pipeline established in Cycle 1, this cycle introduces a significant expansion of the system's capabilities and accessibility. The primary goals are twofold: first, to enhance the scientific sophistication of the backend engine by implementing advanced simulation and sampling techniques; and second, to develop a user-friendly Web-based Graphical User Interface (Web UI) that lowers the barrier to entry and provides an intuitive, interactive alternative to the CLI.

The scientific enhancements will focus on the `Exploration` and `Sampling` modules. The `MDEngine` will be upgraded to a hybrid MD/MC (Monte Carlo) engine, capable of performing atom swaps and other non-standard moves that are critical for efficiently exploring the compositional and configurational space of materials like alloys. It will also incorporate automatic thermodynamic ensemble switching and the ZBL potential for more realistic high-energy collision handling. The `Sampling` module will be upgraded with a Farthest Point Sampling (FPS) implementation, allowing for intelligent, diversity-driven selection of structures, which is vastly superior to the random sampling of Cycle 1.

The second major focus is the creation of a Web UI. This interface will allow users to interactively define a material system, configure simulation parameters, launch a pipeline run as a background process, monitor its progress, and visualize the resulting structures directly in their browser. This feature is crucial for making the MLIP-AutoPipe tool accessible to a wider audience, including students and researchers who may not be comfortable with a purely command-line-based workflow. The successful completion of this cycle will transform the MLIP-AutoPipe from a functional tool for experts into a comprehensive and accessible platform for automated MLIP dataset generation.

## 2. System Architecture

The architecture in Cycle 2 expands upon the existing structure. The core backend logic remains, but it is enhanced with new capabilities. A new, distinct entry point, the `web_ui.py`, is added, which interacts with the established pipeline.

**File Structure:**
Bold entries indicate new or significantly modified files for this cycle.

```
.
├── pyproject.toml
└── src/
    └── mlip_autopipec/
        ├── __init__.py
        ├── cli.py
        ├── **web_ui.py**                 # New entry point for the Web UI (Streamlit)
        ├── pipeline/
        │   ├── __init__.py
        │   └── runner.py
        ├── generators/
        │   ├── __init__.py
        │   ├── base.py
        │   ├── alloy.py
        │   └── factory.py
        ├── explorers/
        │   ├── __init__.py
        │   └── **md_engine.py**          # Modified to include hybrid MD/MC, ZBL, etc.
        ├── sampling/
        │   ├── __init__.py
        │   ├── base.py
        │   ├── random_sampler.py
        │   └── **fps_sampler.py**        # New sampler using Farthest Point Sampling
        ├── storage/
        │   ├── __init__.py
        │   └── ase_db_writer.py
        └── common/
            ├── __init__.py
            ├── atoms_validator.py
            └── exceptions.py
```

**Code Blueprints:**

*   **`web_ui.py`**: This file will be built using the `streamlit` library.
    *   It will contain UI elements (e.g., `st.text_input`, `st.slider`, `st.selectbox`) for the user to define all the necessary configuration parameters.
    *   A "Run Pipeline" button (`st.button`) will trigger the main logic.
    *   On button click, the script will:
        1.  Assemble the user inputs into a dictionary matching the Hydra config structure.
        2.  Write this dictionary to a temporary YAML file.
        3.  Use Python's `subprocess.Popen` to launch the CLI command (`mlip-autopipec run ...`) as a background process, redirecting its stdout and stderr to a log file.
        4.  Display the contents of the log file in the UI in real-time to show progress.
        5.  Once the process completes, it will provide a link to download the database or attempt to visualize the results directly using a library like `py3Dmol`.
*   **`explorers/md_engine.py`**: The `MDEngine` class will be substantially upgraded.
    *   The `run` method will be modified to accept more complex exploration configurations.
    *   It will include logic to perform Monte Carlo moves (e.g., `atom_swap`) periodically during the MD run.
    *   A `_setup_calculator` method will be added to handle the mixing of an MLIP calculator with the ZBL potential for accurate short-range repulsions.
    *   A `_detect_ensemble` method will be implemented to automatically choose between NVT and NPT ensembles based on whether a vacuum slab is detected in the `ase.Atoms` object.
*   **`sampling/fps_sampler.py`**:
    *   `FPSSampler(BaseSampler)`: A new class that implements the `sample` method.
    *   It will first read the entire trajectory.
    *   It will then use a library like `scikit-learn` or a custom implementation to compute SOAP (Smooth Overlap of Atomic Positions) descriptors for each structure.
    *   Finally, it will run the Farthest Point Sampling algorithm on these descriptors to select a structurally diverse subset of atoms.

## 3. Design Architecture

The design in Cycle 2 focuses on managing increased complexity in the backend and creating a clean separation between the UI and the core application logic.

*   **Configuration Schema Expansion (Conceptual Pydantic Models):**
    The configuration models from Cycle 1 will be expanded to include the new features.

    ```python
    # In src/mlip_autopipec/config_models.py (updated)

    class MCCarbonCopy(BaseModel):
        # Configuration for Monte Carlo moves
        enabled: bool = False
        swap_frequency: int = 10
        # ... other MC parameters

    class ExplorationConfig(BaseModel):
        # ... existing fields
        temperature_k: float = Field(gt=0)
        num_steps: int = Field(gt=0)
        calculator: str
        ensemble_switching: bool = True
        use_zbl_potential: bool = True
        mc_moves: MCCarbonCopy = Field(default_factory=MCCarbonCopy)

    class SamplingConfig(BaseModel):
        method: str = "Random" # Can now be "FPS"
        num_samples: int = Field(gt=0)
        # ... FPS specific parameters might be added here

    # ... MainConfig remains largely the same, but now accepts
    # the expanded ExplorationConfig and SamplingConfig
    ```

*   **UI-Backend Interaction:**
    The key design choice here is to **decouple the Web UI from the core pipeline**. The Web UI will not import and run the `PipelineRunner` directly. Instead, it will act as a "configurator" and "process launcher".
    1.  **Producer (Web UI):** The user's interactions with sliders and forms in the browser produce a state. When the "Run" button is clicked, this state is serialized into a YAML configuration file.
    2.  **Consumer (CLI):** The `subprocess.Popen` call executes the standard CLI entry point. This consumer is the exact same, well-tested application from Cycle 1. It knows nothing about the Web UI; it simply consumes the configuration file it is given.
    3.  **Data Flow:** The flow is `Web UI State -> YAML Config File -> CLI Process -> Log File & Database`. The Web UI then reads the `Log File` and `Database` to display the results.
    This decoupled design is robust. The core logic does not need to be modified to support the UI, and the UI can be developed and tested independently. It also ensures that a workflow is reproducible whether it is initiated from the UI or manually via the CLI.

*   **Invariants and Constraints:**
    *   **Process Isolation:** The Web UI process is separate from the computational pipeline process. This is a critical constraint. If the user closes their browser tab, the `subprocess` will continue running on the server, preventing data loss.
    *   **Statelessness (UI):** The Streamlit app should be largely stateless. All the information needed to run a calculation is captured in the generated config file. This makes the system more predictable.
    *   **Scientific Validity:** The new features in `MDEngine` must maintain the physical invariants established in Cycle 1. The ZBL potential prevents atoms from getting too close, reinforcing the `overlap_check`. The atom swap MC moves must be designed to conserve the overall stoichiometry of the system.

## 4. Implementation Approach

The implementation will be phased, focusing first on upgrading the backend capabilities and then on building the user-facing Web UI.

1.  **Upgrade Dependencies:** Add `streamlit`, `scikit-learn` (or a lighter SOAP library), and potentially `py3Dmol` to `pyproject.toml` and run `uv sync`.
2.  **Implement `FPSSampler`:**
    *   Create the `fps_sampler.py` file.
    *   Implement the `FPSSampler` class. This is a scientifically complex task and should be tackled first. It will involve finding a suitable library for SOAP descriptor calculation and implementing the iterative FPS algorithm.
    *   Write unit tests to verify its behavior with a small, predictable set of vectors.
3.  **Enhance `MDEngine`:**
    *   Modify the `explorers/md_engine.py` file.
    *   Add the logic for Monte Carlo moves. This can be a simple `_perform_swap_move` method that is called every `N` steps within the main simulation loop.
    *   Implement the logic for mixing the MLIP calculator with ASE's `ZBL` potential. This usually involves creating a wrapper class or function.
    *   Implement the `_detect_ensemble` logic. This can be a simple check on the `atoms.pbc` (Periodic Boundary Conditions) array.
4.  **Integrate Backend Features:**
    *   Update the `PipelineRunner` and the Pydantic configuration models to accept and pass through the new configuration options for sampling (`method: "FPS"`) and exploration (MC moves, etc.).
    *   Extend the CLI integration test from Cycle 1 to run a small case using the new `FPS` sampler and one of the new MD features to ensure they are correctly plumbed in.
5.  **Develop the Web UI (`web_ui.py`):**
    *   Scaffold the basic Streamlit application layout. Use `st.sidebar` for configuration options.
    *   Add widgets for all key parameters: element selection, composition, temperature, etc.
    *   Write the function that is called when the "Run" button is pressed. This function will collect the state from all widgets, build a dictionary, and save it as a YAML file.
    *   Implement the `subprocess.Popen` logic to execute the `mlip-autopipec run ...` command.
    *   Implement the real-time log display by reading the stdout of the subprocess.
    *   After the process finishes, add a button to download the resulting database file (`st.download_button`).
    *   (Optional Stretch Goal) Implement a simple results visualizer by reading the first few structures from the output database and displaying them with `py3Dmol`.

## 5. Test Strategy

Testing in Cycle 2 must cover both the new backend logic and the new Web UI, requiring different tools and approaches for each.

**Unit Testing Approach (Min 300 words):**
*   **`FPSSampler`**: This is the most critical new component to unit test. A test will be created with a small, pre-defined trajectory file containing a few simple, non-random structures (e.g., a compressed cell, an expanded cell, a sheared cell). The SOAP descriptors for these will be calculated. We will assert that the `FPSSampler`, when asked to pick two samples, correctly chooses the most dissimilar ones (e.g., the compressed and expanded cells). The test should be deterministic.
*   **`MDEngine` Enhancements**: The new logic within the engine will be tested in isolation.
    *   A test for the MC swap move will be created. It will take a two-element `Atoms` object (e.g., one Fe, one Pt), run the swap function, and assert that the positions of the two atoms have been exchanged.
    *   A test for ensemble detection will be created. We will pass an `Atoms` object with `pbc=[True, True, True]` and assert the engine selects the NPT ensemble. We will then pass an object with `pbc=[True, True, False]` and assert it selects NVT.
*   **Web UI Backend Logic**: While much of the UI is visual, the logic that generates the configuration file can be unit tested. A function `create_config_from_state(ui_state: dict) -> dict` can be created and tested. We will provide a sample `ui_state` dictionary and assert that the output dictionary has the correct structure and values, ready to be written to YAML. This tests the data transformation logic without needing to run the full UI.

**Integration Testing Approach (Min 300 words):**
*   **Extended CLI Test**: The integration test from Cycle 1 will be duplicated and modified. The new test will use a configuration file that explicitly enables the advanced features: `sampling.method: "FPS"` and `exploration.mc_moves.enabled: True`. The test will run the full pipeline via the `CliRunner` and perform the same assertions as before (database exists, has the correct number of atoms). This verifies that the new components are correctly integrated into the main `PipelineRunner`.
*   **End-to-End UI Testing**: The Web UI requires a different approach. We will use a browser automation tool like `Playwright`. The test script will:
    1.  Start the Streamlit application as a separate process.
    2.  Use Playwright to open a browser and navigate to the local Streamlit URL.
    3.  Programmatically interact with the UI elements: find the text input for elements and type "Fe, Pt", find the slider for temperature and set its value.
    4.  Click the "Run Pipeline" button.
    5.  Wait for the log output in the UI to contain the "Pipeline completed successfully" message. This confirms the backend process was triggered and completed.
    6.  After the run, the test script will locate the generated database file on disk and, just like the CLI integration test, connect to it and verify its contents. This provides a true end-to-end test of the entire Web UI workflow, from user interaction to final data storage.
