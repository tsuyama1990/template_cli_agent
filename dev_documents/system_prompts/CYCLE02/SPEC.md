# SPECIFICATION: Cycle 2 - Advanced Exploration Engine and Web UI

## 1. Summary

This document provides the detailed technical specification for the second and final development cycle of the MLIP-AutoPipe project. This cycle is designed to build upon the robust command-line foundation established in Cycle 1, introducing the advanced, scientifically sophisticated features that truly define the project's innovative core. Furthermore, it will deliver a user-friendly web interface to make these powerful capabilities accessible to a much broader audience. The centerpiece of this cycle's technical development is the implementation of the hybrid Molecular Dynamics/Monte Carlo (MD/MC) exploration engine. This is not merely an incremental improvement; it is the core intellectual property of the project, designed to facilitate a more intelligent and efficient exploration of the potential energy surface than what is possible with standard simulation techniques. This advanced engine will integrate real, high-fidelity MLIP models (such as MACE) to ensure the accuracy of energy and force calculations. Crucially, it will also co-integrate a classical ZBL potential to correctly model the strong, short-range repulsive interactions between atoms. This mixing of potentials is a critical feature for ensuring the stability of high-temperature simulations, preventing the catastrophic "Coulomb explosion" failures that can plague simulations of dense or highly energetic systems.

The second major pillar of this cycle is the development of a comprehensive, web-based Graphical User Interface (GUI). While the CLI from Cycle 1 provides essential functionality for expert users and automated, large-scale batch processing, the GUI will democratize access to the tool. It will provide an intuitive, interactive environment where users can visually construct their run configurations, launch complex pipeline runs with the click of a button, monitor the progress of those runs in real-time, and, most importantly, visualize the resulting atomic structures. This makes the tool suitable for educational purposes, for experimental scientists who may not be comfortable with the command line, and for rapid prototyping of ideas. In parallel, this cycle will also see the implementation of the Farthest Point Sampling (FPS) algorithm. This replaces the basic random sampler from Cycle 1 with a more intelligent, diversity-driven method for selecting the most informative structures from the vast exploration trajectories. Upon the successful completion of Cycle 2, the MLIP-AutoPipe project will be a feature-complete, scientifically robust, and highly user-friendly computational framework, fully realizing the ambitious vision laid out in the initial architectural plans.

## 2. System Architecture

The architecture in Cycle 2 is an extension and refinement of the foundation laid in Cycle 1. It involves a significant upgrade of existing modules with more complex scientific logic and the introduction of a completely new set of modules to handle the web-based user interface. The core principles of modularity and separation of concerns remain paramount.

**File Structure (Cycle 2 Focus):**

The files and directories to be created or modified in this cycle are marked in **bold**. It is important to note that many files from Cycle 1 (e.g., `md_engine.py`, `config.py`) will be significantly *modified* and enhanced to incorporate the new, advanced logic.

```
.
├── src/
│   └── mlip_autopipec/
│       ├── **main_gui.py**         # The primary entry point for launching the web UI
│       ├── explorers/
│       │   └── **md_engine.py**    # Will be significantly upgraded to the hybrid MD/MC engine
│       ├── samplers/
│       │   ├── **fps.py**          # New file: FarthestPointSampler implementation
│       │   └── **base.py**         # Modified to accommodate the new sampler type
│       ├── common/
│       │   └── **config.py**       # Pydantic models will be updated with many new options
│       └── ui/
│           ├── **__init__.py**
│           ├── **app.py**          # The main Flask/Streamlit application logic
│           ├── **static/**         # Directory for CSS, JS assets
│           └── **templates/**      # Directory for HTML templates (if using Flask)
└── tests/
    ├── **__init__.py**
    ├── **test_explorer.py**      # New test file for the advanced hybrid engine
    ├── **test_sampler.py**       # New test file for the FPS sampler
    └── **test_ui.py**            # New file for E2E tests of the web interface
```

**Code Blueprints:**

These blueprints provide a detailed view of the key classes and their new or modified structures.

*   **`explorers/md_engine.py` - `MDExplorer` (Significantly Upgraded):**
    This class will evolve from a simple placeholder into the complex, intelligent heart of the application.
    ```python
    from mlip_autopipec.common.utils import detect_vacuum # A new utility

    class MDExplorer:
        def __init__(self, config: ExplorationConfig):
            self.config = config
            # ZBL potential for short-range repulsion is initialized once
            self._zbl_potential = self._init_zbl()
            # The main MLIP calculator is now initialized lazily within each worker
            # process to ensure multiprocessing safety and efficiency.

        def _get_calculator(self):
            # Implements the "late-binding" calculator pattern.
            # This method is called inside the worker process.
            mlip_model = load_mlip_model(self.config.mlip_model_path)
            # The MixedCalculator combines the MLIP with ZBL for a robust potential
            return MixedCalculator(mlip_model, self._zbl_potential)

        def _perform_mc_swap(self, atoms: Atoms):
            # Detailed logic for attempting a Monte Carlo atomic swap move.
            # 1. Select two different atoms to swap.
            # 2. Check charge safety constraints for ionic systems.
            # 3. Calculate energy before and after the swap.
            # 4. Accept or reject the move based on the Metropolis criterion.
            pass

        def _run_single_md(self, atoms: Atoms) -> list[Atoms]:
            calculator = self._get_calculator()
            atoms.set_calculator(calculator)

            # Automatically detect if the system has a vacuum slab
            is_surface = detect_vacuum(atoms)
            # Choose the correct statistical ensemble based on the detection
            if is_surface:
                ensemble = NVTEnsemble(atoms, ...)
            else:
                ensemble = NPTEnsemble(atoms, ...)

            # The main MD loop now interleaves MD steps with MC moves
            trajectory = []
            for step in range(self.config.md_steps):
                ensemble.run(self.config.md_step_interval)
                if self.config.use_mc_swaps and step % self.config.mc_swap_frequency == 0:
                    self._perform_mc_swap(atoms)
                trajectory.append(atoms.copy())
            return trajectory

        def run(self, initial_structures: list[Atoms]) -> list[list[Atoms]]:
            # This method will use ProcessPoolExecutor to run _run_single_md
            # in parallel for each of the initial structures.
            # It will handle task distribution and result collection.
            pass
    ```

*   **`samplers/fps.py` - `FarthestPointSampler`:**
    This new class provides the intelligent sampling logic.
    ```python
    from .base import BaseSampler
    # This module will have new dependencies: scikit-learn for PCA/SOAP,
    # and potentially a dedicated library for SOAP descriptors.

    class FarthestPointSampler(BaseSampler):
        def sample(self, trajectories: list[list[Atoms]]) -> list[Atoms]:
            # 1. Flatten all trajectories into a single, large list of structures
            all_structures = self._flatten(trajectories)

            # 2. Compute SOAP descriptors for all structures in parallel.
            # This is a computationally intensive step.
            soap_descriptors = self._compute_soap_in_parallel(all_structures)

            # 3. Perform the iterative Farthest Point Sampling algorithm on the
            # high-dimensional descriptor vectors to find the most diverse set.
            # This involves repeatedly finding the point that is farthest
            # from the current set of selected points.
            selected_indices = self._perform_fps(soap_descriptors, num_to_select=self.config.num_samples)

            # 4. Return the Atoms objects corresponding to the selected indices.
            return [all_structures[i] for i in selected_indices]
    ```

*   **`main_gui.py` and `ui/app.py`:**
    These files will contain the logic for the web application, likely using Streamlit for its rapid development capabilities.
    ```python
    import streamlit as st
    from mlip_autopipec.pipeline.runner import PipelineRunner
    # ... other necessary imports

    def build_config_from_ui_state():
        # This function will gather all the values from the Streamlit widgets
        # (sliders, text inputs, etc.) and construct a validated MainConfig
        # Pydantic object. It will handle type conversions and error checking.
        pass

    def main():
        st.title("MLIP-AutoPipe: Automated Dataset Generation")

        # The configuration UI will be in a persistent sidebar
        with st.sidebar:
            st.header("1. System Configuration")
            # Create UI widgets for all config options found in the Pydantic models
            # e.g., st.multiselect for elements, st.number_input for temperature

        # Main application area
        st.header("2. Run Pipeline")
        if st.button("Start Generation Pipeline"):
            try:
                # 1. Build and validate the config from the UI
                config = build_config_from_ui_state()

                # 2. Instantiate the runner and run in a separate thread
                # to avoid blocking the UI. Update UI with progress.
                st.info("Pipeline started...")
                runner = PipelineRunner(config)
                # This would need a more complex setup with threading and callbacks
                # to provide real-time feedback to the UI.
                runner.run()
                st.success("Pipeline finished successfully!")

            except Exception as e:
                st.error(f"An error occurred: {e}")

        st.header("3. View Results")
        # Add UI elements to select and load a .db file, and then
        # use ASE's plotting capabilities to visualize the structures.
    ```

## 3. Design Architecture

While the core four-stage pipeline architecture from Cycle 1 remains intact, Cycle 2 extends it with significantly more sophisticated internal logic and an entirely new user-facing layer for the UI. The Pydantic schemas defined in `common/config.py` continue to be the single source of truth for all configurable options, but they will be substantially expanded to include parameters for the new, advanced features. For example, the `ExplorationConfig` model will gain new fields such as `mlip_model_path: str`, `use_mc_swaps: bool`, `mc_swap_frequency: int`, and `auto_ensemble_switching: bool`. Similarly, the `SamplingConfig` model will be updated to include a `method: Literal["random", "fps"]` field, allowing the user to choose their desired sampling strategy. This schema-driven approach ensures that even as the system's complexity grows, its configuration remains robust, well-documented, and easy to validate.

The Web UI will be designed as a cleanly separated presentation layer that interacts with the core pipeline but does not contain any of the scientific logic itself. This is a critical design choice for maintainability. The UI's responsibilities are strictly limited to the following three tasks:
1.  **Configuration:** It will programmatically inspect the Pydantic schemas to dynamically generate a user-friendly web form with appropriate widgets (sliders for temperature, dropdowns for literal options, etc.). This means that when a developer adds a new option to the Pydantic model, it can automatically appear in the UI with minimal extra code.
2.  **Execution:** Upon user submission, the UI will gather the values from the form, construct the `MainConfig` Pydantic object (triggering all the built-in validation), and then pass this single object to the `PipelineRunner` to initiate the pipeline in a background process or thread.
3.  **Monitoring & Visualization:** The UI will provide a display area for real-time feedback. This will likely be implemented via the `PipelineRunner` being modified to accept a callback function or a queue, to which it can post status updates (e.g., "Now running exploration for structure 5 of 20..."). After the run is complete, the UI will provide tools to load the output `.db` file and visualize the atomic structures, giving the user immediate visual confirmation of the results.

This strict separation ensures that the core `mlip_autopipec` library remains a powerful, independent tool that can be used via the CLI in HPC environments, while the UI serves as an alternative, user-friendly entry point.

## 4. Implementation Approach

The implementation of Cycle 2 will be tackled in a logical sequence, starting with the backend scientific upgrades and then building the user-facing interface on top of them.

1.  **Update Configuration and Dependencies:**
    *   The first step is to expand the Pydantic models in `common/config.py` to include all the new configuration parameters for the hybrid MD/MC engine and the FPS sampler.
    *   The `pyproject.toml` file will be updated to include the new dependencies required for this cycle: `mace-torch` (or another MLIP library), `scikit-learn` (for its powerful numerical tools, potentially useful for FPS), a library for SOAP calculations, and `streamlit` (or `flask`) for the web UI. A `uv sync` will install these new packages.

2.  **Implement FPS Sampler (`samplers/fps.py`):**
    *   Create the `FarthestPointSampler` class, inheriting from `BaseSampler`.
    *   Implement the logic for calculating SOAP descriptors for a list of `ase.Atoms` objects. This is a non-trivial step that may require an external dependency and should be designed to run in parallel.
    *   Implement the core FPS algorithm itself. This is an iterative process where, at each step, the point with the maximum minimum distance to the set of already-selected points is chosen.
    *   Update the `PipelineRunner` to check the `sampling.method` field in the configuration and instantiate either the `RandomSampler` or the new `FarthestPointSampler` accordingly.

3.  **Upgrade the Explorer (`explorers/md_engine.py`):** This is the most complex and scientifically critical part of the cycle.
    *   Refactor the `MDExplorer` to fully incorporate the "late-binding" calculator pattern. This is essential for preventing issues with pickling large PyTorch models when using `ProcessPoolExecutor`.
    *   Add the logic to load a specified MACE model and combine it with a ZBL potential using a custom `MixedCalculator` class that delegates calls appropriately.
    *   Implement the `detect_vacuum` utility function and use it within the explorer to automatically select the correct NPT (for bulk systems) or NVT (for surfaces) statistical ensemble for the MD simulation.
    *   In the main simulation loop, interleave the MD steps with calls to the new `_perform_mc_swap` method, which will contain the logic for the Monte Carlo atomic swap moves, including the Metropolis acceptance criterion.

4.  **Develop the Web UI (`main_gui.py`, `ui/app.py`):**
    *   The Streamlit framework will be used for its simplicity and rapid development cycle.
    *   The main application layout will be created in `ui/app.py`. A persistent sidebar (`st.sidebar`) will be used to house the entire configuration form.
    *   The UI will be built by creating Streamlit widgets (e.g., `st.text_input`, `st.slider`) that directly correspond to the fields in the Pydantic configuration models.
    *   The "Run Pipeline" button's logic will be implemented. To prevent the UI from freezing during a run, the `PipelineRunner` will need to be executed in a separate thread. A `Queue` object will be used to pass progress messages from the pipeline thread back to the UI thread for display.
    *   A results visualization section will be created. It will use `st.file_uploader` to allow the user to select the output `final_dataset.db` and then use a combination of ASE's plotting tools and Streamlit's `st.pyplot` to render images of the atomic structures.

5.  **Testing (`tests/`):**
    *   New unit tests will be written in `test_explorer.py` for the advanced MD engine. The MLIP model itself will be mocked. We will test the ensemble selection logic and the MC swap acceptance/rejection criteria with deterministic inputs.
    *   Unit tests will be written in `test_sampler.py` for the `FarthestPointSampler`. A small, known set of 2D input vectors will be used to assert that the implementation selects the correct points in the correct order.
    *   E2E tests for the UI will be created in `test_ui.py` using a tool like Playwright. These automated browser tests will launch the web app, simulate a user filling in the configuration values, click the "Run" button, and then verify that the UI updates to show progress and that the expected output file is eventually created.

## 5. Test Strategy

The testing strategy for Cycle 2 must be comprehensive, covering the new scientific logic, the user interface, and the integration between them.

**Unit Testing Approach:**

*   **Explorer:** The complex logic of the hybrid MD/MC engine demands thorough unit testing. We will create specific test cases in `test_explorer.py`. For example, `test_select_ensemble` will be given a mock `Atoms` object representing a bulk crystal and another representing a surface with a vacuum slab, and we will assert that it returns the correct ASE dynamics object (`NPT` vs `NVT`) for each. The `_perform_mc_swap` method will be tested to ensure it correctly applies the Metropolis criterion: a swap that lowers energy should always be accepted, while a swap that raises energy should be accepted with the correct probability.
*   **Sampler:** The `FarthestPointSampler` will be tested in `test_sampler.py` to validate the correctness of the core algorithm. We will not test it on high-dimensional SOAP vectors in the unit test. Instead, we will test it against a simple, known 2D point set where the correct FPS selection sequence is trivial to calculate by hand. We will assert that our implementation's selection sequence matches the expected sequence, thus validating that the algorithm is correct before it is applied to more complex scientific data.

**Integration Testing Approach:**

*   The end-to-end CLI test from Cycle 1 will be upgraded and expanded. A new, more complex configuration file will be created that enables all the advanced features (a real MLIP model, MC swaps, FPS sampling). The test will run the full pipeline and perform more sophisticated scientific checks on the output database. For example, it will calculate the average structural diversity of the final dataset (e.g., by computing the standard deviation of the SOAP descriptor distribution) and assert that this diversity is measurably higher when using the "fps" sampler compared to the "random" sampler. This test moves beyond "did it run?" to "did it produce a scientifically better result?".

**End-to-End (E2E) UI Testing:**

*   A new suite of E2E tests using the Playwright framework will be developed in `test_ui.py` to ensure the GUI works as expected from a user's perspective.
    *   **Test 1 (Happy Path):** This test will automate launching the app, filling in all the configuration fields in the sidebar with valid data, clicking the "Run Pipeline" button, waiting for the run to complete by polling for the "Success" message, and finally verifying that the expected output database file was created in the correct location.
    *   **Test 2 (Input Validation):** This test will automate entering invalid data into the UI (e.g., a non-numeric string into the "Temperature" field). It will then assert that a user-friendly error message appears next to the input field and that the "Run Pipeline" button is disabled, preventing the user from starting a run with an invalid configuration.
    *   **Test 3 (Visualization):** This test will first run the pipeline to generate a database. Then, it will simulate the user uploading this database via the file uploader widget and assert that a plot or image representing an atomic structure is correctly rendered on the page.
