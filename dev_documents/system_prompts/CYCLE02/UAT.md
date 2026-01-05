# User Acceptance Testing (UAT): CYCLE02

This document outlines the User Acceptance Testing (UAT) plan for the second and final development cycle of the MLIP-AutoPipe project. The focus of this UAT is to verify the advanced, scientific features of the application from an end-user perspective. The tests are designed to ensure that the Molecular Dynamics exploration and intelligent sampling capabilities work as expected, produce scientifically valid results, and provide a tangible benefit to the user's workflow. Unlike the unit and integration tests, which focus on technical correctness, these UAT scenarios are designed to answer the question: "Does this feature help the user achieve their research goals in a reliable and intuitive way?" The scenarios are more holistic, often involving the use of external tools like Jupyter notebooks to demonstrate that the output of the pipeline is not just correct, but also immediately useful and interoperable with the standard scientific Python ecosystem. Passing this UAT will signify that MLIP-AutoPipe has successfully transitioned from a functional tool to a genuinely useful scientific application.

## 1. Test Scenarios

The following scenarios are designed to be executed by a user to confirm that the complex functionality delivered in Cycle 2 meets the project's requirements for generating diverse and high-quality datasets. They cover the core scientific features and the new user interface.

| Scenario ID | Description                                                          | Priority |
| :---------- | :------------------------------------------------------------------- | :------- |
| UAT-C2-001  | Successful End-to-End Run with MD Exploration and Trajectory Analysis | High     |
| UAT-C2-002  | Verify Structural Diversity using Farthest Point Sampling (FPS)      | High     |
| UAT-C2-003  | Verify Automatic Ensemble Switching for a Slab System                | Medium   |
| UAT-C2-004  | Interactive Configuration, Execution, and Verification via Web UI    | High     |

---

### **Scenario UAT-C2-001: Successful End-to-End Run with MD Exploration and Trajectory Analysis**

*   **Objective:** To verify that the user can successfully run the full pipeline, including the computationally intensive Molecular Dynamics (MD) exploration step, and to confirm through analysis that the resulting structures have been altered from their initial state in a physically meaningful way. This is the most critical test of the core scientific functionality of the application. It proves that the exploration engine is not just running without errors, but is actively exploring the potential energy surface.
*   **User Story:** As a computational materials scientist, I want to run a short MD simulation on my initial structures to generate a set of thermally-relaxed and diverse configurations. I need to be confident that the simulation actually ran and that my final structures are not the same as my initial ones, so that my training data better represents a real material at a finite temperature.
*   **Preconditions:**
    *   The MLIP-AutoPipe application is installed and functional.
    *   A fast, built-in potential like ASE's EMT (Effective Medium Theory) is available and configured for testing purposes, to ensure the test can run in a reasonable amount of time.
    *   A valid YAML configuration file named `config_md.yaml` has been created.
*   **Test Steps:**
    1.  Create a file named `config_md.yaml` with the following content. This configures a short but meaningful MD run.
        ```yaml
        project_name: uat_test_cu_md
        system:
          elements: ['Cu']
          lattice: 'fcc'
          num_structures: 1
        exploration:
          temperature: 500.0
          md_steps: 100
          mlip_model: 'emt' # Use a fast, built-in potential for testing
        sampling:
          method: 'random' # Use random sampling to focus on the MD part
          fraction: 1.0
        ```
    2.  Create a Jupyter Notebook `verify_md_run.ipynb` to orchestrate the test and perform the analysis.
    3.  In the notebook, add and execute a cell to define and save a known initial structure for later comparison. This establishes a clear baseline.
        ```python
        import os
        from ase.db import connect
        from ase.build import bulk
        import numpy as np

        # Create a known, perfect initial structure for comparison
        initial_atoms = bulk('Cu', 'fcc', a=3.6)
        initial_positions = initial_atoms.get_positions().copy()
        initial_energy = initial_atoms.get_potential_energy() # Using EMT
        print(f"Initial Potential Energy: {initial_energy:.4f} eV")
        ```
    4.  Add a cell to execute the MLIP-AutoPipe pipeline from the command line.
        ```python
        # Run the MLIP-AutoPipe pipeline
        os.system('mlip-autopipec run --config config_md.yaml')
        print("Pipeline execution complete.")
        ```
    5.  Add a cell to load the result from the generated output database and perform a detailed comparison.
        ```python
        # Load the result from the output database
        final_db_path = 'uat_test_cu_md.db'
        db_final = connect(final_db_path)
        final_atoms = db_final.get_atoms(id=1)
        final_positions = final_atoms.get_positions()
        final_energy = final_atoms.get_potential_energy()

        print(f"Final Potential Energy: {final_energy:.4f} eV")

        # Verify that the atomic positions have changed significantly
        position_difference = np.sum((initial_positions - final_positions)**2)
        print(f"Sum of squared differences in positions: {position_difference}")
        assert position_difference > 1e-6, "MD simulation did not change the atomic positions!"
        print("Success: The structure was demonstrably modified by the MD exploration.")
        ```
*   **Expected Result:**
    *   The pipeline runs to completion without any errors or crashes.
    *   The Jupyter Notebook executes all cells successfully.
    *   The potential energy of the final structure will be different from the initial structure, reflecting the system's evolution during the simulation.
    *   The final assertion passes, confirming that the sum of squared differences between the initial and final atomic positions is greater than a small tolerance. This provides incontrovertible proof that the MD simulation was executed and dynamically altered the structure as intended.

---

### **Scenario UAT-C2-022: Verify Structural Diversity using Farthest Point Sampling (FPS)**

*   **Objective:** To verify that the advanced Farthest Point Sampling (FPS) algorithm correctly selects a subset of structures that are more structurally diverse than a simple random sample. This test is crucial for demonstrating the "intelligence" of the sampling stage and its value in creating high-quality, non-redundant datasets.
*   **User Story:** As a researcher, I want to use FPS to intelligently select the most unique structures from my long and often repetitive simulation trajectory. This is important so that I can train my MLIP on the most informative data possible and avoid wasting computational resources on training on thousands of nearly identical structures.
*   **Preconditions:**
    *   The `FPSSampler` component is fully implemented and all its dependencies (e.g., a library for computing SOAP descriptors) are installed.
*   **Test Steps:**
    1.  Create a Jupyter Notebook `verify_fps.ipynb` to perform a controlled experiment.
    2.  In the notebook, create a small, artificial, and highly controlled dataset of `ase.Atoms` objects where some structures are deliberately made to be very similar, and one is made to be very different.
        ```python
        from ase import Atoms
        import numpy as np

        # Create 3 structures that are very similar to each other (just slightly rattled)
        atoms1 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        atoms2 = atoms1.copy(); atoms2.rattle(stdev=0.01, seed=42)
        atoms3 = atoms1.copy(); atoms3.rattle(stdev=0.01, seed=123)

        # Create 1 structure that is very different (a stretched bond)
        atoms4 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 2.0]])

        manual_structures = [atoms1, atoms2, atoms3, atoms4]
        ```
    3.  Import the `FPSSampler` from the application's code and run it on this manual list, instructing it to select the 2 "best" structures.
        ```python
        from mlip_autopipec.sampling.samplers import FPSSampler

        # Initialize the sampler to select 2 samples
        sampler = FPSSampler(num_samples=2)
        selected_structures = sampler.sample(manual_structures)
        print(f"Selected {len(selected_structures)} structures.")
        ```
    4.  Add a final cell to verify that the selected structures are, as expected, the two most different ones from the original set.
        ```python
        selected_distances = [s.get_distance(0, 1) for s in selected_structures]
        print("Interatomic distances in selected structures:", selected_distances)

        # The two most different structures are the original one and the highly stretched one.
        # FPS should pick these two and discard the other two rattled versions.
        expected_distances = [0.74, 2.0]

        # Rounding to handle minor floating point variations
        assert set(np.round(selected_distances, 2)) == set(np.round(expected_distances, 2)), \
            "FPS did not select the most diverse structures!"
        print("Success: FPS correctly selected the most structurally diverse configurations.")
        ```
*   **Expected Result:**
    *   The `FPSSampler` runs without error, successfully computing the necessary internal descriptors.
    *   The final assertion in the notebook passes. The sampler should intelligently select the original `H2` molecule (with a bond length of 0.74 Angstrom) and the highly stretched `H2` molecule (with a bond length of 2.0 Angstrom), as these are the two points that are geometrically "farthest" from each other in the dataset. It should correctly identify that `atoms2` and `atoms3` are redundant and discard them.

---

### **Scenario UAT-C2-004: Interactive Configuration, Execution, and Verification via Web UI**

*   **Objective:** To verify that the user can successfully launch the web UI, graphically configure a simple pipeline run, execute it, and see the results, including a 3D visualization of a generated structure. This test validates the usability and accessibility of the tool for users who may be less comfortable with the command line.
*   **User Story:** As a new user, I want a simple graphical interface to learn how the pipeline works and to run my first few test cases. This will help me understand the key parameters and see the results of my choices visually, so that I can build confidence before moving on to more complex, scripted workflows.
*   **Preconditions:**
    *   All dependencies for the web UI (e.g., Streamlit) are installed correctly.
*   **Test Steps:**
    1.  In a terminal, run the command provided in the documentation to launch the web UI application: `streamlit run src/mlip_autopipec/web_ui.py`.
    2.  Open the URL provided by Streamlit (usually `http://localhost:8501`) in a modern web browser.
    3.  Interact with the input widgets on the web page to configure a simple run (e.g., for Silicon, 'Si', with a small number of structures).
    4.  Click the "Run Pipeline" button.
    5.  Wait for the pipeline to finish, observing any progress indicators displayed on the page.
    6.  Observe the results section of the UI after the run is complete.
*   **Expected Result:**
    *   The web UI launches successfully and displays a clean, well-organized set of configuration options.
    *   After the user clicks the "Run Pipeline" button, the UI should provide clear feedback that the process has started, for example, by showing a spinner or a progress bar.
    *   Upon completion, the UI should display a clear success message (e.g., "Pipeline completed successfully!").
    *   A 3D visualization of one of the generated Silicon structures should be displayed on the page, which the user can interact with (rotate, zoom). This provides immediate visual confirmation of the output.

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behaviour of the advanced features in a formal, unambiguous manner.

**Feature: Molecular Dynamics Exploration Engine**

As a user of MLIP-AutoPipe, I want the system to be able to perform Molecular Dynamics simulations using a specified interatomic potential to generate realistic atomic configurations that have been subjected to thermal motion.

**Scenario: Running a pipeline with the MD exploration step enabled**

*   **GIVEN** I have a valid configuration file for a single crystal of Copper (Cu).
*   **AND** the configuration specifies running a 100-step MD simulation at a temperature of 500K using the 'emt' potential.
*   **WHEN** I execute the `run` command with this configuration file.
*   **THEN** the application should successfully initialize the MD engine and complete the simulation without crashing or raising errors.
*   **AND** the final atomic positions in the structures that are stored in the output database should be measurably and significantly different from the perfect lattice positions in the initial crystal structure, indicating that the system has evolved dynamically.

---

**Feature: Intelligent Data Sampling with FPS**

As a user, I want the system to use Farthest Point Sampling to automatically create a diverse and non-redundant dataset from a potentially large and repetitive simulation trajectory.

**Scenario: Running a pipeline with FPS sampling on a trajectory that contains many duplicate structures**

*   **GIVEN** the exploration stage of the pipeline generates a set of 100 structures.
*   **AND** within this set, 50 of the structures are nearly identical to each other (e.g., from a long period where the simulation was stable in an energy minimum).
*   **AND** the other 50 structures are all very different from each other and from the first group (e.g., from a high-temperature part of the simulation).
*   **AND** the user's configuration specifies using the `fps` sampling method to select a final dataset of 25 structures.
*   **WHEN** the sampling stage is executed.
*   **THEN** the final set of 25 structures selected by the `FPSSampler` should contain only one representative structure from the group of 50 nearly identical ones.
*   **AND** the other 24 structures in the final dataset should be selected from the set of 50 diverse structures.
