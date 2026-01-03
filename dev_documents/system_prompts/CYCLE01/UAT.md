# CYCLE 01: USER ACCEPTANCE TESTING (UAT)

This document details the User Acceptance Testing (UAT) plan for the first development cycle of the MLIP-AutoPipe project. The fundamental goal of this UAT is to provide a clear, verifiable demonstration that the core, non-interactive pipeline functions as specified. It is designed to confirm that the system can successfully take a minimal, high-level user input for a simple material and autonomously execute the entire, complex workflow of structure generation, DFT-based labeling, and model training to produce a tangible, usable output. This process serves as the ultimate validation of the architectural foundation built during Cycle 01.

Crucially, this UAT is designed with a dual purpose. While its primary function is to test the system, it is also intended to serve as the user's first hands-on tutorial. It provides a guided, practical example of how to use the software and what to expect from its outputs. To achieve this in a user-friendly and interactive manner, the central artifact for this UAT will be a Jupyter Notebook, named `C01_UAT_Notebook.ipynb`. This notebook will contain a step-by-step walkthrough of the entire process, including the initial setup, the command execution, and, most importantly, a series of code cells that allow the user to introspect the results of each stage of the pipeline. This approach empowers the user to not just trust that the system worked, but to see for themselves the data that was generated and the model that was created, thus building confidence and understanding.

## 1. Test Scenarios

### Scenario ID: UAT-C01-001
**Priority:** High
**Title:** End-to-End MLIP Generation for a Simple Crystalline Solid (Silicon)

**Description:**
This foundational test scenario is designed to validate the entire linear workflow implemented in Cycle 01 using a canonical, well-understood test case: crystalline Silicon (Si). Silicon is selected as the ideal initial test material for several reasons. Firstly, it is a simple, single-element system, which simplifies the interpretation of results. Secondly, it is non-magnetic, which avoids the complexities of magnetic calculations that will be handled by the system but are not the primary focus of this initial test. Thirdly, its ground state is the simple and highly symmetric diamond crystal structure. The physical properties of Silicon, as calculated by DFT, are extensively documented in scientific literature, providing a clear and unambiguous benchmark against which the system's results can be compared. This allows us to rigorously verify the correctness of our `LabelingEngine` implementation.

The user's journey through this scenario will be entirely contained within the `C01_UAT_Notebook.ipynb`. The notebook is structured to be both a test script and a learning tool.
1.  **Environment Setup**: The first section of the notebook will ensure a clean and reproducible test environment. It will contain shell commands to create a new, empty working directory for the test run, ensuring that the results are not contaminated by previous runs. This step emphasizes the self-contained nature of the pipeline.
2.  **Minimal Configuration**: The notebook will then guide the user to create the `input.yaml` file. The content will be deliberately minimal, containing only the essential chemical information (`elements: ["Si"]`). This step is designed to directly test the `ConfigExpander` and to highlight one of the core value propositions of the system: simplicity for the user.
3.  **Pipeline Execution**: The user will execute a single code cell that runs the main pipeline command (e.g., `!uv run mlip-pipe run --input-path input.yaml`). The output of the command will be displayed in the notebook, providing the user with a real-time view of the pipeline's progress through its various stages (Structure Generation, Labeling, Training).
4.  **Results Verification and Introspection**: This is the most detailed section of the notebook. After the pipeline execution is complete, the user will be guided through a series of verification steps using supplied Python code snippets.
    *   **Configuration Inspection**: The user will load the auto-generated `exec_config_dump.yaml` file and the code will pretty-print key sections. The user will verify that the `ConfigExpander` has correctly applied its heuristics, for instance, by selecting the appropriate SSSP pseudopotential for Silicon and setting a physically reasonable plane-wave energy cutoff.
    *   **Database Querying**: The notebook will provide code to connect to the newly created SQLite database (`mlip_auto_pipe.db`). The user will execute queries to confirm that a set of initial structures were generated. The tests will include checking the total number of structures and visualizing one of the structures (the perfect diamond lattice) to confirm its correctness.
    *   **DFT Result Validation**: The user will then execute a query to retrieve the calculated DFT results for all structures. The notebook will include code to check that the forces on the atoms in the perfect, undisturbed diamond lattice structure are negligibly small (close to zero), which is a strong indicator that the DFT calculation was performed correctly. It will also plot a histogram of the energies of the generated structures.
    *   **Final Output Confirmation**: The final step is to confirm that the `TrainingEngine` has successfully produced the final output. The user will execute a command to check for the existence of the trained model file (e.g., `final_model.ace`) in the working directory.

The successful execution of every cell in this notebook without errors, and the confirmation that the outputs match the expected results described, will constitute a formal pass for this UAT scenario. It provides comprehensive proof that the core, automated workflow is fully functional.

## 2. Behavior Definitions

The following Gherkin-style behavior definitions describe the specific, observable outcomes that will be verified during the UAT. Each feature corresponds to a major stage in the Cycle 01 pipeline, and the scenarios provide the precise conditions and expected results that will be inspected within the UAT Jupyter Notebook.

---

**Feature: Heuristic Configuration Expansion**
As a non-expert user, I want the system to automatically determine the complex, low-level computational parameters from my simple, high-level input, so that I can use the tool without needing a deep background in quantum mechanics.

**Scenario:** A user provides a minimal configuration file for crystalline Silicon.
**GIVEN** a file named `input.yaml` exists in the current directory.
**AND** the content of `input.yaml` is:
```yaml
system:
  elements: ["Si"]
  composition: "Si"
```
**WHEN** the user executes the pipeline command: `uv run mlip-pipe run --input-path input.yaml`.
**THEN** a new file named `exec_config_dump.yaml` must be created in the same directory.
**AND** this new file must be a valid YAML file.
**AND** when parsed, the `exec_config_dump.yaml` file must contain a `dft_compute` section.
**AND** the `dft_compute` section must contain the key `pseudopotentials` with the value `"SSSP_1.3_PBE_precision"`, which is the correct, heuristically determined protocol.
**AND** the `dft_compute` section must contain the key `ecutwfc` with a floating-point value greater than 50.0, representing a physically sound energy cutoff for Silicon.

---

**Feature: Automated Initial Structure Generation**
As a user, I want the system to automatically generate a varied set of initial atomic structures for training, so that the resulting potential is exposed to a diverse range of configurations and is therefore more robust.

**Scenario:** The pipeline initiates and runs the structure generation phase for Silicon.
**GIVEN** the pipeline is executed with a valid, fully expanded configuration for Silicon.
**WHEN** the "Structure Generation" stage of the pipeline completes successfully.
**THEN** a new SQLite database file named `mlip_auto_pipe.db` must be created in the working directory.
**AND** when queried, the 'structures' table in the database must contain at least 10 rows.
**AND** at least one of the structures in the database must be identifiable as a perfect, 8-atom conventional cell of the Silicon diamond lattice.
**AND** other structures in the database must be variations of this lattice, for example, having different cell volumes (representing hydrostatic strain) or slight random perturbations of the atomic positions.

---

**Feature: Robust and Automated DFT Labeling**
As a user, I want the system to automatically perform and manage the complex DFT calculations for all the generated structures, so that I can obtain a complete and accurately labelled dataset without any manual intervention or error handling.

**Scenario:** The pipeline processes the "uncalculated" structures in the database.
**GIVEN** the database file `mlip_auto_pipe.db` contains a set of at least 10 structures, all marked with the status "uncalculated".
**WHEN** the "DFT Labeling" stage of the pipeline runs and completes successfully.
**THEN** when queried, all structures in the database must now have the status "calculated".
**AND** each of these structures in the database must now have a corresponding valid floating-point value for its `energy`.
**AND** each structure must have a corresponding `forces` array, and the dimensions of this array must be Nx3, where N is the number of atoms in that specific structure.
**AND** for the specific structure representing the perfect Silicon diamond lattice, the magnitude of the force vector on every atom must be less than 1.0e-4 eV/Å, confirming a correct relaxation to the ground state.

---

**Feature: Automated Model Training and Output**
As a user, I want the system to automatically use the newly generated dataset to train a machine learning model, so that the final output of the pipeline is a single, usable potential file.

**Scenario:** The pipeline initiates the training phase after all data has been successfully labeled.
**GIVEN** the `mlip_auto_pipe.db` database contains a complete set of structures, all of which have been successfully labeled with DFT results.
**WHEN** the "Model Training" stage of the pipeline is executed and completes successfully.
**THEN** a new file, `final_model.ace` (or a similarly named file depending on the training engine), must be created in the working directory.
**AND** this file must not be empty.
**AND** the main pipeline process must terminate gracefully with an exit code of 0, indicating success.
