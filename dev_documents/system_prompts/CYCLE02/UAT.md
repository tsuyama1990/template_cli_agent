# UAT.md: Cycle 2 - Advanced Simulation, Sampling, and UI

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for Cycle 2 of the MLIP-AutoPipe project. These tests are designed to validate the advanced features built upon the core framework from Cycle 1. The focus is on ensuring that the sophisticated simulation and sampling techniques are functional, produce scientifically plausible results, and are accessible to the user. The Web UI proof of concept will also be tested for basic usability. A dedicated Jupyter Notebook (`UAT_Cycle2_Validation.ipynb`) will be the primary tool for executing these scenarios, allowing for the configuration of advanced features and in-depth analysis of the resulting data.

| Scenario ID | Scenario Name                                      | Priority |
| :---------- | :------------------------------------------------- | :------- |
| UAT-C2-001  | Successful Pipeline Run with MACE Potential        | High     |
| UAT-C2-002  | Verification of Farthest Point Sampling (FPS)      | High     |
| UAT-C2-003  | Validation of Hybrid MD/MC Exploration             | Medium   |
| UAT-C2-004  | Successful Generation of Ionic Structures          | Medium   |
| UAT-C2-005  | Basic Web UI Interaction                           | Low      |

### **UAT-C2-001: Successful Pipeline Run with MACE Potential** (Priority: High)

**Description:**
This scenario is a critical test of the system's ability to integrate with and utilize a real, state-of-the-art Machine Learning Interatomic Potential (MLIP). The user will configure the pipeline to use a pre-trained MACE model as the calculator for the exploration phase. A successful run validates several key technical features: the "late-binding" of the calculator in the parallel environment (avoiding memory and pickling issues), the correct parsing of the new configuration options, and the ability of the exploration engine to compute energies and forces using the complex model. The output database will be inspected to ensure that the results are physically reasonable, for example, by checking that the computed energies are within an expected range for the given material. This test confirms that the tool is ready for practical scientific use.

### **UAT-C2-002: Verification of Farthest Point Sampling (FPS)** (Priority: High)

**Description:**
This scenario tests the effectiveness of the intelligent sampling module. The goal is to verify that the Farthest Point Sampling (FPS) algorithm produces a more structurally diverse dataset compared to random sampling. The user will run two separate pipelines on the same initial structure, one configured to use "random" sampling and the other "fps". The Jupyter Notebook will then guide the user through a comparative analysis of the two resulting databases. This analysis will involve computing the structural similarity (e.g., using RDF or SOAP descriptor distances) between all pairs of structures in each dataset. The expectation is that the average similarity in the FPS-generated dataset will be significantly lower than in the randomly sampled one, providing a quantitative validation of the FPS implementation's benefit.

### **UAT-C2-003: Validation of Hybrid MD/MC Exploration** (Priority: Medium)

**Description:**
This test validates the hybrid MD/MC exploration feature, specifically the atom swap move. The user will configure a pipeline for a binary alloy (e.g., CuAu) and enable the `mc_moves` option. The key verification step is to analyze the resulting trajectory or the final sampled database to find evidence that atom swaps have occurred. In the Jupyter Notebook, the user will load the output structures and can, for instance, track the nearest neighbors of a specific atom to see if they change from Cu to Au over the course of the simulation. A successful test will demonstrate that the system's composition is evolving, which is crucial for exploring the full configurational space of alloys.

### **UAT-C2-004: Successful Generation of Ionic Structures** (Priority: Medium)

**Description:**
This scenario validates the new `IonicGenerator`. The user will configure the system to generate initial structures for a simple ionic compound, such as NaCl or MgO. The primary success criterion is to verify that the generated structures are charge-neutral. The UAT notebook will guide the user to inspect the output of the generation stage, assigning formal charges (oxidation states) to each atom and summing them to ensure the total charge of the simulation cell is zero. This confirms that the generator is correctly applying chemical rules, a prerequisite for running meaningful simulations of ionic materials.

### **UAT-C2-005: Basic Web UI Interaction** (Priority: Low)

**Description:**
This scenario provides a basic usability test for the proof-of-concept Web UI. The user will launch the web server and navigate to the provided URL in their browser. They will interact with the web form to input a simple pipeline configuration. Upon submission, the user should receive immediate feedback that the job has started. The test is successful if the user can launch a pipeline run from the UI and can see the `results.db` file appear in the project directory after a short time. This test is not intended to be a comprehensive UI/UX evaluation but rather a confirmation that the front-end is correctly communicating with the back-end and can trigger the core pipeline.

## 2. Behavior Definitions

---

### **Scenario: Successful Pipeline Run with MACE Potential**

*   **GIVEN** The user has a valid YAML configuration file specifying a MACE model for the `calculator`.
*   **AND** The specified MACE model is available and installed.
*   **WHEN** The user executes the command `mlip-autopipec run` with the MACE configuration.
*   **THEN** The pipeline should run to completion without any errors related to model loading or calculation.
*   **AND** The application should exit with a status code of 0.
*   **AND** A `results.db` file should be created.
*   **AND** Upon inspection, the structures within the database should have energy and force values computed.

---

### **Scenario: Verification of Farthest Point Sampling (FPS)**

*   **GIVEN** The user has a configuration file with `sampling.method` set to `"fps"`.
*   **WHEN** The user runs the pipeline to generate a database named `results_fps.db`.
*   **AND** The user runs a second pipeline with an identical configuration, but with `sampling.method` set to `"random"`, to generate `results_random.db`.
*   **THEN** Both pipelines should complete successfully.
*   **AND** When a structural similarity analysis is performed on both databases, the average similarity score for `results_fps.db` should be measurably lower than the score for `results_random.db`.

---

### **Scenario: Validation of Hybrid MD/MC Exploration**

*   **GIVEN** The user has a configuration file for a binary alloy with `exploration.mc_moves` set to `true`.
*   **WHEN** The user runs the pipeline.
*   **THEN** The pipeline completes successfully.
*   **AND** When the output trajectory or database is analyzed, there is clear evidence of atom species swapping positions during the simulation. For example, the list of atom species for a given set of coordinates changes between different frames.

---

### **Scenario: Successful Generation of Ionic Structures**

*   **GIVEN** The user has a configuration file specifying the generation of an ionic compound (e.g., NaCl).
*   **AND** The generator is configured as `ionic`.
*   **WHEN** The user runs the generation stage of the pipeline.
*   **THEN** The generated `ase.Atoms` objects should be charge-neutral (e.g., contain an equal number of Na and Cl atoms).

---

### **Scenario: Basic Web UI Interaction**

*   **GIVEN** The user has started the web server via the appropriate command (e.g., `mlip-autopipec web-ui`).
*   **WHEN** The user navigates to the web application in their browser.
*   **AND** The user fills out the configuration form and clicks "Submit".
*   **THEN** The web page should display a confirmation message indicating the job has started.
*   **AND** A pipeline process should be initiated on the server.
*   **AND** After a reasonable amount of time, a `results.db` file corresponding to the web configuration should appear in the filesystem.
