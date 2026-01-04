# CYCLE01 User Acceptance Test: Core CLI Pipeline

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) for the deliverables of Cycle 1. The primary goal of this cycle is to produce a functional Command-Line Interface (CLI) that can orchestrate the entire dataset generation pipeline for a simple binary alloy. The UAT is designed to be a hands-on, positive user experience, acting as both a validation tool and an introductory tutorial to the system's core functionality.

| Scenario ID | Priority | Description                                                                                                                               |
| :---------- | :------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
| UAT-C1-001  | High     | **Successful End-to-End Pipeline Run for an Alloy.** A user can successfully configure, execute, and inspect the output of a complete pipeline run for a standard Iron-Platinum (FePt) alloy system using the CLI. |

**UAT Approach:**
The UAT will be conducted via a single, comprehensive Jupyter Notebook (`UAT_CYCLE01.ipynb`). This format is chosen for several reasons:
*   **Interactivity:** It allows the user to execute commands and view the outputs directly within a single, narrative document.
*   **Clarity:** It combines explanations (in Markdown cells) with code execution (in code cells), making the process easy to follow.
*   **Verification:** It provides a clear, step-by-step process for the user to verify that the software behaves as expected. The notebook will guide the user in checking for the existence of output files and inspecting the contents of the final database.
*   **Educational Value:** This notebook will serve as the primary tutorial for new users to understand the basic workflow of the MLIP-AutoPipe tool. A successful run through this notebook will give them the confidence to start using the tool for their own research.

The notebook will be structured to guide the user through three main phases:
1.  **Configuration:** Setting up the necessary YAML configuration file for the FePt system.
2.  **Execution:** Running the main pipeline from the command line (emulated within the notebook).
3.  **Verification:** Inspecting the output files and programmatically checking the contents of the generated ASE database.

## 2. Behavior Definitions

This section provides a Gherkin-style definition of the expected behavior for the primary test scenario, UAT-C1-001.

---

### **Scenario: UAT-C1-001 - Successful End-to-End Pipeline Run for an Alloy**

**Feature:** Core CLI Data Generation Pipeline

**As a** Materials Science Researcher,
**I want to** use a single command-line tool to generate a diverse set of atomic structures for a binary alloy,
**So that** I can create a foundational training dataset for my Machine Learning Interatomic Potential.

---

**GIVEN** a clean project environment with the `mlip-autopipec` package installed.

**AND** I have created a Hydra configuration directory named `conf`.

**AND** inside `conf`, I have created a YAML file named `config.yaml` with the following content, defining an FePt alloy system and a simple simulation workflow:
```yaml
# conf/config.yaml
defaults:
  - _self_

# System Configuration: Defines the material to be generated
system:
  type: "Alloy"
  elements: ["Fe", "Pt"]
  composition: {"Fe": 0.5, "Pt": 0.5}
  crystal_structure: "fcc" # Face-Centered Cubic
  num_initial_structures: 2

# Exploration Configuration: Defines the simulation parameters
exploration:
  temperature_k: 300.0
  num_steps: 10 # A short run for testing purposes
  calculator: "EMT" # Use ASE's fast Effective Medium Theory potential

# Sampling Configuration: Defines how to select structures from the simulation
sampling:
  method: "Random"
  num_samples: 5

# Output Configuration
db_path: "results/fept_dataset.db"
```

**WHEN** I execute the following command in my terminal:
`mlip-autopipec run --config-path conf --config-name config`

**THEN** the program should execute without raising any errors and display a series of log messages indicating its progress through the four stages:
1.  "Starting Stage 1: Generation..."
2.  "Starting Stage 2: Exploration..."
3.  "Starting Stage 3: Sampling..."
4.  "Starting Stage 4: Storage..."
5.  "Pipeline completed successfully."

**AND** a directory named `results` should be created in the current working directory.

**AND** inside the `results` directory, a file named `fept_dataset.db` should exist.

**AND** I can write a Python script (or use a Jupyter Notebook cell) to connect to the database using `ase.db.connect("results/fept_dataset.db")`.

**AND** when I query the database for the total number of entries (`db.count()`), the result should be exactly 5, matching the `num_samples` parameter in the configuration file.

**AND** when I retrieve one of the entries from the database (`row = db.get(1)`), it should contain the expected atomic data, including keys for `energy` and `forces`.

---
