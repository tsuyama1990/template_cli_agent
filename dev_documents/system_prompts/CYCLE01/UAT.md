# UAT.md for CYCLE01: Core Data Pipeline & Manual Workflow

## 1. Test Scenarios

This document outlines the User Acceptance Tests (UAT) for CYCLE01 of the MLIP-AutoPipe project. The central goal of this development cycle is to deliver a functional, reliable, and manually-operated pipeline. As such, the UAT is sharply focused on verifying that a user can successfully perform the two fundamental tasks of this cycle: first, labeling a given atomic structure with DFT-calculated properties, and second, training a valid MLIP from the data accumulated in the database. These tests are designed to be performed from the perspective of an end-user interacting solely with the command-line interface. The outputs of these commands—specifically the created database and the final potential file—must be valid, usable, and conform to expected standards. To make this process as clear and accessible as possible, the primary test scenario is presented in the form of a detailed Jupyter Notebook. This notebook serves a dual purpose: it is both a rigorous, executable test script and a hands-on tutorial that guides a new user through the core functionality of the system, demonstrating the intended workflow and the expected outcomes at each step. This approach ensures that we are not only testing for functional correctness but also evaluating the usability and clarity of the user experience.

*   **Scenario ID**: UAT-C01-01
*   **Priority**: High
*   **Title**: Manual End-to-End Pipeline for a Simple Crystalline System (Silicon)
*   **Description**: This test scenario provides a comprehensive verification of the entire manual workflow using a simple, well-understood, and computationally inexpensive material: crystalline Silicon (Si). The user will begin with a single, standard structure file for a 2-atom conventional Si unit cell in the diamond structure. They will first use the `label` command provided by the CLI. This command should orchestrate the process of reading the structure, generating the appropriate Quantum Espresso input, executing the DFT calculation, parsing the results, and storing them in a new SQLite database file. Subsequently, the user will execute the `train` command, pointing to the same database. This command should then read the labeled data and create a basic Atomic Cluster Expansion (ACE) potential. The ultimate success of this test is determined by the successful, error-free creation of both the database and the final potential file, and critically, the ability to load and use that generated potential file with a standard, external tool like the Atomic Simulation Environment's (ASE) calculator interface. Although a potential trained on a single data point is physically meaningless for real simulations, its successful creation and subsequent usability are sufficient to prove the functional integrity of the entire pipeline. This test validates that all the core components developed in CYCLE01—the CLI, orchestrator, config parser, labeling engine, database wrapper, and training engine—are correctly integrated and are producing valid, interoperable data artifacts. It is the fundamental "smoke test" for the entire system.

*   **Scenario ID**: UAT-C01-02
*   **Priority**: Medium
*   **Title**: System Robustness and User Feedback for Invalid Configuration
*   **Description**: This test scenario is designed to ensure that the system is robust against common user errors and provides clear, helpful feedback when it encounters invalid input. A key feature of the architecture is the use of Pydantic for rigorous, schema-based configuration validation. This test verifies that this feature is working correctly and is exposed to the user in a helpful way. The user will intentionally create and attempt to use an invalid configuration YAML file. This file will contain one or more clear errors, such as a negative value for a physically positive-definite quantity (like the plane-wave cutoff energy, `ecutwfc`) or a misspelled key. The user will then attempt to run the `label` command with this invalid configuration. The test is considered successful if, and only if, the application immediately refuses to run any expensive computations, exits with a non-zero status code indicating an error, and prints a user-friendly error message to the console that clearly and precisely indicates which configuration parameter is invalid and why (e.g., "Validation error for DFTSettings: ecutwfc must be positive"). This test validates the robustness of the Pydantic-based configuration validation system and confirms that we are improving the overall user experience by preventing opaque, hard-to-debug errors from occurring deep within a long-running DFT calculation. It ensures the system "fails fast" and provides actionable feedback.

## 2. Behavior Definitions

The expected behavior of the system for each scenario is defined below using the Gherkin-style (GIVEN/WHEN/THEN) syntax, which provides a clear, unambiguous, and business-readable description of the system's intended functionality.

---

### **Scenario: UAT-C01-01 - Successful Manual Pipeline Run**

**GIVEN** a user is in a clean working directory on their local machine
**AND** a valid, executable Quantum Espresso `pw.x` binary is available in the system's PATH
**AND** the user has prepared a valid configuration file named `config.yaml`, which contains all necessary and physically reasonable settings for a Silicon DFT calculation (e.g., appropriate cutoffs, k-points, etc.)
**AND** the user has prepared a valid atomic structure file named `si_unit_cell.cif`, which contains the atomic positions for a standard 2-atom Silicon conventional unit cell.

**WHEN** the user executes the following command in their terminal:
`mlip-pipe label --config-file config.yaml --atoms-file si_unit_cell.cif --db-file si_data.db`

**THEN** the command should execute without printing any error messages to the standard error stream
**AND** the command should complete successfully, exiting with a status code of 0.
**AND** upon completion, a new SQLite database file named `si_data.db` should have been created in the working directory.
**AND** this database file should contain exactly one row (or entry).
**AND** querying this entry should reveal that it contains the correct atomic structure for the 2-atom Si cell.
**AND** this entry's data section should contain correctly parsed, non-null floating-point values for the calculated `energy`, `forces`, and `stress`.

**WHEN** the user, having successfully generated the database, subsequently executes the training command in their terminal:
`mlip-pipe train --config-file config.yaml --db-file si_data.db`

**THEN** the command should execute without printing any error messages to the standard error stream
**AND** the command should complete successfully, exiting with a status code of 0.
**AND** upon completion, a new potential file named `si_model.ace` (or a similarly named default) should have been created in the working directory.
**AND** this `si_model.ace` file must be a valid, well-formed ACE potential file that can be successfully loaded and used by external third-party tools, specifically the ASE `ase.calculators.ace.ACE` calculator class.

---

### **Tutorial Notebook: `UAT-C01-Manual-Pipeline-Tutorial.ipynb`**

To make the UAT process concrete, repeatable, and user-friendly, the following demonstrates the detailed content of a Jupyter Notebook that will be created and distributed to users. This notebook will be the primary vehicle for executing the test and demonstrating the core functionality of the CYCLE01 deliverable.

```python
# UAT-C01: A Step-by-Step Tutorial for the Manual MLIP-AutoPipe Pipeline

import os
import subprocess
from ase.io import read, cif
from ase.db import connect
from ase import build
from ase.calculators.ace import ACE

# --- 1. SETUP: PREPARING THE INPUTS ---
# In this step, we will programmatically create all the necessary input files that a user
# would normally prepare manually. This makes the test self-contained and reproducible.

print("--- Step 1: Setting up input files ---")

# a) Create the Silicon (Si) unit cell CIF file. We use ASE's built-in tools for this.
# This represents a typical, well-defined starting structure for a material.
atoms = build.bulk('Si', 'diamond', a=5.43)
cif.write_cif('si_unit_cell.cif', atoms)
print("Successfully created si_unit_cell.cif")

# b) Create the detailed configuration YAML file. In this cycle, the user is expected
# to provide this file with all necessary details.
config_yaml = """
system:
  elements: ["Si"]

dft_compute:
  # NOTE: The user must ensure this command points to their local Quantum Espresso installation.
  command: "pw.x -in QE_input.in > QE_output.out"
  ecutwfc: 40.0
  ecutrho: 320.0
  kpoints_density: 0.2
  smearing: "mv"
  degauss: 0.01
  magnetism: "none" # Silicon is not a magnetic material.

mlip_training:
  model_type: "ace"
  r_cut: 5.5
  delta_learning: false # We disable delta learning for this simple test case.
  base_potential: "none"
  loss_weights: {energy: 1.0, force: 100.0, stress: 10.0}
"""
with open("config.yaml", "w") as f:
    f.write(config_yaml)
print("Successfully created config.yaml")
print("Setup complete. All input files are ready.")

# --- 2. EXECUTE THE LABELING COMMAND ---
# This is the first major action a user takes. We invoke the pipeline's `label` command
# to perform the DFT calculation and populate our database.

print("\n--- Step 2: Running the 'label' command via the CLI ---")

label_command = [
    "mlip-pipe", "label",
    "--config-file", "config.yaml",
    "--atoms-file", "si_unit_cell.cif",
    "--db-file", "si_data.db"
]

# We assume the user has a working Quantum Espresso installation in their PATH.
# For automated CI/CD testing, this `subprocess.run` call would be mocked to avoid
# dependency on a complex external binary.
result = subprocess.run(label_command, capture_output=True, text=True, check=True)

print("Label command executed successfully.")

# --- 3. VERIFY THE DATABASE INTEGRITY ---
# After the labeling command finishes, we must verify that it produced the correct output.
# The primary output is the database file, which should now contain our labeled structure.

print("\n--- Step 3: Verifying the database content ---")

assert os.path.exists("si_data.db"), "Assertion Failed: Database file 'si_data.db' was not created!"

db = connect("si_data.db")
count = db.count()
assert count == 1, f"Assertion Failed: Expected 1 entry in the database, but found {count}."

row = db.get(id=1)
assert "energy" in row.data, "Assertion Failed: 'energy' key not found in database entry's data."
assert "forces" in row.data, "Assertion Failed: 'forces' key not found in database entry's data."
assert "stress" in row.data, "Assertion Failed: 'stress' key not found in database entry's data."
print("Database verified successfully. It contains one labeled structure with all required data.")


# --- 4. EXECUTE THE TRAINING COMMAND ---
# With data in our database, we can now proceed to the second major user action: training the MLIP.

print("\n--- Step 4: Running the 'train' command via the CLI ---")

train_command = [
    "mlip-pipe", "train",
    "--config-file", "config.yaml",
    "--db-file", "si_data.db",
    "--output-file", "si_model.ace" # We assume the CLI supports specifying an output file.
]

result = subprocess.run(train_command, capture_output=True, text=True, check=True)
print("Train command executed successfully.")


# --- 5. VERIFY THE FINAL POTENTIAL FILE ---
# The final and most important verification. We check that the potential file was created
# and, critically, that it is a valid potential that can be used by standard simulation tools.

print("\n--- Step 5: Verifying the generated ACE potential file ---")

assert os.path.exists("si_model.ace"), "Assertion Failed: ACE potential file 'si_model.ace' was not created!"

# The ultimate test: can an ASE calculator be instantiated with this potential file?
try:
    calc = ACE("si_model.ace")
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    print(f"Successfully loaded potential with ASE. Calculated energy for the test structure: {energy:.4f} eV")
except Exception as e:
    assert False, f"FATAL: Failed to load the generated ACE potential with ASE. The file may be corrupt or invalid. Error: {e}"

print("\nUAT-C01-01 PASSED: All steps completed successfully and all artifacts were validated!")

# --- 6. CLEANUP ---
# In a real run, a user would keep these files. For a clean test, we can remove them.
# os.remove("si_unit_cell.cif")
# os.remove("config.yaml")
# os.remove("si_data.db")
# os.remove("si_model.ace")
# print("\nCleanup complete.")

```

---

### **Scenario: UAT-C01-02 - Invalid Configuration**

**GIVEN** a user is in a clean working directory
**AND** the user has created a valid structure file `si_unit_cell.cif`
**AND** the user has created an intentionally *invalid* configuration file named `bad_config.yaml` where the `ecutwfc` key has been assigned a physically nonsensical negative number. The content is:
```yaml
# bad_config.yaml
system:
  elements: ["Si"]
dft_compute:
  ecutwfc: -40.0 # Invalid value
  ecutrho: 320.0
  kpoints_density: 0.2
mlip_training:
  # ... other settings
```

**WHEN** the user executes the command with the invalid configuration:
`mlip-pipe label --config-file bad_config.yaml --atoms-file si_unit_cell.cif --db-file test.db`

**THEN** the command should fail immediately, without attempting to run any external programs.
**AND** the command should exit with a non-zero status code, signaling an error to the shell.
**AND** the application must print a clear and informative error message to the standard error stream.
**AND** this error message must unambiguously state that the `dft_compute.ecutwfc` field failed validation because its value must be positive.
**AND** no new files, specifically no database file named `test.db`, should be created in the working directory.
