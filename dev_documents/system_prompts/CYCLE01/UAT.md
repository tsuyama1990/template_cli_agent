# UAT: CYCLE 01 - Core Engine and CLI

This document outlines the User Acceptance Testing (UAT) plan for the deliverables of Cycle 01. The primary goal is to verify that the command-line tool can perform a full, end-to-end data generation workflow for an alloy system, and that the results are scientifically plausible and meet the user's requirements.

This UAT will be performed using a Jupyter Notebook (`UAT_Cycle01.ipynb`) that will guide the user through running the CLI and analyzing the output.

## 1. Test Scenarios

| Scenario ID | Test Scenario | Priority |
|---|---|---|
| UAT-C01-001 | **End-to-End Pipeline Execution (Happy Path)** | High |
| UAT-C01-002 | **Verification of Output Database Content (Scientific Validation)** | High |
| UAT-C01-003 | **Verification of Physical Constraints** | High |
| UAT-C01-004 | **Handling of Invalid Scientific Configuration** | Medium |

---

### **Scenario UAT-C01-001: End-to-End Pipeline Execution (Happy Path)**

*   **Summary:** This scenario tests the main functionality. The user will execute the `mlip-autopipec` command with a valid configuration file for a simple binary alloy (e.g., Copper-Gold, CuAu). The test is successful if the command completes without errors and produces the final ASE database.
*   **Expected Outcome:** The CLI command runs to completion and exits with a status code of 0. An SQLite database file is created at the specified location.

---

### **Scenario UAT-C01-002: Verification of Output Database Content (Scientific Validation)**

*   **Summary:** This is the core scientific validation. The user, via the Jupyter Notebook, will open the generated ASE database and perform a series of checks on the data. The purpose is to confirm that the generated structures are physically meaningful and align with the input configuration.
*   **Expected Outcome:**
    1.  The number of structures in the database exactly matches the `num_samples` parameter from the configuration file.
    2.  **Scientific Check 1 (Composition):** The chemical composition of every structure in the database must match the `composition` defined in the configuration (e.g., 50% Cu, 50% Au).
    3.  **Scientific Check 2 (Energy Plausibility):** The potential energies of the structures are within a physically reasonable range. The notebook will plot a histogram of the energies to verify the distribution.
    4.  **Scientific Check 3 (Structural Integrity):** Visual inspection of a few sample structures from the database (plotted within the notebook) shows plausible atomic arrangements for a crystalline alloy.

---

### **Scenario UAT-C01-003: Verification of Physical Constraints**

*   **Summary:** This scenario verifies that the initial structure generation respects physical limits. While we cannot directly observe the rejected structures, we can infer the check is working by analyzing the output from a valid run.
*   **Expected Outcome:** The Jupyter Notebook will calculate the minimum interatomic distance for all atoms in all structures in the final database. It will assert that **no distance is smaller than the configured overlap threshold (e.g., 0.7 Å)**. This confirms that the atomic overlap check was enforced correctly.

---

### **Scenario UAT-C01-004: Handling of Invalid Scientific Configuration**

*   **Summary:** This scenario tests the tool's robustness when provided with a scientifically invalid configuration. The user will be instructed to run the CLI with a configuration where the composition fractions do not sum to 1.0.
*   **Expected Outcome:** The tool should not crash. Instead, it should gracefully exit and display a clear, user-friendly error message from the Pydantic validator, indicating that the compositions are invalid.

## 2. Behavior Definitions

**GIVEN** a user has a valid configuration file `config.yaml` for a CuAu alloy, with `num_samples: 10`
**WHEN** the user executes the command `mlip-autopipec run --config config.yaml`
**THEN** the program should complete successfully
**AND** a database file `output.db` must be created.

---

**GIVEN** the `output.db` file has been created
**WHEN** the user analyzes the database in a Jupyter Notebook
**THEN** the database must contain exactly 10 structures.
**AND** every structure in the database must have a chemical composition of 50% Cu and 50% Au.
**AND** the minimum interatomic distance found in any structure must be greater than 0.7 Å.

---

**GIVEN** a user has an invalid configuration file `invalid_config.yaml` where composition is `[0.5, 0.6]`
**WHEN** the user executes `mlip-autopipec run --config invalid_config.yaml`
**THEN** the program must exit with a non-zero exit code.
**AND** a clear error message must be printed, stating that the composition fractions must sum to 1.0.
