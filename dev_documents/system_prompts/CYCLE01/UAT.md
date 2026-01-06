# UAT.md - CYCLE 01

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 01 of the **MLIP-AutoPipe** project. The focus of this cycle is to verify the core functionality of the basic structure generation and storage pipeline. The tests are designed to be run from the perspective of an end-user, ensuring that the command-line tool is intuitive and produces the expected, physically valid outputs.

The primary UAT will be conducted using a Jupyter Notebook (`UAT_CYCLE01.ipynb`). This format allows a user to easily define a configuration, run the command-line tool, and then immediately load and inspect the results from the generated database, providing an interactive and transparent testing experience.

| Scenario ID | Test Scenario                                     | Priority |
| :---------- | :------------------------------------------------ | :------- |
| UAT-C01-001 | Generate a valid set of binary alloy structures   | High     |
| UAT-C01-002 | Tool rejects configuration with invalid composition | High     |
| UAT-C01-003 | Tool rejects configuration with atomic overlap    | Medium   |
| UAT-C01-004 | Tool handles non-existent config file gracefully  | High     |

---

### **Scenario ID: UAT-C01-001**

**Description (Min 300 words):**
This is the primary success-path test case for Cycle 01. Its purpose is to verify that a user can successfully generate a database of atomic structures for a simple binary alloy system using a valid configuration file. This test confirms that all the core components of the initial pipeline—configuration parsing, structure generation, physical validation, and database storage—are working together correctly. The user will create a simple YAML file specifying a Silicon-Germanium (SiGe) alloy, requesting a small number of structures (e.g., 10). They will then execute the `mlip-autopipec` command-line tool, pointing to this configuration and a desired output database path. The successful execution of this test will be determined by the tool exiting without error and the creation of a valid ASE database file. The user will then use the provided Jupyter Notebook to open this database, inspect its contents, and verify that the structures within it match the requested parameters. This includes checking that the database contains exactly 10 structures, and that each structure is composed of only Silicon and Germanium atoms in the correct proportions. Furthermore, the notebook will include a function to check the interatomic distances in a sample structure, confirming that the physical validation step (preventing atomic overlap) was successfully applied. This scenario provides end-to-end validation of the core functionality and gives the user confidence in the tool's basic capabilities.

### **Scenario ID: UAT-C01-002**

**Description (Min 300 words):**
This test case focuses on the robustness and validation capabilities of the configuration parsing system. It is designed to ensure that the tool provides clear and immediate feedback to the user if they provide a physically incorrect or logically inconsistent configuration. For this scenario, the user will create a YAML configuration file where the `composition` values do not sum to 1.0 (e.g., `{'Si': 0.6, 'Ge': 0.6}`). This represents a common user error. The expectation is that the Pydantic-based validation logic will detect this inconsistency before any computationally expensive structure generation begins. When the user runs the `mlip-autopipec` CLI with this invalid configuration, the tool should immediately exit with a non-zero status code, indicating an error. Crucially, it must print a user-friendly error message to the console that clearly explains what is wrong (e.g., "Error: Composition ratios must sum to 1.0"). This test is vital for user experience, as it prevents silent failures and helps the user correct their input without confusion. The Jupyter Notebook for this UAT will include a cell that demonstrates running the CLI with this bad configuration and captures the expected error output, confirming that the validation is working as designed. This ensures the tool is robust against common configuration mistakes.

## 2. Behavior Definitions

### **GIVEN/WHEN/THEN Definitions**

**Scenario: UAT-C01-001 - Generate a valid set of binary alloy structures**
```gherkin
GIVEN a valid configuration file "config_sige.yaml" exists with the following content:
  """
  system:
    elements: ['Si', 'Ge']
    composition: {'Si': 0.5, 'Ge': 0.5}
    num_structures: 10
  generator:
    name: 'alloy'
    min_atomic_distance: 1.8
  """
AND a clean workspace with no database file named "sige_structures.db"

WHEN the user executes the command:
  `mlip-autopipec --config-path config_sige.yaml --db-path sige_structures.db`

THEN the command should execute successfully and exit with a status code of 0
AND a file named "sige_structures.db" should be created in the workspace
AND when inspected, the "sige_structures.db" database should contain exactly 10 structures
AND each structure in the database should be composed of only Si and Ge atoms
AND the composition of each structure should be approximately 50% Si and 50% Ge
AND for any structure checked from the database, the minimum distance between any two atoms should be greater than 1.8 Angstroms.
```

---

**Scenario: UAT-C01-002 - Tool rejects configuration with invalid composition**
```gherkin
GIVEN a configuration file "invalid_composition.yaml" exists with the following content:
  """
  system:
    elements: ['Si', 'Ge']
    composition: {'Si': 0.6, 'Ge': 0.6} # Does not sum to 1.0
    num_structures: 10
  generator:
    name: 'alloy'
    min_atomic_distance: 1.8
  """

WHEN the user executes the command:
  `mlip-autopipec --config-path invalid_composition.yaml --db-path test.db`

THEN the command should fail and exit with a non-zero status code
AND an informative error message should be printed to the console, such as "Composition ratios must sum to 1.0"
AND no database file named "test.db" should be created.
```

---

**Scenario: UAT-C01-003 - Tool rejects configuration with atomic overlap**
```gherkin
GIVEN a valid configuration file "tight_config.yaml" exists with an impossibly large `min_atomic_distance`:
  """
  system:
    elements: ['Si', 'Ge']
    composition: {'Si': 0.5, 'Ge': 0.5}
    num_structures: 10
  generator:
    name: 'alloy'
    min_atomic_distance: 5.0 # A very large value that will always fail
  """

WHEN the user executes the command:
  `mlip-autopipec --config-path tight_config.yaml --db-path tight.db`

THEN the command may exit with a success or failure code, but must print a warning message indicating that 0 valid structures were generated
AND if a database file "tight.db" is created, it must contain 0 structures.
```

---

**Scenario: UAT-C01-004 - Tool handles non-existent config file gracefully**
```gherkin
GIVEN that no file named "non_existent_config.yaml" exists in the workspace

WHEN the user executes the command:
  `mlip-autopipec --config-path non_existent_config.yaml --db-path test.db`

THEN the command should fail and exit with a non-zero status code
AND an informative error message should be printed to the console, such as "Configuration file not found at path: non_existent_config.yaml"
AND no database file named "test.db" should be created.
```
