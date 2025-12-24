# Cycle 2 User Acceptance Testing (UAT)

This document outlines the User Acceptance Testing (UAT) for Cycle 2 of the MLIP-AutoPipe project. The focus of this cycle is on the user-facing input system and the automated generation of initial data, significantly improving the system's usability and autonomy.

## 1. Test Scenarios

These scenarios test the new capabilities of generating structures and configurations from a minimal user input. The "user" is now a materials scientist who may not be an expert in DFT or MLIPs.

| Scenario ID | Scenario Name                                | Priority |
|-------------|----------------------------------------------|----------|
| UAT-C2-001  | Generate Full Config from Minimal Input      | High     |
| UAT-C2-002  | Generate Initial Structures for an Alloy     | High     |
| UAT-C2-003  | Generate Initial Structures for a Molecule   | High     |
| UAT-C2-004  | Handle Invalid User Input Gracefully         | Medium   |

---

**Scenario ID: UAT-C2-001**
- **Description:** This test verifies that the `Config Expander` can correctly process a simple, valid `input.yaml` file and generate a complete, well-formed `exec_config_dump.yaml`. This is the core of the new two-tier configuration strategy. The user provides the minimal input, and the system should fill in all other required parameters using its internal heuristics.
- **Priority:** High
- **Preconditions:**
    - The `mlip-pipe` command-line tool is installed.
    - An `input.yaml` file exists and contains a valid, minimal system definition (e.g., `system: {elements: ["Si"], composition: "Si"}`).

---

**Scenario ID: UAT-C2-002**
- **Description:** This test verifies that the system can generate a set of initial atomic structures for an alloy system. The user provides a minimal `input.yaml` for a binary alloy (e.g., FePt). The system should correctly identify it as an alloy and run the SQS and strain generation procedures.
- **Priority:** High
- **Preconditions:**
    - The `mlip-pipe` command-line tool is installed.
    - An `input.yaml` file exists for a binary alloy system.

---

**Scenario ID: UAT-C2-003**
- **Description:** This test verifies that the system can generate a set of initial atomic structures for a molecular system. The user provides a minimal `input.yaml` for a molecule (e.g., H2O), possibly including an initial structure file. The system should identify it as a molecule and run the Normal Mode Sampling procedure.
- **Priority:** High
- **Preconditions:**
    - The `mlip-pipe` command-line tool is installed.
    - An `input.yaml` file exists for a molecule, and if necessary, a corresponding coordinate file (e.g., `h2o.xyz`).

---

**Scenario ID: UAT-C2-004**
- **Description:** This test ensures the system provides clear feedback when the user provides an invalid `input.yaml`. For example, specifying an element that does not exist. The system should not crash with an unhandled exception but should provide a user-friendly error message.
- **Priority:** Medium
- **Preconditions:**
    - The `mlip-pipe` command-line tool is installed.
    - An `input.yaml` file is created with an invalid entry (e.g., `elements: ["Fe", "Xx"]`).

## 2. Behavior Definitions

---

**Scenario: UAT-C2-001 - Generate Full Config from Minimal Input**

```gherkin
GIVEN a user has created an 'input.yaml' file with the following content:
  """
  system:
    elements: ["Ni", "Al"]
    composition: "NiAl"
  """
AND the user is in a terminal in the same directory

WHEN the user runs the command "uv run mlip-pipe input.yaml"

THEN the command should execute successfully
AND a new file named 'exec_config_dump.yaml' should be created in the current directory
AND the 'exec_config_dump.yaml' file should contain a 'system' section with 'structure_type: "alloy"'
AND the 'exec_config_dump.yaml' file should contain a 'dft_compute' section with populated parameters like 'ecutwfc' and 'pseudopotentials'
```

---

**Scenario: UAT-C2-002 - Generate Initial Structures for an Alloy**

```gherkin
GIVEN a user has an 'input.yaml' configured for an "FePt" alloy
AND the system has been run to the point where the 'StructureGenerator' is invoked

WHEN the 'StructureGenerator' module is executed

THEN it should produce a list of 'ase.Atoms' objects in memory
AND the list should contain more than 10 structures
AND each structure in the list should contain only Iron (Fe) and Platinum (Pt) atoms
AND the chemical composition of each structure should be approximately 50% Fe and 50% Pt
```

---

**Scenario: UAT-C2-003 - Generate Initial Structures for a Molecule**

```gherkin
GIVEN a user has an 'input.yaml' configured for a "H2O" molecule, referencing a starting 'h2o.xyz' file
AND the system has been run to the point where the 'StructureGenerator' is invoked

WHEN the 'StructureGenerator' module is executed

THEN it should produce a list of 'ase.Atoms' objects
AND the list should contain multiple structures
AND each structure should contain exactly two Hydrogen atoms and one Oxygen atom
AND the atomic positions in the generated structures should be different from the initial 'h2o.xyz' file, representing vibrational modes
```

---

**Scenario: UAT-C2-004 - Handle Invalid User Input Gracefully**

```gherkin
GIVEN a user has created an 'input.yaml' file with the following content:
  """
  system:
    elements: ["Si", "BogusElement"]
    composition: "Si"
  """
AND the user is in a terminal in the same directory

WHEN the user runs the command "uv run mlip-pipe input.yaml"

THEN the command should fail with a non-zero exit code
AND the command should print a clear error message to the console, such as "Error: Invalid element 'BogusElement' found in input.yaml."
AND no 'exec_config_dump.yaml' file should be created
```
