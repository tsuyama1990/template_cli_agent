# Cycle 1 User Acceptance Testing (UAT)

This document outlines the User Acceptance Testing (UAT) plan for the Core Engine of the MLIP-AutoPipe system. The focus of this cycle is to verify that the fundamental components for DFT calculation and MLIP training are functioning correctly.

## 1. Test Scenarios

The following scenarios are designed to test the core functionality developed in Cycle 1 from an end-user's perspective. For this cycle, the "user" is a developer or a technical user who can manually prepare the necessary inputs.

| Scenario ID | Scenario Name                                | Priority |
|-------------|----------------------------------------------|----------|
| UAT-C1-001  | Successful DFT Calculation for a Simple System | High     |
| UAT-C1-002  | Successful MLIP Training from Labeled Data   | High     |
| UAT-C1-003  | Handling of a Failed DFT Calculation         | Medium   |

---

**Scenario ID: UAT-C1-001**
- **Description:** This test verifies that the `LabelingEngine` can correctly take a simple atomic structure, run a Quantum Espresso calculation, and parse the results. The user provides a valid atomic structure and a standard configuration for the DFT calculation. The expected outcome is that the system successfully produces the energy, forces, and stress for that structure without any errors.
- **Priority:** High
- **Preconditions:**
    - A valid Quantum Espresso installation is available and its path is correctly configured.
    - An `ase.Atoms` object representing a simple system (e.g., a silicon crystal unit cell) has been created.
    - A configuration file with valid DFT parameters (pseudopotentials, cutoffs) is available.

---

**Scenario ID: UAT-C1-002**
- **Description:** This test verifies that the `TrainingEngine` can successfully train an MLIP from a set of pre-calculated, labeled structures. This is the primary goal of the cycle. The user provides a small dataset of structures, each with its corresponding DFT energy and forces. The system is expected to run the training process and output a valid, usable MLIP model file.
- **Priority:** High
- **Preconditions:**
    - A list of at least 10 `ase.Atoms` objects, each containing pre-calculated DFT energy and forces, is available.
    - The MLIP framework (e.g., ACE) is correctly installed.
    - A configuration file specifying the MLIP model type and the Delta Learning settings is available.

---

**Scenario ID: UAT-C1-003**
- **Description:** This test ensures that the system behaves gracefully if the external DFT calculation fails. For example, the calculation may fail to converge. The user provides an input that is known to be problematic (e.g., an invalid pseudopotential setting). The system should not crash; instead, it should detect the error from the DFT run and report it clearly to the user.
- **Priority:** Medium
- **Preconditions:**
    - A Quantum Espresso installation is available.
    - An `ase.Atoms` object and a configuration file with parameters that are known to cause a QE failure (e.g., an SCF convergence failure) are provided.

## 2. Behavior Definitions

The following behavior definitions are described in Gherkin style to provide clear, unambiguous steps for testing each scenario.

---

**Scenario: UAT-C1-001 - Successful DFT Calculation for a Simple System**

```gherkin
GIVEN a developer has a Python script for testing
AND the script defines an 'ase.Atoms' object for a 2-atom silicon unit cell
AND the script defines a configuration dictionary specifying the path to the 'pw.x' executable and valid SSSP parameters

WHEN the developer executes a method that calls the 'LabelingEngine' with the silicon 'Atoms' object and the configuration

THEN the method should return an 'ase.Atoms' object
AND the returned object's 'info' dictionary should contain a key 'energy' with a floating-point value
AND the returned object's 'arrays' dictionary should contain a key 'forces' with a numpy array of shape (2, 3)
```

---

**Scenario: UAT-C1-002 - Successful MLIP Training from Labeled Data**

```gherkin
GIVEN a developer has a Python script for testing
AND the script loads a list of 15 'ase.Atoms' objects from a file
AND each 'Atoms' object in the list has pre-populated 'energy' and 'forces' data
AND the script defines a configuration dictionary for the 'TrainingEngine', specifying the 'ACE' model and a 'ZBL' reference potential for Delta Learning

WHEN the developer executes a method that calls the 'TrainingEngine' with the list of labeled 'Atoms' objects and the configuration

THEN the method should complete without raising any exceptions
AND the method should return a valid trained model object or a path to a model file
AND the model file, if created, should exist on the filesystem
```

---

**Scenario: UAT-C1-003 - Handling of a Failed DFT Calculation**

```gherkin
GIVEN a developer has a Python script for testing
AND the script defines an 'ase.Atoms' object
AND the script defines a configuration dictionary with a very low 'ecutwfc' (cutoff energy) that is known to cause a convergence error in Quantum Espresso

WHEN the developer executes a method that calls the 'LabelingEngine' with the 'Atoms' object and the erroneous configuration

THEN the method should raise a specific 'DFTRuntimeError' exception
AND the exception message should contain informative text, such as "Quantum Espresso SCF calculation did not converge."
```
