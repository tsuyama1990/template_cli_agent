# User Acceptance Testing (UAT): CYCLE 01 - Core Engine

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 01 of the MLIP-AutoPipe project. The goal of these tests is to verify that the core functionalities delivered in this cycle meet the primary user requirements and provide a solid foundation for future development. The tests are designed from the user's perspective, focusing on realistic use cases and ensuring a positive and intuitive user experience.

The recommended method for conducting these UATs is through a Jupyter Notebook (`UAT_Cycle01.ipynb`). This approach allows for an interactive and educational experience, where the user can not only run the tests but also visualize the results (e.g., the atomic structures) and inspect the outputs directly. The notebook will serve as both a validation tool and a tutorial for new users, demonstrating how to configure and run the pipeline and how to analyze the generated data.

| Scenario ID | Scenario Name                       | Priority |
| :---------- | :---------------------------------- | :------- |
| UAT-C01-001 | Generate a Binary Alloy Dataset     | High     |
| UAT-C01-002 | Generate an Ionic Crystal Dataset   | High     |
| UAT-C01-003 | Handle Invalid Configuration        | Medium   |

---

### **Scenario UAT-C01-001: Generate a Binary Alloy Dataset**

**Description (Min 300 words):**
This is the most fundamental use case for the MLIP-AutoPipe tool. The primary goal of a user is to generate a set of atomic structures for a specific material system to train a Machine Learning Interatomic Potential (MLIP). This test scenario will validate the complete end-to-end workflow for a common binary alloy, such as Copper-Gold (CuAu). The user will start by creating a simple YAML configuration file that defines the system (elements, composition) and the exploration parameters (temperature, number of steps). They will then execute the pipeline using the command-line interface. The expected outcome is the creation of an ASE database file containing both the initial, randomly generated seed structures and the corresponding trajectories from the molecular dynamics exploration. This test is critical because it verifies the integration of all core components developed in Cycle 01: the configuration parser, the alloy generator, the exploration engine, and the database writer. The user experience should be seamless; a single command should orchestrate the entire complex workflow. The Jupyter Notebook for this UAT will guide the user through creating the `config.yml`, running the `mlip-autopipec run` command, and then loading the resulting `mlip_data.db` file. The notebook will include cells to query the database, retrieve the generated structures, and use a visualization library (like `nglview`) to display the 3D atomic configurations. This provides immediate, tangible feedback that the tool is working as expected and producing physically plausible structures. It validates that the software is not just running without errors, but that it is producing a scientifically valuable output.

---

### **Scenario UAT-C01-002: Generate an Ionic Crystal Dataset**

**Description (Min 300 words):**
This scenario tests the flexibility and extensibility of the generator framework by targeting a different class of materials: ionic crystals. Unlike alloys, ionic crystals have strong constraints on their structure due to the need for charge neutrality. This UAT will verify the system's ability to generate valid initial structures for a simple ionic compound like Sodium Chloride (NaCl). The process will be similar to the alloy test: the user will define an appropriate `config.yml`, specifying the system type as `ionic` and providing the elements and composition. The pipeline will then be executed. The key verification point for this test is the output of the `IonicGenerator`. The generated structures must not only be physically plausible (no overlapping atoms) but should also be chemically correct. The Jupyter Notebook for this test will have specific cells designed to analyze the generated ionic structures. For example, it could include a simple analysis to check the coordination number of the ions or to verify that the stoichiometry is maintained. This test is important for demonstrating that the factory pattern for generators is working correctly and that the pipeline can be easily extended to support different types of physical systems with unique constraints. It gives the user confidence that the tool is not a one-trick pony but a versatile framework for materials science research. A successful run of this scenario proves that the core architecture is sound and capable of handling the diverse requirements of the materials science domain.

---

### **Scenario UAT-C01-003: Handle Invalid Configuration**

**Description (Min 300 words):**
A robust tool must not only work correctly with valid input but must also fail gracefully and informatively when given invalid input. This UAT scenario focuses on the user experience of error handling. The test will involve creating a deliberately incorrect configuration file. Examples of invalid configurations include a composition that does not sum to 1.0, a negative simulation temperature, or specifying an unsupported element symbol. The user will then attempt to run the pipeline with this invalid configuration. The expected behavior is that the program should not crash with an obscure traceback. Instead, it should exit cleanly and print a clear, user-friendly error message explaining exactly what is wrong with the configuration file. This test directly validates the effectiveness of the Pydantic-based configuration system. Pydantic's validation exceptions should be caught by the CLI, and the error messages should be presented to the user in a helpful manner. The Jupyter Notebook for this UAT will demonstrate this by showing a sample invalid `config.yml`, using a `try...except` block to capture the expected exit from the CLI command, and then printing the captured error output. This test is crucial for building user trust. A tool that provides helpful feedback when things go wrong is much more likely to be adopted. It shows that the system is designed with the user in mind and is robust enough to handle common mistakes without causing frustration.

## 2. Behavior Definitions

This section provides the detailed, Gherkin-style behavior definitions for each of the test scenarios.

---

### **Scenario: Generate a Binary Alloy Dataset**

```gherkin
GIVEN a user wants to generate a training dataset for a Copper-Gold (CuAu) alloy
AND the user has created a valid YAML configuration file named "cuau_config.yml" with the following properties:
  """
  system:
    type: "alloy"
    elements: ["Cu", "Au"]
    composition: {"Cu": 0.5, "Au": 0.5}
    num_structures: 5
  exploration:
    temperature_k: 500
    num_steps: 100
  db_path: "cuau_alloy.db"
  """

WHEN the user executes the command: `mlip-autopipec run --config cuau_config.yml`

THEN the command should execute successfully and exit with a status code of 0
AND a new file named "cuau_alloy.db" should be created in the current directory
AND the "cuau_alloy.db" file should be a valid ASE database
AND the database should contain a group of 5 structures under the key "initial"
AND the database should contain 5 groups of explored trajectories, named "explored_traj_0" through "explored_traj_4"
AND each explored trajectory should contain 101 structures (initial frame + 100 steps)
```

---

### **Scenario: Generate an Ionic Crystal Dataset**

```gherkin
GIVEN a user wants to generate a training dataset for Sodium Chloride (NaCl)
AND the user has created a valid YAML configuration file named "nacl_config.yml" with the following properties:
  """
  system:
    type: "ionic"
    elements: ["Na", "Cl"]
    composition: {"Na": 0.5, "Cl": 0.5}
    num_structures: 3
  exploration:
    temperature_k: 400
    num_steps: 50
  db_path: "nacl_ionic.db"
  """

WHEN the user executes the command: `mlip-autopipec run --config nacl_config.yml`

THEN the command should execute successfully and exit with a status code of 0
AND a new file named "nacl_ionic.db" should be created in the current directory
AND the "nacl_ionic.db" file should be a valid ASE database
AND the database should contain a group of 3 structures under the key "initial"
AND the structures in the "initial" group should be charge-neutral
AND the database should contain 3 groups of explored trajectories
```

---

### **Scenario: Handle Invalid Configuration**

```gherkin
GIVEN a user has created an invalid YAML configuration file named "invalid_config.yml" where the composition fractions do not sum to 1.0
  """
  system:
    type: "alloy"
    elements: ["Fe", "Pt"]
    composition: {"Fe": 0.6, "Pt": 0.5} # Sum is 1.1
    num_structures: 1
  exploration:
    temperature_k: 300
    num_steps: 100
  db_path: "should_not_be_created.db"
  """

WHEN the user executes the command: `mlip-autopipec run --config invalid_config.yml`

THEN the command should fail and exit with a non-zero status code
AND the command should print a user-friendly error message to the console
AND the error message should clearly state that the composition fractions must sum to 1.0
AND NO database file named "should_not_be_created.db" should be created
```
