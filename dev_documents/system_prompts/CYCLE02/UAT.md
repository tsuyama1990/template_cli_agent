# UAT.md - CYCLE 02

## 1. Test Scenarios

This document provides the User Acceptance Tests for Cycle 02 of the **MLIP-AutoPipe** project. This cycle's UAT focuses on verifying the advanced, user-facing features: the full four-stage data generation pipeline (Generate -> Explore -> Sample -> Store). The tests are designed to confirm that the Molecular Dynamics (MD) exploration engine and the intelligent sampling modules are working correctly and are well-integrated into the workflow.

As in the previous cycle, the primary testing tool will be a Jupyter Notebook (`UAT_CYCLE02.ipynb`). This allows the user to define a more complex configuration, execute the full pipeline, and then perform a detailed analysis of the output, including visualizing the trajectory and the final sampled structures.

| Scenario ID | Test Scenario                                                   | Priority |
| :---------- | :-------------------------------------------------------------- | :------- |
| UAT-C02-001 | Successfully run the full pipeline with MD exploration and random sampling | High     |
| UAT-C02-002 | Successfully run the full pipeline using FPS sampling           | High     |
| UAT-C02-003 | Tool correctly uses a user-provided MLIP model for exploration  | High     |
| UAT-C02-004 | Verify that the auto-ensemble switching logic works for a slab system | Medium   |

---

### **Scenario ID: UAT-C02-001**

**Description (Min 300 words):**
This is the cornerstone UAT for Cycle 02, designed to validate the successful integration and operation of the complete four-stage pipeline. The user will construct a configuration file that specifies a simple system (e.g., bulk Silicon) and enables the exploration stage with a short MD simulation. To ensure the test runs quickly, a fast classical potential like EMT will be used instead of a full MLIP. The configuration will also specify the "random" sampling method to select a small number of final structures. The user will then execute the `mlip-autopipec` CLI with this configuration. A successful run will result in the command completing without errors and producing a final ASE database. The core of this UAT lies in the verification step within the Jupyter Notebook. The user will load the final database and confirm that it contains the correct number of sampled structures. They will also inspect the intermediate trajectory files created during the exploration stage. The notebook will include functions to load these trajectories and plot the potential energy over time. A successful test will show a fluctuating but stable energy profile, confirming that the MD simulation ran correctly. This end-to-end test provides high confidence that the data flow between the four stages is correct and that the exploration engine is fundamentally working.

### **Scenario ID: UAT-C02-002**

**Description (Min 300 words):**
This scenario tests the intelligent sampling capability of the pipeline, specifically the Farthest Point Sampling (FPS) method. The goal is to verify that the FPS implementation produces a more structurally diverse dataset compared to a simple random sample. The user will start with the same setup as in `UAT-C02-001`, running a short MD simulation to generate a trajectory. However, the configuration file will be modified to set the sampling `method` to "fps". After the pipeline finishes, the user will use the `UAT_CYCLE02.ipynb` notebook to analyze the output. The notebook will contain a dedicated section for this test, which will load the 5-10 structures selected by the FPS algorithm. It will then compute the SOAP descriptors for these structures and visualize their distribution, for example, using a dimensionality reduction technique like PCA or by plotting a similarity matrix. For comparison, the notebook will also load a random sample of structures from the same trajectory and generate the same plots. A successful UAT will show that the structures selected by FPS are more spread out in the descriptor space than the randomly selected ones, demonstrating that the FPS algorithm is working as intended to maximize diversity. This confirms that one of the key "intelligent" features of the tool is delivering its intended value.

## 2. Behavior Definitions

### **GIVEN/WHEN/THEN Definitions**

**Scenario: UAT-C02-001 - Successfully run the full pipeline with MD exploration and random sampling**
```gherkin
GIVEN a valid configuration file "config_full_pipeline.yaml" exists, specifying an MD exploration and random sampling:
  """
  system:
    elements: ['Si']
    composition: {'Si': 1.0}
    num_structures: 2 # Generate 2 initial seeds
  exploration:
    model_path: 'emt' # Use the fast classical EMT potential for testing
    temperature: 300
    n_steps: 500
    save_interval: 50
  sampling:
    method: 'random'
    n_select: 5
  """
AND a clean workspace.

WHEN the user executes the command:
  `mlip-autopipec --config-path config_full_pipeline.yaml --db-path full_pipeline.db`

THEN the command should execute successfully and exit with a status code of 0.
AND intermediate trajectory files should be created in a temporary directory.
AND a final database file "full_pipeline.db" should be created.
AND the "full_pipeline.db" should contain exactly 5 structures.
AND when the trajectory files are inspected, they should contain a total of (500 / 50) * 2 = 20 frames.
AND the potential energy of the structures in the database should be valid numbers (not NaN or infinity).
```

---

**Scenario: UAT-C02-002 - Successfully run the full pipeline using FPS sampling**
```gherkin
GIVEN a valid configuration file "config_fps.yaml" exists, specifying an MD exploration and FPS sampling:
  """
  system:
    elements: ['Si']
    composition: {'Si': 1.0}
    num_structures: 1
  exploration:
    model_path: 'emt'
    temperature: 500 # Higher temperature for more diversity
    n_steps: 1000
    save_interval: 20
  sampling:
    method: 'fps'
    n_select: 10
    fps_soap_config: # Basic config for SOAP
      species: ['Si']
      r_cut: 4.0
      n_max: 6
      l_max: 4
  """

WHEN the user executes the command:
  `mlip-autopipec --config-path config_fps.yaml --db-path fps_pipeline.db`

THEN the command should execute successfully and exit with a status code of 0.
AND a final database file "fps_pipeline.db" should be created.
AND the database should contain exactly 10 structures.
AND when the structures in the database are analyzed, their structural fingerprints (SOAP descriptors) should be demonstrably more diverse than a random sample of 10 structures taken from the same source trajectory.
```

---

**Scenario: UAT-C02-004 - Verify that the auto-ensemble switching logic works for a slab system**
```gherkin
GIVEN a configuration file "config_slab.yaml" that generates a slab structure (which will have a vacuum):
  """
  system:
    elements: ['Pt']
    composition: {'Pt': 1.0}
    num_structures: 1
    # Generator settings would need to support creating a slab
  exploration:
    model_path: 'emt'
    temperature: 300
    n_steps: 100
  ...
  """
AND the application's logging is configured to output the chosen ensemble type.

WHEN the user executes the command with the slab configuration.

THEN the command should run successfully.
AND the application's log output should contain a message indicating that the "NVT" ensemble was chosen for the simulation.
AND the z-dimension of the simulation cell should remain constant throughout the simulation trajectory, which can be verified by inspecting the trajectory file.
```
