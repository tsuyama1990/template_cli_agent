# User Acceptance Testing (UAT): CYCLE 02 - Advanced Features

## 1. Test Scenarios

This document outlines the User Acceptance Tests for Cycle 02 of the MLIP-AutoPipe project. These tests are designed to validate the new, advanced features introduced in this cycle, ensuring they not only function correctly but also provide significant value to the end-user. The focus is on the hybrid exploration engine, the intelligent sampling stage, and the new web-based user interface.

As in the previous cycle, the primary method for conducting these UATs will be through a dedicated Jupyter Notebook (`UAT_Cycle02.ipynb`). This interactive format is ideal for testing scientific software, as it allows the user to execute the pipeline, inspect the results programmatically, and visualize the outputs, providing a comprehensive and intuitive validation experience.

| Scenario ID | Scenario Name                       | Priority |
| :---------- | :---------------------------------- | :------- |
| UAT-C02-001 | Validate Farthest Point Sampling (FPS) | High     |
| UAT-C02-002 | Verify Hybrid MD/MC Exploration    | High     |
| UAT-C02-003 | Interact with the Web UI            | Medium   |

---

### **Scenario UAT-C02-001: Validate Farthest Point Sampling (FPS)**

**Description (Min 300 words):**
The introduction of the sampling stage, particularly Farthest Point Sampling (FPS), is one of the most scientifically significant features of this cycle. The purpose of FPS is to select a subset of atomic structures that are maximally diverse, which is crucial for training robust and accurate MLIPs. This UAT scenario is designed to provide the user with tangible proof that FPS is working as intended and is superior to simple random sampling. The test will involve running the full pipeline on a system known to have distinct structural motifs (e.g., a simulation that includes both ordered and disordered phases). The pipeline will be configured to generate a large number of frames during exploration but to sample only a small subset for the final dataset using FPS. The core of this UAT will be the analysis performed in the Jupyter Notebook. After the pipeline run is complete, the notebook will load the final, FPS-sampled dataset from the database. It will then calculate a measure of diversity for this dataset. For comparison, it will also load the full trajectory and perform random sampling to create a dataset of the same size. By comparing the diversity scores of the FPS-sampled set and the randomly-sampled set, the user will be able to quantitatively verify that FPS is producing a more diverse dataset. The notebook will also include visualizations of the structures, allowing the user to qualitatively see the wider variety of configurations selected by FPS. This scenario is critical for building user trust in the "intelligence" of the pipeline.

---

### **Scenario UAT-C02-002: Verify Hybrid MD/MC Exploration**

**Description (Min 300 words):**
This scenario focuses on validating the advanced capabilities of the new hybrid MD/MC exploration engine. The key feature to be tested is the Monte Carlo atom swap move, which is essential for efficiently exploring the compositional and configurational space of materials like alloys. A simple MD simulation can get stuck in a local energy minimum, but the addition of MC swaps allows the system to jump between different configurations, leading to a much broader exploration. This UAT will involve setting up a simulation for a binary alloy and enabling the `enable_atom_swap` option in the configuration. The user will run the pipeline and then analyze the output trajectory to see the effect of the swaps. The Jupyter Notebook for this test will load the full exploration trajectory from the database. It will then perform an analysis to track the local environment of specific atoms throughout the simulation. For example, it could calculate the radial distribution function (RDF) at different points in the trajectory. A successful test would show that the RDF changes significantly, indicating that the atomic environments are evolving and the system is not static, which is a clear sign that the atom swaps are effectively rearranging the structure. This provides the user with strong evidence that the hybrid engine is more powerful than a simple MD engine. The test gives the user confidence that the pipeline is capable of exploring complex energy landscapes and generating the rich, varied data needed for high-quality MLIPs.

---

### **Scenario UAT-C02-003: Interact with the Web UI**

**Description (Min 300 words):**
While the backend scientific features are the core of the application, usability is key to broad adoption. This UAT scenario is focused on the user experience of the new web-based interface. The goal is to verify that the UI is an intuitive, helpful, and error-free tool for creating configurations and visualizing results. The test is interactive and exploratory. The user will launch the UI by running the `mlip-autopipec ui` command. They will then interact with the various widgets in the configuration panel. They should be able to easily set all the parameters for a pipeline run (e.g., select elements, adjust temperature with a slider, choose the sampling method from a dropdown). As they make changes, they should see the YAML output in the text area update in real-time. This validates the configuration builder functionality. The second part of the test involves the visualizer. The user will take an `.xyz` file generated from a previous pipeline run and upload it using the UI's file uploader. The expected outcome is a smooth, interactive 3D rendering of the atomic structure. The user should be able to rotate, zoom, and inspect the structure without any lag or graphical glitches. This scenario is important because it tests the application from a non-expert user's perspective. A positive experience with the UI can significantly lower the barrier to entry for new users and make the tool more accessible to a wider scientific audience.

## 2. Behavior Definitions

This section provides the detailed, Gherkin-style behavior definitions for each of the test scenarios in Cycle 02.

---

### **Scenario: Validate Farthest Point Sampling (FPS)**

```gherkin
GIVEN a user has run a long exploration on a binary alloy, generating 1000 total frames
AND the user has configured the pipeline to use the "fps" sampling method
AND the configuration specifies selecting 50 samples for the final dataset:
  """
  sampling:
    method: "fps"
    num_samples: 50
  """

WHEN the user executes the pipeline and it completes successfully
AND the user loads the final dataset of 50 structures from the "final_dataset" group in the database
AND the user also creates a comparative dataset by randomly selecting 50 structures from the full 1000-frame trajectory

THEN the final dataset should contain exactly 50 atomic structures
AND the diversity score (e.g., the average pairwise distance of SOAP descriptors) of the FPS-sampled dataset should be measurably higher than the diversity score of the randomly-sampled dataset
```

---

### **Scenario: Verify Hybrid MD/MC Exploration**

```gherkin
GIVEN a user wants to explore the structure of a Ni-Al alloy
AND the user has created a valid configuration file with the hybrid MD/MC exploration enabled:
  """
  system:
    type: "alloy"
    elements: ["Ni", "Al"]
    # ...
  exploration:
    temperature_k: 800
    enable_atom_swap: true
    mc_interval: 5
    num_steps: 1000
  """

WHEN the user executes the pipeline and it completes successfully
AND the user loads the full exploration trajectory from the database

THEN the trajectory should show evidence of significant structural rearrangement over time
AND an analysis of the atom positions should confirm that atoms of different species have swapped positions during the simulation
AND the radial distribution function (RDF) calculated from the end of the trajectory should be different from the RDF calculated from the beginning
```

---

### **Scenario: Interact with the Web UI**

```gherkin
GIVEN the MLIP-AutoPipe package is correctly installed

WHEN the user executes the command: `mlip-autopipec ui`
AND a web browser window opens displaying the application

THEN the UI should display a sidebar with widgets for all configuration parameters
AND when the user changes a value in a widget (e.g., moves a temperature slider)
AND the text area displaying the YAML configuration should update instantly to reflect the change

WHEN the user clicks the "Upload Structure" button
AND selects a valid ".xyz" file from their local filesystem

THEN a 3D visualization of the atomic structure from the file should appear in the main panel
AND the user should be able to interactively rotate and zoom the 3D model with their mouse
```
