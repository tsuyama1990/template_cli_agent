# User Acceptance Testing (UAT): CYCLE04 - Surrogate-Based Exploration (MACE Integration)

## 1. Test Scenarios

This UAT is focused on validating the pipeline's new ability to perform rapid, large-scale explorations using the MACE surrogate potential. The user, a computational scientist, wants to ensure that the system can take a few static initial structures and generate a rich, diverse set of configurations by simulating them at various temperatures. The emphasis is on the successful execution of these simulations and the generation of trajectory data.

| Scenario ID | Test Scenario                                                     | Priority |
| :---------- | :---------------------------------------------------------------- | :------- |
| UAT-C04-01  | Successful execution of MACE-driven MD for an initial structure   | High     |
| UAT-C04-02  | Verification of exploration at multiple temperatures              | High     |
| UAT-C04-03  | Correct handling of GPU/CPU device selection                      | Medium   |

---

### **Scenario UAT-C04-01: Successful execution of MACE-driven MD for an initial structure**

**(Min 300 words)**
This scenario is the core test of this cycle. The user wants to confirm that the system can correctly load the pre-trained MACE model and use it to run a stable Molecular Dynamics (MD) simulation. This validates the integration of the `mace-torch` library and ASE's dynamics engine.

The user will start with a database containing a single, well-defined structure, for example, a 32-atom SQS cell of FePt generated in the previous cycle. They will configure the `input.yaml` to run a short exploration, for instance, 1000 MD steps at a single temperature of 300 K. When the user runs the pipeline, they will monitor the system's logs. The first acceptance criterion is to see log messages confirming that the MACE-MP-0 model is being downloaded (on the first run) and successfully loaded as a calculator. They should see logs indicating the start of the MD simulation for the initial structure at 300 K.

The most important criterion is the successful completion of the MD run without any errors from the MACE library or ASE. The system should log a message like "MD simulation for structure 1 at 300 K completed." After the run, the user will inspect the output. The system, for this UAT, will be configured to save the resulting trajectory to a file (e.g., `trajectory_1_300K.traj`). The user will load this trajectory file and visualize it. They will verify that it contains 10 frames (assuming a save interval of 100 steps for a 1000-step run). By playing the trajectory, they should see the atoms vibrating around their lattice sites, consistent with a solid at 300 K. This visual confirmation proves that the surrogate potential is correctly calculating forces and propagating the dynamics.

---

### **Scenario UAT-C04-02: Verification of exploration at multiple temperatures**

**(Min 300 words)**
This scenario tests the orchestrator's ability to manage a more complex exploration task. The user wants to generate data that covers a range of thermal conditions, which is crucial for building a robust potential. They will verify that the system runs separate MD simulations for all the temperatures specified in the configuration.

The user will start with the same single FePt structure as in the previous scenario. This time, they will configure the `input.yaml` (or the expanded `exec_config_dump.yaml`) to include multiple temperature steps, for example: `simulation: {temperature_steps: [300, 650, 1000]}`. When they execute the pipeline, the acceptance criteria will be based on the logs and the output files. The user must see distinct log messages indicating the start and completion of three separate MD simulations: one at 300 K, one at 650 K, and one at 1000 K.

After the run, the user will check the output directory. They should find three distinct trajectory files: `trajectory_1_300K.traj`, `trajectory_1_650K.traj`, and `trajectory_1_1000K.traj`. By visualizing these trajectories, the user should be able to observe qualitative differences in the dynamics. The 300 K trajectory should show small atomic vibrations. The 650 K trajectory should show much larger amplitude vibrations. The 1000 K trajectory should show significant atomic motion, possibly including some diffusive events or local atomic rearrangements, indicating the system is approaching a disordered or molten state. Seeing this clear temperature-dependent behavior provides strong evidence that the exploration module is correctly sampling different thermodynamic states.

---

### **Scenario UAT-C04-03: Correct handling of GPU/CPU device selection**

**(Min 300 words)**
This scenario is a technical UAT to ensure the system makes efficient use of available hardware. MACE simulations are significantly faster on a GPU. The user will test the system's ability to automatically detect and use a GPU when available, and to fall back gracefully to the CPU otherwise.

This test will be run twice. First, on a machine with a compatible NVIDIA GPU and CUDA installed. The user will run a short exploration and monitor the logs. The key acceptance criterion is to find a log message from the `Explorer` module stating "CUDA device detected, using GPU for MACE calculations." Additionally, during the MD run, the user can run `nvidia-smi` in a separate terminal and should see the Python process utilizing the GPU.

Second, the user will run the exact same test on a machine without a GPU, or by setting the `CUDA_VISIBLE_DEVICES` environment variable to an empty string to hide the GPU. When the pipeline runs this time, the user must see a log message like "No CUDA device found, using CPU for MACE calculations." The MD simulation should still run to completion without any errors, albeit more slowly. The successful completion of both tests confirms that the hardware detection logic is robust and the system is optimized for performance on different platforms without requiring manual user configuration.

## 2. Behavior Definitions

### GIVEN/WHEN/THEN Definitions

**(Min 500 words)**

**Scenario: UAT-C04-01 & UAT-C04-02 - Successful execution of MACE-driven MD at multiple temperatures**

*   **GIVEN** a project database `asedb.db` containing at least one initial structure (e.g., a 32-atom FePt cell).
*   **AND** the `FullConfig` specifies an exploration using the `mace_mp` surrogate model.
*   **AND** the `FullConfig` specifies `md_steps_per_structure: 2000`.
*   **AND** the `FullConfig` specifies `temperature_steps: [300, 1000]`.
*   **WHEN** the user runs the pipeline and the `WorkflowOrchestrator` calls the `Explorer` module.
*   **THEN** the system should log that it is starting exploration for the initial structure.
*   **AND** the system should log that it is loading the `mace_mp` model.
*   **AND** the system should log that it is starting an MD simulation at 300 K for 2000 steps.
*   **AND** the MD simulation at 300 K should complete successfully.
*   **AND** the system should then log that it is starting a second MD simulation at 1000 K for 2000 steps.
*   **AND** the MD simulation at 1000 K should also complete successfully.
*   **AND** two trajectory files, one for each temperature, should be saved to the output directory.
*   **AND** the trajectory for 300 K should show solid-like atomic vibrations when visualized.
*   **AND** the trajectory for 1000 K should show much larger, liquid-like atomic motion when visualized.
*   **AND** the collection of all frames from both trajectories should be passed to the next stage of the workflow.

---

**Scenario: UAT-C04-03 - Correct handling of GPU/CPU device selection**

*   **GIVEN** a system with a compatible NVIDIA GPU available.
*   **AND** a project configured to run a MACE exploration.
*   **WHEN** the `Explorer` module initializes the MACE calculator.
*   **THEN** the system should log "CUDA device detected, using GPU for MACE calculations."
*   **AND** the subsequent MD simulation should show GPU utilization when monitored with `nvidia-smi`.
*
*   **GIVEN** a system with no GPU, or where the GPU is masked.
*   **AND** the same project is configured to run a MACE exploration.
*   **WHEN** the `Explorer` module initializes the MACE calculator.
*   **THEN** the system should log "No CUDA device found, using CPU for MACE calculations."
*   **AND** the MD simulation should run to completion successfully using only CPU resources.
