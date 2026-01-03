# CYCLE02/UAT.md

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for Cycle 2 of the MLIP-AutoPipe project. This cycle introduces the advanced Exploration and Sampling stages and a Web UI. The UATs are designed to verify that the new features are functional, produce scientifically meaningful results, and are accessible to the user. The primary testing tool will again be a Jupyter Notebook (`CYCLE02_UAT.ipynb`) for the backend features, which will guide the user through running a full pipeline and analyzing the results. A separate set of manual checks will be defined for the Web UI.

| Scenario ID | Priority | Test Scenario Description                                                                |
|-------------|----------|------------------------------------------------------------------------------------------|
| UAT-C2-001  | High     | **End-to-End MD Exploration**: User wants to run the full four-stage pipeline on a simple system (e.g., bulk Silicon). The user will start with a single seed structure, run a short MD simulation at high temperature to generate diverse configurations, sample a few structures using random sampling, and verify the final database. This tests the successful integration of all four pipeline stages. |
| UAT-C2-002  | High     | **FPS Diversity Sampling**: User wants to compare the results of Random Sampling vs. Farthest Point Sampling (FPS). The user will run a simulation that generates a trajectory with some obvious structural variations. They will then run the Sampling stage twice: once with `method: random` and once with `method: fps`. The user will visually inspect the two resulting sets of structures to confirm that the FPS set is more structurally diverse (e.g., contains both solid-like and more disordered frames). |
| UAT-C2-003  | Medium   | **Automatic Ensemble Switching**: User wants to verify that the system correctly identifies a surface slab and applies the NVT ensemble. The user will generate a seed structure for a material surface (e.g., a Pt(111) slab with a vacuum layer), run the exploration pipeline, and inspect the logs or output to confirm that the NVT ensemble was used for the simulation, preventing the vacuum layer from collapsing. |
| UAT-C2-004  | High     | **Web UI - Pipeline Execution**: User wants to use the Web UI to configure and launch a pipeline run. The user will navigate to the web interface, fill in the form fields to define a simple system (similar to UAT-C2-001), click the "Run" button, and monitor the output log in the UI to see the pipeline progress and complete successfully. This is a critical test of the UI's core functionality. |
| UAT-C2-005  | Medium   | **Web UI - Configuration Validation**: User enters an invalid parameter in the Web UI (e.g., a negative temperature). The UI should display an immediate validation error message next to the field, preventing the user from submitting the invalid form. This ensures the UI provides a good user experience and prevents backend errors. |

## 2. Behavior Definitions

The following Gherkin-style definitions describe the expected behavior for each UAT scenario.

---

**Scenario: UAT-C2-001 - Successful End-to-End Pipeline Run with MD**

*   **GIVEN** a user has a trained MACE model for Silicon (`Si.model`).
*   **AND** the user creates a `config_si_md.yaml` file specifying a single Si diamond seed structure, an exploration phase of 100 steps at 1500K, and a random sampling of 5 final structures.
*   **WHEN** the user executes the command `mlip-autopipec generate config_si_md.yaml` from the terminal.
*   **THEN** the command should complete successfully, showing log output for all four stages (Generation, Exploration, Sampling, Storage).
*   **AND** an output database file should be created.
*   **AND** the database should contain exactly 5 structures.
*   **AND** when the user inspects the energies of the structures in the database, they should have different potential energy values, reflecting the MD simulation at high temperature.

---

**Scenario: UAT-C2-002 - Verification of FPS Diversity**

*   **GIVEN** a trajectory from a simulation that includes both ordered and melted structures. (This can be generated and saved in a step in the UAT notebook).
*   **AND** the user has two configuration files, `config_random.yaml` and `config_fps.yaml`, that are identical except for the sampling `method`.
*   **WHEN** the user runs the pipeline with `config_random.yaml`.
*   **AND** the user runs the pipeline with `config_fps.yaml`.
*   **THEN** both runs should complete successfully, producing `random_samples.db` and `fps_samples.db`.
*   **AND** when the user visually inspects the structures in `random_samples.db`, the structures may look visually similar.
*   **AND** when the user visually inspects the structures in `fps_samples.db`, the set of structures should contain a wider variety of configurations, including some that look solid and some that are clearly disordered or liquid-like, confirming the diversity goal of FPS.

---

**Scenario: UAT-C2-003 - Correct Ensemble Selection for a Surface Slab**

*   **GIVEN** a user has a configuration file `config_pt_surface.yaml` that generates a Pt(111) surface slab with 20 Angstroms of vacuum.
*   **AND** the configuration specifies an MD exploration phase.
*   **WHEN** the user runs the pipeline with this configuration.
*   **THEN** the command should complete successfully.
*   **AND** the log output for the Exploration stage should contain a line explicitly stating "Vacuum detected, using NVT ensemble."
*   **AND** when the user visualizes the final structures from the output database, the vacuum layer should still be present and the slab should not have collapsed or expanded unnaturally into the vacuum.

---

**Scenario: UAT-C2-004 - Basic Pipeline Execution via Web UI**

*   **GIVEN** the user has launched the Web UI by running `mlip-autopipec web`.
*   **AND** the user opens the provided URL in their web browser.
*   **WHEN** the user fills the web form with parameters for a simple run (e.g., generating 5 random AuCu alloy structures, without the exploration stage for speed).
*   **AND** the user clicks the "Submit" or "Run Pipeline" button.
*   **THEN** a progress indicator or log output should appear on the page.
*   **AND** the log should show the pipeline starting and completing successfully.
*   **AND** an output database file should be created in the application's working directory.
*   **AND** the database should contain 5 AuCu structures.

---

**Scenario: UAT-C2-005 - Real-time Validation in Web UI**

*   **GIVEN** the user is viewing the Web UI in their browser.
*   **WHEN** the user enters "-100" into the "Temperature" input field.
*   **THEN** a red error message should appear immediately next to the field, stating "Temperature must be positive."
*   **AND** the "Submit" button should be disabled.
*   **WHEN** the user corrects the temperature to a valid value like "300".
*   **THEN** the error message should disappear.
*   **AND** the "Submit" button should become enabled.
