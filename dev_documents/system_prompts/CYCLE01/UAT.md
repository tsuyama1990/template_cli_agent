# CYCLE01 USER ACCEPTANCE TEST (UAT): Core Pipeline and CLI Foundation

This document outlines the User Acceptance Testing (UAT) plan for the initial release of the MLIP-AutoPipe tool. The goal of this UAT is to verify that the core command-line pipeline is functional, reliable, and capable of performing a complete data generation workflow for a simple alloy system. This will establish confidence in the foundational architecture before more complex features are added in subsequent cycles.

## 1. Test Scenarios

The following scenarios are designed to be executed by a user from the command line. They represent the primary use cases for the functionality delivered in this cycle. A successful test run will involve setting up a configuration file, running the main CLI command, and inspecting the output to ensure it is correct and physically plausible. For simplicity and to facilitate user verification, we will use a Jupyter Notebook (`UAT_CYCLE01.ipynb`) as the primary test environment. This will allow the user to define the configuration, execute the command-line tool, and then use Python scripts to inspect the results all in one place.

---

**Scenario ID:** UAT-C01-001
**Priority:** High
**Title:** End-to-End Pipeline Run for a Binary Alloy System

**Description (Min 300 words):**
This is the most critical test scenario for Cycle 1. Its purpose is to validate the entire workflow from start to finish for a common use case: generating a small dataset for a binary alloy. The user will configure the system to generate structures for Copper-Gold (CuAu), a well-understood alloy. The test will verify that the tool can correctly interpret the configuration, generate initial random alloy structures, run a short Molecular Dynamics (MD) simulation using a fast, built-in potential, sample frames from the resulting trajectory, and store the final structures in a new database file.

The user will start by creating a YAML configuration file specifying the elements (`Cu`, `Au`), the composition (50/50), and the number of initial structures to create. They will also configure the MD exploration phase to run for a small number of steps at a constant temperature. The sampling will be set to `random`. The user will then execute the main pipeline command from their terminal, pointing to this configuration file. The primary success criterion is the creation of an ASE database file (`results.db`) in the output directory. The user will then inspect this database to confirm that it contains the expected number of atomic structures. Further verification will involve checking that the saved structures have the correct chemical composition and that their calculated energies are physically reasonable (i.e., not infinite or `NaN`). This scenario confirms that all four pipeline stages are correctly integrated and that the data flows between them as designed.

---

**Scenario ID:** UAT-C01-002
**Priority:** Medium
**Title:** Verify Pipeline Restart Capability

**Description (Min 300 words):**
A key feature of the pipeline is its ability to be stopped and restarted without losing progress. This scenario tests this checkpointing capability. Real-world exploration runs can take hours or days, and interruptions are common, so this feature is essential for robustness. The user will simulate an interruption after the initial structure generation phase is complete.

First, the user will run the pipeline as in the previous scenario. However, immediately after the command prints a log message indicating that the 'Generation' stage is complete, the user will manually stop the process (e.g., using `Ctrl+C`). At this point, a directory containing the initial seed structures (`.xyz` files) should exist. The user will then re-run the exact same command. Instead of generating new structures, the tool should detect the existing files and immediately proceed to the 'Exploration' stage, using the structures that were created in the first run. The user will verify this by observing the log output, which should indicate that it is skipping the generation step. The test is successful if the pipeline completes and the final database is created, just as in the first scenario. This confirms that the state isolation between the pipeline stages is working correctly and that the tool can recover from interruptions, saving significant computation time and improving reliability.

## 2. Behavior Definitions

The following definitions describe the expected system behavior in a Gherkin-style format. These will be implemented and verified within the `UAT_CYCLE01.ipynb` notebook.

---

**Feature:** Core Data Generation Pipeline

**As a** Materials Scientist
**I want to** run a command-line tool with a simple configuration file
**So that** I can automatically generate a database of atomic structures for an alloy system.

---

**Scenario:** Successful end-to-end run for a CuAu alloy

**GIVEN** I have a configuration file named `config.yaml` specifying:
  - A generator for an `alloy` system with `50% Cu` and `50% Au`.
  - The generation of `2` initial structures.
  - An exploration phase using `MD` for `50` steps with the `EMT` potential.
  - A `random` sampler that selects `10` structures.
  - An output database file named `CuAu_dataset.db`.

**AND** I have a clean working directory with no existing database or structure files.

**WHEN** I execute the command `mlip-autopipec-run --config-name=config.yaml` in my terminal.

**THEN** the command should execute without errors and exit with a success code.
**AND** I should see log messages indicating the start and successful completion of the `Generation`, `Exploration`, `Sampling`, and `Storage` stages.
**AND** a new database file named `CuAu_dataset.db` must exist in the output directory.
**AND** when I connect to the `CuAu_dataset.db` database, it must contain exactly `10` atomic structures.
**AND** each structure in the database must contain an equal number of `Cu` and `Au` atoms.
**AND** each structure in the database must have a finite, non-null value for its energy.

---

**Scenario:** Resuming an interrupted pipeline run

**GIVEN** I have the same `config.yaml` file as in the previous scenario.

**AND** I have previously run the pipeline and interrupted it after the `Generation` stage was complete.

**AND** a directory named `outputs/.../initial_structures` containing `.xyz` files already exists.

**WHEN** I re-execute the exact same command `mlip-autopipec-run --config-name=config.yaml`.

**THEN** the command should execute without errors and exit with a success code.
**AND** I should see a log message indicating that the `Generation` stage is being **skipped**.
**AND** the `Exploration` stage should start immediately, using the pre-existing structures.
**AND** a new database file named `CuAu_dataset.db` must exist in the output directory.
**AND** the final database must contain the correct number of structures, as specified in the configuration.
