# CYCLE 01: UAT.md

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for Cycle 01 of the Autonomous Development Environment (AC-CDD) project. The focus of this cycle is on the core CLI and the project initialisation process. The following scenarios are designed to verify that the `init` command functions as expected from a user's perspective.

| Scenario ID | Description                                                                 | Priority |
|-------------|-----------------------------------------------------------------------------|----------|
| UAT-001     | Successful initialisation of a new project.                                 | High     |
| UAT-002     | Handling of existing configuration files with user confirmation to overwrite. | Medium   |
| UAT-003     | Handling of existing configuration files with user declining to overwrite.  | Medium   |

**UAT-001: Successful initialisation of a new project.**
This scenario tests the "happy path" where a user runs the `init` command in a new, unconfigured project. The test will verify that the user is prompted for the required API keys and that the `.env` and `ac_c-dd_config.py` files are created with the correct content. This is the most critical scenario as it represents the primary use case for the `init` command.

**UAT-002: Handling of existing configuration files with user confirmation to overwrite.**
This scenario is designed to test the system's behaviour when the `init` command is run in a project that has already been configured. The test will verify that the user is warned that existing configuration files will be overwritten and that the system proceeds with the initialisation process only after the user confirms their intention to do so. This is an important scenario for ensuring that users do not accidentally lose their existing configurations.

**UAT-003: Handling of existing configuration files with user declining to overwrite.**
This scenario tests the case where the user chooses not to overwrite the existing configuration files. The test will verify that the `init` command exits gracefully without making any changes to the filesystem. This ensures that the user's decision is respected and that the system behaves in a predictable and non-destructive manner.

## 2. Behavior Definitions

The following behavior definitions are written in Gherkin-style (GIVEN/WHEN/THEN) to provide a clear and unambiguous description of the expected system behavior for each UAT scenario.

### UAT-001: Successful initialisation of a new project

**GIVEN** a new project directory without any existing `.env` or `ac_cdd_config.py` files.
**WHEN** the user runs the `init` command.
**AND** the user provides valid API keys when prompted.
**THEN** the system should create a `.env` file in the project's root directory.
**AND** the `.env` file should contain the API keys provided by the user.
**AND** the system should create an `ac_cdd_config.py` file in the project's root directory.
**AND** the `ac_cdd_config.py` file should contain the default configuration settings.
**AND** the system should display a success message confirming that the initialisation is complete.

### UAT-002: Handling of existing configuration files with user confirmation to overwrite

**GIVEN** a project directory that already contains an `.env` file and an `ac_cdd_config.py` file.
**WHEN** the user runs the `init` command.
**THEN** the system should display a warning message indicating that the existing configuration files will be overwritten.
**AND** the system should prompt the user to confirm whether they want to proceed.
**WHEN** the user confirms that they want to proceed.
**AND** the user provides new, valid API keys when prompted.
**THEN** the system should overwrite the existing `.env` file with the new API keys.
**AND** the system should overwrite the existing `ac_cdd_config.py` file with the default configuration settings.
**AND** the system should display a success message confirming that the initialisation is complete.

### UAT-003: Handling of existing configuration files with user declining to overwrite

**GIVEN** a project directory that already contains an `.env` file and an `ac_cdd_config.py` file.
**WHEN** the user runs the `init` command.
**THEN** the system should display a warning message indicating that the existing configuration files will be overwritten.
**AND** the system should prompt the user to confirm whether they want to proceed.
**WHEN** the user declines to proceed.
**THEN** the system should not make any changes to the existing `.env` file.
**AND** the system should not make any changes to the existing `ac_cdd_config.py` file.
**AND** the system should display a message indicating that the operation was cancelled.
**AND** the system should exit gracefully.
