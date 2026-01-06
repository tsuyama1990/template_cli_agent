# CYCLE 02 UAT.md: User Acceptance Testing for Advanced Features

## 1. Test Scenarios

This document outlines the User Acceptance Testing (UAT) scenarios for the successful completion of Cycle 02. The primary focus of this UAT is to provide end-users with a high degree of confidence in the full, end-to-end functionality of the MLIP-AutoPipe. The tests are designed to be performed by a user to verify that the newly implemented advanced features—exploration, sampling, and storage—are not only functional and integrated but also produce scientifically meaningful and high-quality results. The central verification tool for this UAT will be an interactive Jupyter Notebook, `UAT_CYCLE02_Advanced_Workflow.ipynb`. This notebook will serve as a comprehensive, hands-on guide, allowing the user to execute a complete data generation workflow, inspect the final database, and perform analyses that provide tangible proof of the software's capabilities. This approach is designed to be both a rigorous testing procedure and an educational experience, demonstrating the value and power of the completed application.

| Scenario ID | Scenario Name                               | Priority |
| :---------- | :------------------------------------------ | :------- |
| UAT-02.01   | Successful Execution of the End-to-End Pipeline | High     |
| UAT-02.02   | Detailed Verification of Final Database Contents | High     |
| UAT-02.03   | Tangible Verification of Structural Diversity (FPS vs. Random) | Medium   |

### Scenario UAT-02.01: Successful Execution of the End-to-End Pipeline

This scenario represents the ultimate user workflow, demonstrating the full, integrated power of the tool from a simple configuration to a final, rich database. The user will configure and execute a complete pipeline run, encompassing all four stages from initial structure generation to final database storage. To make this feasible within an interactive session, the Jupyter Notebook will provide a pre-configured `config.yaml` that uses a computationally inexpensive, fast-running potential (like ASE's built-in EMT potential) in place of a full-scale MLIP. The notebook will guide the user to execute the main `mlip-autopipec run` command from a shell cell. The primary success criterion for this test is the successful completion of this command without any errors or crashes. The user should see log messages from the application indicating the start and successful completion of each of the four stages. This scenario provides the user with the gratifying experience of seeing their high-level configuration be seamlessly transformed into a complete, structured dataset with a single command, which serves as the ultimate confirmation that all the complex services developed in this cycle are correctly integrated and fully operational. It is the "it just works" test that builds fundamental trust in the application's reliability.

### Scenario UAT-02.02: Detailed Verification of Final Database Contents

The ultimate product of the MLIP-AutoPipe is the final, curated database. The integrity and correctness of this database are therefore of the utmost importance. This scenario focuses on a detailed, programmatic inspection of the generated database to ensure it contains the correct data in the correct, accessible format. Following the successful run in the previous scenario, the Jupyter Notebook will guide the user in connecting to the newly created SQLite database file using the standard ASE DB API. The user will then execute a series of Python cells that perform specific, targeted queries on the database. These queries will programmatically verify key properties of the dataset. For example, one cell will count the total number of rows in the database and assert that it exactly matches the `num_samples` parameter from the configuration file. Another cell will retrieve a few specific atomic structures from the database and print their associated metadata, allowing the user to confirm that essential information like potential energy and atomic forces has been correctly calculated and stored. This scenario gives the user high confidence that the data has been persisted correctly, is complete, and is immediately ready for the next stage of their research, which would typically be the training of an MLIP model.

### Scenario UAT-02.03: Tangible Verification of Structural Diversity (FPS vs. Random)

This scenario is designed to provide a compelling, data-driven demonstration of the scientific value added by the Farthest Point Sampling (FPS) algorithm. A key promise of the tool is not just to create data, but to create *diverse* and therefore *efficient* data. To test this, the Jupyter Notebook will instruct the user to run two separate, small-scale pipeline runs on the same initial system. The first run will be configured to use simple `Random` sampling, while the second will be configured to use the more advanced `FPS` method. This will result in two distinct databases: `random_set.db` and `fps_set.db`. The notebook will then provide a suite of visualization and analysis tools to directly compare the contents of these two datasets. For instance, it will guide the user to extract the potential energy of every structure from both databases and plot them as overlaid histograms. The test is considered successful if the visualizations clearly and unambiguously show that the dataset generated with FPS covers a wider range of the potential energy surface (i.e., has a broader energy distribution) compared to the dataset generated by random sampling. This provides a tangible, intuitive, and scientifically meaningful demonstration of the power of the intelligent sampling feature and its direct benefit for creating more effective and efficient training datasets for MLIPs.

## 2. Behavior Definitions

The expected behavior of the system for each advanced UAT scenario is defined below with precision and clarity using the Gherkin syntax (GIVEN/WHEN/THEN).

### Behavior for UAT-02.01: Successful Execution of the End-to-End Pipeline

**GIVEN** a user has created a valid and complete `config.yaml` file that specifies a full, four-stage pipeline run for a Silicon (Si) system.
**AND** the configuration specifies the use of the fast "EMT" potential for the exploration stage to ensure a timely execution.
**AND** the configuration details are: generation of 2 initial structures, an MD exploration of 50 steps, and the final sampling of 10 structures.
**AND** the configuration specifies the final output database path as `Si_dataset.db`.

**WHEN** the user executes the command `mlip-autopipec run --config-path config.yaml` from their terminal.

**THEN** the command should execute the full pipeline—Generation, Exploration, Sampling, and Storage—in sequence and complete successfully with an exit code of 0.
**AND** a new SQLite database file named `Si_dataset.db` must be created in the current working directory.
**AND** the console output during the run must clearly log the progress, printing messages that indicate the start and successful completion of each of the four pipeline stages.

### Behavior for UAT-02.02: Detailed Verification of Final Database Contents

**GIVEN** a database file named `Si_dataset.db` has been successfully generated by a full pipeline run.

**WHEN** a user connects to this database using the standard ASE DB API, for example, by executing the Python code: `from ase.db import connect; db = connect('Si_dataset.db')`.

**THEN** the connection to the database should be successful and return a valid connection object.
**AND** querying the total number of rows in the database (e.g., via `len(db)`) must return exactly 10, precisely matching the `num_samples` parameter from the original configuration.
**AND** iterating through the database rows (e.g., via `for row in db.select(): ...`) must yield row objects that contain valid, readable `ase.Atoms` data.
**AND** each of these row objects must contain the expected metadata keys, including `energy` and `forces`, with plausible, non-null numerical values.

### Behavior for UAT-02.03: Tangible Verification of Structural Diversity

**GIVEN** a user has executed the pipeline twice for the identical initial system, resulting in two distinct database files:
1.  `random_set.db`, which was generated using the configuration `sampling: {method: Random}`.
2.  `fps_set.db`, which was generated using the configuration `sampling: {method: FPS}`.

**WHEN** the user runs a provided analysis script (e.g., in a Jupyter Notebook) that extracts the potential energy of every structure from both databases and plots the distributions as histograms.

**THEN** the resulting plot must visually demonstrate that the histogram for the `fps_set.db` data has a broader distribution (i.e., a larger standard deviation) of energies compared to the histogram for the `random_set.db` data.
**AND** this visible difference provides tangible evidence that the FPS method has successfully sampled a more diverse set of configurations, including the more scientifically interesting high-energy structures, from the raw simulation trajectory.
