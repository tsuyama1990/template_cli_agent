# SPEC.md: Cycle 02 - Initial Structure Generation & Configuration

## 1. Summary

Cycle 02 represents a pivotal evolution for the MLIP-AutoPipe project, transforming it from a developer-driven, manually-operated core engine into a genuinely user-friendly and automated system. The primary, overarching goal of this cycle is to completely eliminate the need for hardcoded, script-based inputs. Instead, we will empower the end-user to control the entire pipeline through a simple, high-level, human-readable configuration file. This ambitious goal will be realised through the development of two key, tightly interconnected, and architecturally significant components: the **Structure Generator (Module A)** and the **Config Expander**. The latter is the technical implementation of the project's core usability principle, the "Two-Tier Configuration" strategy.

The first major deliverable, the Structure Generator, is responsible for algorithmically creating the initial, diverse set of atomic configurations that will serve as the crucial seed for the entire learning process. This module is a critical step towards fulfilling the project's central philosophy of "removing the human expert from the loop," as it explicitly replaces the subjective, manual selection of starting structures with a physics-aware, deterministic, and automated procedure. The module will be designed with built-in intelligence to first classify the input material based on its fundamental chemical nature (e.g., identifying it as a metallic alloy, a discrete molecule, or an ionic crystal) and then dynamically dispatch to the most appropriate and physically sound generation technique. For instance, it will employ the Special Quasirandom Structures (SQS) method for alloys or Normal Mode Sampling (NMS) for molecules.

Complementing the structure generation is the development of the Config Expander, which is the heuristic engine that forms the heart of the system's usability and accessibility. This component is designed to ingest a minimalistic, user-provided `input.yaml` file—a file that contains only the most essential information, such as the chemical elements involved and the target composition. It then processes this high-level input and expands it into a comprehensive, verbose, and explicit `exec_config_dump.yaml` file. This generated file will contain every detailed parameter required for every subsequent stage of the pipeline, from fine-grained DFT settings like plane-wave cutoff energies and k-point mesh densities, to the specific hyperparameters for the MLIP training. By embedding a significant amount of expert knowledge and physics-based heuristics directly into this engine, the system will be capable of making intelligent, automated decisions, thereby freeing the user from the immense burden of complex, low-level configuration. Upon the successful completion of Cycle 02, the MLIP-AutoPipe will have matured into a truly usable system, capable of initiating a full, complex data generation workflow from a single, simple user command, marking a significant milestone in its journey towards full autonomy.

## 2. System Architecture

The architecture of the MLIP-AutoPipe system undergoes a significant and strategic enhancement in Cycle 02. Two major new components, the `Config Expander` and the `Structure Generator (Module A)`, are introduced at the very beginning of the workflow, effectively prepending the core engine that was developed and validated in Cycle 01. This architectural evolution means that the workflow is no longer initiated programmatically with a hardcoded ASE Atoms object. Instead, the entire process is now triggered by a user pointing the system's command-line interface to a single, user-defined `input.yaml` file. This change marks the transition from a purely backend process to a user-driven application.

The updated architectural workflow proceeds as follows:
1.  **Initiation**: The user launches the pipeline from their command-line terminal, providing the path to their minimal `input.yaml` file as the primary argument. This is the new, formal entry point for the entire application.
2.  **Configuration Expansion**: The `Orchestrator`, upon starting, immediately invokes the `Config Expander`. This dedicated module reads and validates the minimal `input.yaml`. It then applies a sophisticated set of pre-programmed, physics-based rules and heuristics to infer and generate a complete set of execution parameters. The result is the creation of a comprehensive `exec_config_dump.yaml` file. This full configuration is then loaded by the Orchestrator to govern all subsequent steps.
3.  **Structure Generation Request**: Armed with the validated and complete configuration, the Orchestrator proceeds to the next logical step, calling the `Structure Generator (Module A)`.
4.  **Bond Type Analysis**: The first action of the Structure Generator is to analyse the chemical composition provided in the configuration (e.g., "FePt", "H2O") to classify the material's primary bonding character. This classification (e.g., 'alloy', 'molecule') is a critical decision point that determines the subsequent generation strategy.
5.  **Structure Creation**: Based on the determined bond type, the module executes the most appropriate and physically sound generation algorithm (e.g., SQS for an alloy, NMS for a molecule). This results in a list of diverse ASE Atoms objects being created in memory.
6.  **Data Persistence**: The list of newly generated ASE Atoms objects is then passed to the `AseDB` wrapper, where they are saved to the database. Each structure is assigned an initial state, such as 'initial_generated', to mark its origin and position in the workflow.
7.  **Hand-off to Core Engine**: From this point forward, the workflow seamlessly transitions to the process defined and validated in Cycle 01. The Orchestrator will query the database for structures that are ready to be labelled, and the `Labelling Engine (Module C)` will begin processing the structures that were just generated by Module A.

This updated architecture demonstrates a clean, layered design. The new user-facing configuration and automated data seeding components are placed at the very beginning of the pipeline, acting as the data-provisioning stage. The `Config Expander` serves as a crucial "translation" layer, converting high-level human intent into a low-level, machine-operable execution plan. The `Structure Generator` then acts as the first data production step, populating the database with the essential raw material for the learning process. A key success of this design is that this significant new functionality is added without requiring any modification to the core Labelling and Training Engines from Cycle 01, powerfully demonstrating the extensibility and modularity of the initial architectural choices.

## 3. Design Architecture

The primary focus of Cycle 02 is the introduction of new classes and modules dedicated to configuration handling and structure generation. This will also necessitate updates to the `Orchestrator` and the nascent command-line interface (CLI) to accommodate the new, user-driven entry point. The design continues to emphasise a strong separation of concerns and leverages Pydantic for robust data validation.

**New and Updated Classes/APIs:**

1.  **`mlip_autopipec.config.models`** (New Pydantic Models)
    *   **`MinimalConfig(pydantic.BaseModel)`**: This model will define the strict but simple schema for the user-provided `input.yaml`. Its purpose is to be as simple as possible.
        *   `system: Dict`: A dictionary containing the most essential system definition, including `elements` (a list of chemical symbols, e.g., `["Fe", "Pt"]`) and `composition` (a string representation, e.g., `"FePt"`).
        *   `simulation: Dict`: An optional section for the user to specify high-level simulation goals, such as the target `temperature` range. This allows the heuristic engine to make more informed decisions.
    *   **`FullConfig(pydantic.BaseModel)`**: This model defines the comprehensive schema for the generated `exec_config_dump.yaml`. It will be a large, nested structure, using other Pydantic models for clarity and organisation.
        *   `system: SystemConfig`: A sub-model containing both user-provided and inferred system properties, such as `structure_type` and an estimated `melting_point_guess`.
        *   `dft_compute: DFTConfig`: A dedicated sub-model containing every parameter required by the Labelling Engine, such as the `command` for `pw.x`, paths to `pseudopotentials`, `ecutwfc`, `kpoints_density`, magnetism settings, etc.
        *   `mlip_training: MLIPConfig`: A sub-model containing all parameters for the Training Engine, including `model_type`, `r_cut`, `delta_learning` flag, and model-specific hyperparameters.

2.  **`mlip_autopipec.config.expander.ConfigExpander`** (New Class)
    *   **Purpose**: To house the entire heuristic engine and the logic for the minimal-to-full configuration expansion.
    *   **Public API**:
        *   `expand(minimal_config: MinimalConfig) -> FullConfig`: The main public method, which is designed to be a pure function for testability.
    *   **Key Internal Methods**:
        *   `_determine_structure_type(elements: List[str]) -> str`: Implements the core classification logic, likely based on analysing the electronegativity differences and known properties of the constituent elements to classify the material as 'alloy', 'ionic', 'covalent', or 'molecule'.
        *   `_get_dft_defaults(elements: List[str]) -> DFTConfig`: This method encapsulates a significant amount of expert knowledge. It will contain the logic to select the best default pseudopotentials (e.g., from a lookup table for a specific SSSP protocol), recommend appropriate plane-wave cutoff energies based on those potentials, and set other robust DFT parameters for convergence.
        *   `_get_mlip_defaults(...) -> MLIPConfig`: Provides a set of sensible default hyperparameters for the ACE training process, which can serve as a starting point for more advanced optimisation.

3.  **`mlip_autopipec.modules.structure_generator.StructureGenerator`** (New Class)
    *   **Purpose**: To algorithmically generate the initial, diverse set of atomic structures, acting as the seed for the learning process.
    *   **Public API**:
        *   `__init__(config: FullConfig)`: The generator is initialised with the full system configuration, from which it will draw all necessary parameters.
        *   `generate_initial_structures() -> List[ase.Atoms]`: The main public method that returns a list of generated structures.
    *   **Key Internal Methods**:
        *   `_generate_sqs(composition: Dict, cell_size: int) -> ase.Atoms`: A method to generate Special Quasirandom Structures for alloys. This will likely be a well-encapsulated wrapper around a specialized third-party library like `icet`.
        *   `_generate_nms(atoms: ase.Atoms) -> List[ase.Atoms]`: A method for performing Normal Mode Sampling, primarily for molecular structures. This will involve calculating a Hessian matrix (potentially with a cheap classical potential) and displacing the atoms along the resulting normal modes.
        *   `_generate_airss(...)` and `_generate_rattle(...)`: These will be implemented as placeholder methods for future expansion to handle ionic and covalent materials, respectively.

4.  **`mlip_autopipec.orchestrator.Orchestrator`** (Updated)
    *   **New Method**: `run_initialization_phase(config_path: str)`: This method will replace the old entry point. Its responsibility is to orchestrate the entire setup phase: it will first call the `ConfigExpander`, save the resulting full configuration, then instantiate and call the `StructureGenerator`, and finally save the generated structures to the database.
    *   The old `run_cycle01_workflow` will be refactored into a more specific method like `run_labelling_phase`, which will be called by the main pipeline logic after the initialization is complete.

5.  **`mlip_autopipec.main.py`** (Updated)
    *   The `click` command will be updated to be the formal user entry point, accepting a file path to the `input.yaml`.
    *   The function signature will be updated to `@click.argument('config_file', type=click.Path(exists=True))`.
    *   The function body will now be responsible for instantiating the `Orchestrator` and calling the new `run_initialization_phase` with the provided config file path.

This updated design ensures a clean, logical, and user-centric flow. The configuration is handled comprehensively at the beginning of the process, creating an immutable and explicit set of parameters for the entire run. This is followed by the automated generation of the initial data. Only after this setup is complete does the system proceed to the labelling and training stages that were developed in the previous cycle.

## 4. Implementation Approach

The implementation of Cycle 02 will be approached methodically, focusing on building the configuration system first, as it provides the necessary parameters for all other new components. This will be followed by the development of the structure generation module and its final integration into the main orchestration workflow.

1.  **Pydantic Model Definition**:
    *   The first task is to define the `MinimalConfig` and `FullConfig` Pydantic models in `config/models.py`. We will start with a few key fields for each and then iteratively expand them as the implementation of the `ConfigExpander` and other modules progresses. This "schema-first" approach provides a clear, validated, and self-documenting structure for all configuration data from the outset.

2.  **Config Expander Development**:
    *   With the models defined, the `ConfigExpander` class will be created in `config/expander.py`.
    *   The heuristic logic will be implemented method by method, starting with the highest-level decision: `_determine_structure_type`. This function is the primary fork in the logic and will be implemented first.
    *   For `_get_dft_defaults`, we will create a simple, maintainable lookup table (e.g., a dictionary stored in a JSON or YAML file) that maps chemical elements to their recommended SSSP pseudopotential file names, cutoff energies, and other critical DFT parameters. This separates the heuristic data from the code.
    *   The main public `expand` method will be implemented last. It will orchestrate the calls to these internal, private methods, progressively building and finally returning the validated `FullConfig` object.

3.  **Structure Generator Development (Module A)**:
    *   The `StructureGenerator` class will be created in `modules/structure_generator.py`.
    *   The main `generate_initial_structures` method will be implemented with a primary conditional block (e.g., `if self.config.system.structure_type == 'alloy': ...`) that acts as a dispatcher, calling the appropriate private generation function based on the configuration.
    *   For the `_generate_sqs` method, we will integrate a well-tested third-party library like `icet`. The implementation will focus on writing a robust wrapper that handles the conversion from our internal data structures to the input format required by `icet`, and the subsequent conversion of its output back into a list of `ase.Atoms` objects.
    *   For `_generate_nms`, the initial implementation may use a simple, fast classical potential (like UFF or GFN-FF) to get an approximate Hessian, avoiding a costly DFT calculation at this early stage. The focus will be on correctly applying the displacements along the normal modes.

4.  **Orchestrator and CLI Integration**:
    *   The `Orchestrator` class will be refactored to include the new `run_initialization_phase` method, which will now serve as the primary entry point for a new run.
    *   The `main.py` CLI script will be updated to accept the `input.yaml` file path as a required argument.
    *   The logic within the main CLI function will be implemented to handle reading the YAML file (using a library like `PyYAML`), parsing its contents into the `MinimalConfig` Pydantic model for validation, and then passing this object to the orchestrator.
    *   The orchestrator will then be responsible for the sequence: call the expander, save the full config to `exec_config_dump.yaml` for user inspection and provenance, instantiate the Structure Generator with this new config, and finally save the generator's output to the database.

5.  **YAML Dependency Management**:
    *   The `pyyaml` library will be formally added to the project's dependencies in the `pyproject.toml` file to ensure it is available in the environment.

This detailed implementation plan ensures that the most user-facing and foundational parts of the system are built first, followed by the data generation logic that feeds into the existing core engine from Cycle 01, providing a logical, layered, and highly testable development progression.

## 5. Test Strategy

The test strategy for Cycle 02 is critically focused on validating the correctness of the new automated setup and data seeding phases. This involves extensive testing of the heuristic logic within the `ConfigExpander` and the structural outputs of the `StructureGenerator`. The strategy combines isolated unit tests for the internal logic with a comprehensive integration test for the new user-facing workflow.

**Unit Testing Approach (Min 600 words):**

The unit tests for this cycle will be designed to be highly specific, targeting the individual responsibilities of the new classes and their methods. All external library dependencies will be mocked to ensure the tests are fast and deterministic.

*   **`ConfigExpander`**: This component is perfectly suited for unit testing because its `expand` method is designed as a pure function: the same minimal configuration input should always produce the exact same full configuration output.
    *   A comprehensive test suite will be created with a wide variety of `MinimalConfig` inputs. This will include test cases for different material types: a binary alloy (`{"elements": ["Fe", "Pt"]}`), a simple molecule (`{"elements": ["H", "O"]}`), a covalent material (`{"elements": ["Si"]}`), and an ionic solid (`{"elements": ["Na", "Cl"]}`).
    *   For each test case, the test will execute the `expand` method and then perform a deep and detailed set of assertions on the output `FullConfig` object. For the FePt case, we will assert that `structure_type` was correctly identified as `'alloy'`, that the `dft_compute` section contains the correct SSSP pseudopotential file names for both Fe and Pt, that the `ecutwfc` is set to a physically sensible value (e.g., greater than a known minimum), and that the `magnetism` flag is appropriately set.
    *   We will also create specific tests for edge cases to ensure robustness. This includes providing an unsupported element in the `elements` list or a syntactically incorrect composition string. The tests will assert that the expander fails gracefully in these cases by raising a specific, informative exception, rather than producing an invalid configuration.

*   **`StructureGenerator`**: The tests for the generator will focus on ensuring that its methods produce structurally valid and diverse `ase.Atoms` objects.
    *   The main `generate_initial_structures` method will be tested with a mocked `FullConfig` object. By simply changing the `structure_type` field in this mock config between test runs, we can direct the method to call a specific internal generation function, allowing us to test the dispatcher logic.
    *   For the `_generate_sqs` method, the external `icet` library will be completely mocked. The purpose of this test is not to test `icet` itself, but our integration with it. The test will verify that our code correctly constructs the input objects required by the `icet` API and that it correctly parses the mock output from `icet` back into a list of `ase.Atoms` objects.
    *   For `_generate_nms`, we will provide a simple input molecule (like H2O) and assert that the output is a list of `ase.Atoms` objects. We will verify that each object in the list has the same number and type of atoms as the input, and that their positions are slightly and uniquely perturbed, consistent with the concept of atomic displacements along normal modes.

**Integration Testing Approach (Min 300 words):**

The primary integration test for Cycle 02 is crucial. It will verify the entire new workflow, starting from the user's command-line invocation and ending with the successful population of the database with initial structures. This test connects the CLI, the `ConfigExpander`, the `StructureGenerator`, and the `AseDB` in a single, representative test case.

The test setup will be carefully orchestrated to be self-contained and reproducible:
1.  A temporary directory will be created at the start of the test run to isolate all file I/O.
2.  Inside this directory, a simple `input.yaml` file will be created programmatically (e.g., for elemental Silicon, "Si").
3.  The test will then use the `click.testing.CliRunner` to invoke the main `run` command defined in `main.py`, pointing it to the newly created temporary `input.yaml` file.

After invoking the command, the test will perform a series of verifications to confirm the success of the entire initialization phase:
1.  **Configuration Expansion Verification**: The test will first assert that a new file named `exec_config_dump.yaml` was created in the temporary directory. It will then load the content of this file, parse it using `PyYAML`, and validate it against the `FullConfig` Pydantic model to ensure it is a valid and complete configuration file. It will also perform spot-checks on the content, asserting that it is consistent with the "Si" input (e.g., `structure_type` should be 'covalent').
2.  **Structure Generation Call**: To keep the test fast and independent of complex scientific libraries, the `StructureGenerator`'s actual generation methods (like `_generate_sqs`) can be mocked. The primary goal here is to verify that the `generate_initial_structures` method was called by the `Orchestrator`, confirming that the control flow is correct.
3.  **Database State Verification**: The test will instantiate the `AseDB` to use a temporary database file within the isolated directory. After the CLI command has finished, the test will connect to this database. It will then execute a query and assert that a number of new structures for Silicon have been successfully added to the database, and that their `state` is correctly set to 'initial_generated'.

This single integration test provides a very high degree of end-to-end confidence in the new initialization phase. It confirms that a user can successfully launch a new project, have the configuration automatically and correctly determined, and see the database populated with the first set of structures, all without any manual intervention beyond the creation of the initial, minimal input file.
