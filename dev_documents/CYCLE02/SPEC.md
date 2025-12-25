# Cycle 2 Specification: Structure Generation and Configuration System

**Version:** 1.0.0
**Status:** Final

## 1. Summary

This document provides the detailed technical specification for Cycle 2 of the MLIP-AutoPipe project. Building directly upon the core labelling and training engine established in Cycle 1, this cycle introduces the first layer of true automation and user-centric design. The primary objective is to liberate the user from the burden of data preparation, a critical step towards fulfilling the project's core philosophy of "removing the human expert from the loop." This cycle addresses two major components that work in tandem: first, the implementation of **Module A (Structure Generator)**, a sophisticated suite of algorithms designed to produce physically plausible initial training structures from scratch. Second, the development of the **Two-Tier Configuration System**, an intelligent front-end that translates a minimal, human-friendly user input into the exhaustive, explicit configuration required by the pipeline's backend.

The successful completion of Cycle 2 will mark a significant milestone in the project's usability. The system will no longer require the user to provide a directory of pre-made, expertly curated atomic structures. Instead, the user will interact with the system via a simple `input.yaml` file, specifying little more than the material's chemical composition and the desired simulation conditions. The `ConfigExpander`, the heuristic engine at the heart of the two-tier system, will parse this minimal input and automatically deduce the necessary details. It will classify the material type (e.g., alloy, ionic solid), select the most appropriate strategy for initial structure generation, and populate all the necessary DFT and simulation parameters with sensible, physically-grounded defaults. This expanded configuration will then drive Module A, the `StructureGenerator`, which will execute the selected strategy—be it generating Special Quasirandom Structures (SQS) for a random alloy, applying Normal Mode Sampling (NMS) to a molecule, or running a simulated Melt-Quench protocol for a covalent material—to produce a diverse and high-quality initial dataset. This dataset is then passed seamlessly to the core engine developed in Cycle 1 for labelling and training. This cycle effectively builds the "brain" and "hands" of the initial data generation process, transforming the pipeline from a powerful but manual tool into a genuinely automated, user-friendly solution.

## 2. System Architecture

The architecture of Cycle 2 extends the linear pipeline of Cycle 1 by prepending two new crucial stages: the Config Expander and Module A. This shifts the user's entry point from providing low-level structural data to providing a high-level, declarative description of their goal. The system architecture now reflects a more complete "input-to-model" workflow, where the initial stages are dedicated to interpreting the user's intent and automatically preparing the necessary data for the core engine. The `exec_config_dump.yaml` file plays a pivotal role, acting as the bridge between the high-level user input and the low-level execution details, while also serving as a critical artifact for scientific reproducibility.

```mermaid
graph TD
    A[User Input: input.yaml] --> B{Config Expander (Heuristic Engine)};
    B -- Analyzes Input & Applies Heuristics --> C[Full Config: exec_config_dump.yaml];
    C -- Is Saved for Provenance --> D((Output Directory));
    C -- Configures & Drives --> E[Module A: Structure Generator];
    E -- Generates Diverse Structures --> F[List of Initial Structures];
    F --> G[Module C: Labelling Engine];
    C -- Also Configures --> G;
    G --> H[ASE Database];
    H --> I[Module D: Training Engine];
    C -- Also Configures --> I;
    I --> J[Trained MLIP Model];

    subgraph "New Components in Cycle 2"
        direction LR
        B;
        E;
    end
```

**Detailed Workflow Description:**

1.  **High-Level User Input:** The process now begins with a user creating a purposefully minimal `input.yaml` file. For instance, to generate a potential for amorphous silicon, the user might only provide:
    ```yaml
    system:
      composition: "Si"
      tags: ["amorphous", "covalent"]
    simulation:
      temperature: [300, 1500]
    ```
2.  **Configuration Expansion:** The pipeline is launched with this `input.yaml`. The very first step is the invocation of the `ConfigExpander`. This heuristic engine performs a sequence of intelligent operations:
    a.  **Parsing and Validation:** It first parses the YAML and validates it against a strict Pydantic model (`MinimalConfig`) to ensure the user has provided the necessary information and used correct key names.
    b.  **Material Analysis:** It analyzes the `system` block. It identifies the elements, and based on composition and optional tags, it determines the material's fundamental type (e.g., 'covalent').
    c.  **Heuristic Decision-Making:** Based on this analysis, it makes a series of expert-level decisions. It selects the most appropriate structure generation strategy from Module A (e.g., 'melt_quench' for amorphous Si). It infers sensible defaults for all unspecified parameters, such as the number of atoms for the simulation cell, the parameters for the melt-quench protocol (e.g., peak temperature, cooling rate), and all the required DFT parameters from the SSSP protocol, as in Cycle 1.
    d.  **Full Configuration Generation:** It populates a comprehensive Pydantic model (`FullConfig`) with all these explicit parameters. This complete configuration is then serialized into the `exec_config_dump.yaml` file and saved to the output directory, creating an immediate and permanent record of the run's exact settings.
3.  **Automated Structure Generation (Module A):** The `StructureGenerator` is then instantiated, its behavior dictated entirely by the contents of the `exec_config_dump.yaml`. It reads the selected strategy ('melt_quench') and its associated parameters. It then executes this strategy, which might involve creating an initial crystal, running a short, high-temperature MD simulation with a classical potential to melt it, and then rapidly cooling it to generate a variety of amorphous configurations. The output is a list of `ASE.Atoms` objects, representing the diverse initial dataset.
4.  **Downstream Core Pipeline:** From this point forward, the workflow is identical to that of Cycle 1. The list of newly generated structures is handed off to the `LabellingEngine` (Module C), which calculates their DFT properties. The results are stored in the ASE database, and finally, the `TrainingEngine` (Module D) is called to produce the trained MLIP model. The key difference is that the entire process was initiated from a simple, high-level request, with the crucial and complex data generation step handled completely automatically.

## 3. Design Architecture

Cycle 2 introduces two major new software components, the `ConfigExpander` and the `StructureGenerator`, and necessitates the formalization of the configuration system using Pydantic data models. The design focuses on creating a clear separation between interpreting user intent (config expansion) and acting on that intent (structure generation).

**Key Classes and Modules:**

-   **`src/mlip_autopipe/config/expander.py`:**
    -   **`ConfigExpander` class:** This is the heuristic engine, the "brain" of the automated setup process.
        -   `__init__(self, minimal_config: dict)`: Takes the raw dictionary parsed from the user's `input.yaml`.
        -   `expand(self) -> FullConfig`: The main public method that orchestrates the expansion process and returns a fully validated `FullConfig` Pydantic object.
        -   `_determine_bond_type(self, composition: str, tags: list) -> str`: A key private method that uses a combination of chemical heuristics (leveraging libraries like `pymaten`) and user-provided tags to classify the material into 'alloy', 'ionic', 'covalent', or 'molecular'. This classification is the primary driver for many subsequent decisions.
        -   `_select_generation_strategy(self, bond_type: str) -> dict`: Based on the bond type, this method returns a dictionary containing the name of the generation algorithm to use (e.g., `sqs`) and a set of sensible default parameters for that algorithm.
        -   `_get_dft_defaults(self, elements: list) -> dict`: This method, inherited and expanded from Cycle 1, is now a formal part of the heuristic engine, responsible for populating the entire `dft_compute` section of the full configuration.

-   **`src/mlip_autopipe/config/models.py`:** This new file is central to Cycle 2.
    -   **`MinimalConfig(BaseModel)`:** A Pydantic model that defines the limited, user-facing schema for `input.yaml`. It will use `ConfigDict(extra='forbid')` to strictly prevent users from providing unknown or misspelled configuration keys, providing immediate and clear feedback on typos. The fields will be simple and declarative (e.g., `composition: str`).
    -   **`FullConfig(BaseModel)`:** A comprehensive, nested Pydantic model that defines the entire set of explicit parameters required for a run. This model serves as the single source of truth for the rest of the application. Its strictness ensures that no part of the pipeline can be run with an incomplete or invalid configuration.

-   **`src/mlip_autopipe/modules/structure_generator.py`:**
    -   **`StructureGenerator` class:** This class acts as a factory and dispatcher for various structure generation algorithms.
        -   `__init__(self, config: dict)`: Takes the `structure_generation` section of the `FullConfig` as its input.
        -   `generate(self) -> List[Atoms]`: The main public method. It inspects its configuration to determine which generation algorithm was selected (e.g., `sqs`) and then calls the corresponding private method, passing along the required parameters.
        -   `_generate_sqs(self, ...)`: A private method that implements the SQS generation for alloys. This will likely involve using a dedicated external library like `icet`. The method will be responsible for preparing the input for this library and parsing its output into a list of `ASE.Atoms` objects.
        -   `_generate_nms(self, ...)`: A private method for Normal Mode Sampling. This will involve calculating the Hessian and vibrational modes of an input molecule (using ASE or `pymaten`) and then creating new structures by displacing atoms along these modes.
        -   `_generate_airss(self, ...)`: A private method implementing a simplified random structure search, useful for discovering polymorphs of ionic crystals.
        -   `_generate_rattle_and_quench(self, ...)`: A private method for generating amorphous or defected structures. This method will require a simple, pre-existing classical potential (like Lennard-Jones) to run the fast MD simulations, and the details of this potential must also be part of its configuration.

## 4. Implementation Approach

The implementation of Cycle 2 will be tackled by first defining the data contracts (the Pydantic models) and then building the two main components, the `ConfigExpander` and `StructureGenerator`, in parallel, with continuous unit testing.

1.  **Data Contract Definition:** The first and most critical step is the implementation of the `MinimalConfig` and `FullConfig` Pydantic models in `src/mlip_autopipe/config/models.py`. This step is foundational because these models define the precise inputs and outputs for the `ConfigExpander` and the configuration required by all other modules. Getting this data architecture right is essential for a clean and maintainable codebase.

2.  **Configuration Expander Implementation:**
    a.  The `ConfigExpander` class will be created. Its development will start with the `_determine_bond_type` method, as this is the primary decision point. This will necessitate adding `pymaten` as a new project dependency in `pyproject.toml`.
    b.  The subsequent heuristic methods (`_select_generation_strategy`, `_get_dft_defaults`, etc.) will be implemented one by one. Each method will be responsible for populating a specific section of the `FullConfig` object.
    c.  A comprehensive suite of unit tests will be developed in parallel. We will create a set of diverse `input.yaml` examples and write tests that run the `ConfigExpander` on each one, asserting that the resulting `FullConfig` object is populated with the correct, expected values. This test-driven approach is crucial for ensuring the reliability of the complex heuristic logic.

3.  **Structure Generator Implementation:**
    a.  The `StructureGenerator` class will be implemented, starting with the main `generate` dispatch method.
    b.  Each of the private generation methods (`_generate_sqs`, `_generate_nms`, etc.) will be implemented in turn. For methods that wrap external tools, the implementation will focus on robustly managing the subprocess, preparing its inputs, and parsing its outputs. For native methods like `_generate_rattle_and_quench`, the logic will be implemented directly using the ASE library's MD capabilities.
    c.  Each generation method will be accompanied by its own dedicated unit tests. For example, the SQS test will assert that the generated structure has the correct stoichiometry and cell size. The NMS test will assert that the generated structures correctly represent displacements along calculated vibrational modes. These tests are essential for ensuring the scientific validity of the generated data.

4.  **Integration into the Main Workflow:**
    a.  The main CLI entry point, defined in Cycle 1, will be modified. It will no longer take a directory of structures as input, but rather the path to the minimal `input.yaml` file.
    b.  The main workflow orchestrator will be updated to reflect the new architecture. It will first instantiate and run the `ConfigExpander`. It will then save the resulting `FullConfig` to disk as `exec_config_dump.yaml`. This same `FullConfig` object will then be used to initialize and run the `StructureGenerator`.
    c.  The list of `ASE.Atoms` objects produced by the `StructureGenerator` will then be passed as input to the `LabellingEngine` from Cycle 1, seamlessly connecting the new front-end to the existing core engine.
    d.  A new end-to-end integration test will be created. This test will be crucial for verifying that all components, new and old, work together. It will start with a minimal `input.yaml` for a simple alloy (e.g., CuAu), and will assert that the entire pipeline runs without error, producing a final, trained MLIP model file.

## 5. Test Strategy

The test strategy for Cycle 2 is focused on verifying the correctness of the new, complex heuristic logic and ensuring the seamless integration of the new modules with the existing core pipeline.

**Unit Testing Approach:**
Unit tests for this cycle are critical for ensuring the reliability of the automated decision-making processes.
-   **Config Expander:** The `ConfigExpander` will be the most heavily unit-tested component of this cycle. Its logic is complex and full of decision points, making it a prime candidate for bugs. We will create a dedicated test suite with a variety of sample `input.yaml` files, each designed to test a specific heuristic path. For example, we will have `Si_covalent.yaml`, `FePt_alloy.yaml`, `NaCl_ionic.yaml`, and `H2O_molecular.yaml`. For each of these inputs, a test will run the `expand()` method and then perform a deep assertion on the resulting `FullConfig` object. For the `NaCl_ionic.yaml` case, the test will assert that `bond_type` is correctly identified as 'ionic', that the selected generation strategy is 'airss', and that the `dft_compute` section contains the correct SSSP parameters for both Sodium and Chlorine. We will also explicitly test failure modes. For instance, we will test an `input.yaml` with a misspelled key and assert that the Pydantic validation cleanly raises a `ValidationError` with a user-friendly message.

-   **Structure Generator:** Each of the generation algorithms within the `StructureGenerator` will be tested in isolation to verify its scientific correctness. For the `_generate_sqs` method, the unit test will provide a simple composition (e.g., A0.5B0.5) and a supercell size, and then assert that the returned `ASE.Atoms` object contains the correct total number of atoms and the correct number of atoms of type A and B. For the `_generate_nms` method, the test will start with a known, simple molecule like water. It will then call the method and assert that the returned list of structures includes displacements that correspond to the known symmetric stretch, asymmetric stretch, and bending modes of the water molecule. These tests are crucial for ensuring that the initial data provided to the DFT engine is physically meaningful.

**Integration Testing Approach:**
Integration tests will verify that the newly created modules correctly communicate with each other and with the existing pipeline.
-   **Config to Generation Flow:** A key integration test will be designed to verify the link between the `ConfigExpander` and the `StructureGenerator`. This test will start with a minimal `input.yaml` for an alloy. It will then run the `ConfigExpander` to produce the `FullConfig`. This `FullConfig` object will then be passed to the `StructureGenerator`. The test will use `pytest-mock` to "spy" on the `StructureGenerator` instance, mocking out the actual (and potentially slow) `_generate_sqs` method. The assertion will be that the `_generate_sqs` method was called exactly once, confirming that the dispatcher logic correctly interpreted the configuration and invoked the appropriate algorithm.

-   **Full Pipeline Smoke Test:** To ensure that the new front-end integrates seamlessly with the Cycle 1 backend, a broader "smoke test" will be implemented. This test provides the highest level of confidence that the overall system is functional. It will use a real, minimal `input.yaml` for a very simple system like crystalline Carbon. It will then execute the entire pipeline from start to finish: `ConfigExpander` -> `StructureGenerator` -> `LabellingEngine` -> `TrainingEngine`. To ensure this test runs quickly enough for a CI environment, the `LabellingEngine` will be configured to use a mock "dummy" DFT code—a simple script that takes a QE input file and immediately returns a predefined, correctly formatted output file without performing any real calculation. The success criterion for the test is simple: the entire pipeline must complete without any unhandled exceptions, and a final model file must be created. This test confirms that the data—from the Pydantic config object to the list of `Atoms` to the labelled database—flows correctly through all the integrated components.
