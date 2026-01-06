# SPEC: CYCLE 01 - Core Engine and CLI

## 1. Summary

This document provides the detailed technical specification for the first development cycle of MLIP-AutoPipe. The primary objective of Cycle 01 is to build the foundational backend and command-line interface (CLI). This cycle will deliver a fully functional, end-to-end pipeline capable of generating a database of atomic structures for alloy-type materials. The focus is on creating a robust, modular, and scientifically valid architecture. Key deliverables include the `PipelineRunner` orchestrator, the `AlloyGenerator` with built-in physical constraint checking, a core Molecular Dynamics (MD) engine, a `RandomSampler`, and an `AseDBWrapper`.

## 2. System Architecture

The architecture for Cycle 01 is a linear, four-stage pipeline executed from the command line.

**File Structure (Cycle 01 Focus):**

Files/directories to be created or modified are in **bold**.

```
.
├── **pyproject.toml**
└── src
    └── mlip_autopipec
        ├── **__init__.py**
        ├── **cli.py**
        ├── **config.py**
        ├── **errors.py**         # Custom exception types
        ├── pipeline
        │   └── **runner.py**
        ├── generators
        │   ├── **base.py**
        │   └── **alloy.py**
        ├── explorers
        │   └── **md_engine.py**
        ├── sampling
        │   ├── **base.py**
        │   └── **random.py**
        └── storage
            └── **db_wrapper.py**
```

## 3. Implementation Tasks with Scientific Validation

1.  **Setup `pyproject.toml`**: Add dependencies: `typer`, `pydantic`, `pyyaml`, `ase`, `mace-torch`.
2.  **Implement `config.py`**: Define all Pydantic models. Include validators to ensure scientific sense (e.g., composition fractions must sum to 1.0).
3.  **Implement `errors.py`**: Define custom exceptions like `PhysicsViolationError` to be used for scientific constraint checks.
4.  **Implement `storage/db_wrapper.py`**: Create the `AseDBWrapper` class.
5.  **Implement `generators/base.py` and `generators/alloy.py`**:
    *   **Task 5a (Interface):** Define the `BaseStructureGenerator` ABC.
    *   **Task 5b (Implementation):** Implement the `AlloyGenerator`.
    *   **Task 5c (Scientific Validation):** Inside the `AlloyGenerator`, implement a **strict atomic overlap check**. This method will calculate the distance between all pairs of atoms in a generated structure and raise a `PhysicsViolationError` if any distance is less than a configurable threshold (e.g., 0.7 Å). This is a critical task to ensure only physically plausible structures are passed to the next stage.
6.  **Implement `explorers/md_engine.py`**: Implement the core parallel MD logic.
7.  **Implement `sampling/base.py` and `sampling/random.py`**: Define the `BaseSampler` ABC and implement the `RandomSampler`.
8.  **Implement `pipeline/runner.py`**: Implement the `PipelineRunner` to orchestrate the workflow.
9.  **Implement `cli.py`**: Create the `typer`-based user entrypoint.
10. **Write Tests**: Concurrently write unit and integration tests as per the strategy below.

## 4. Detailed Test Strategy

Testing for Cycle 01 will be rigorous, with a strong emphasis on validating the scientific correctness of the outputs.

**Unit Testing Approach:**
*   **`config.py`**: Test that Pydantic models raise `ValidationError` for scientifically invalid configurations (e.g., negative temperature, composition fractions not summing to 1.0).
*   **`generators/alloy.py`**:
    *   **Test 1 (Composition):** Assert that the chemical formula of every generated `ase.Atoms` object is correct.
    *   **Test 2 (Physical Plausibility):** Write a test that programmatically creates a structure with two overlapping atoms. Assert that the generator's validation method correctly raises a `PhysicsViolationError` for this structure. This directly tests the implementation of **Task 5c**.
*   **`storage/db_wrapper.py`**:
    *   **Test 3 (Data Integrity):** Write `ase.Atoms` objects with known energy and forces to a temporary database. Read the data back and assert that the retrieved physical properties are identical to the originals.

**Integration Testing Approach:**
The primary integration test will validate the entire pipeline's data flow and scientific coherence.
*   **End-to-End Scientific Verification Test**:
    1.  A simple YAML configuration will be created for a binary system (e.g., Cu2).
    2.  A fast, deterministic classical potential (like Lennard-Jones from ASE) will be used to make the test fast and reproducible.
    3.  The `typer.testing.CliRunner` will execute the pipeline.
    4.  **Post-run assertions will validate the scientific output:**
        *   The process must exit with code 0.
        *   The final database must be created.
        *   The number of structures in the database must match the `num_samples` setting.
        *   **A query will be run on the database to check that the chemical composition of every structure is correct.**
        *   **The potential energies of the structures in the database will be checked to ensure they are within a physically reasonable range for the simple potential used.**
