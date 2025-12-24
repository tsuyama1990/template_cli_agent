# Cycle 1 Specification: The Core Engine

## 1. Summary

This first development cycle marks the establishment of the foundational components of the MLIP-AutoPipe system. The primary objective of Cycle 1 is to create a functional, end-to-end workflow that can take a predefined set of atomic structures, calculate their physical properties using Density Functional Theory (DFT), and then use this data to train a Machine Learning Interatomic Potential (MLIP). This cycle focuses exclusively on building the backend computational machinery, specifically the **Labeling Engine (Module C)** and the **Training Engine (Module D)**.

At the completion of this cycle, we will have a system capable of producing an MLIP from manually provided inputs. The workflow will not yet be autonomous; it will require a user to supply the initial structures and to manually trigger the calculation and training steps. However, this cycle is critical as it lays the groundwork for all future automation. It will validate our core architectural choices for interacting with external simulation codes and for integrating MLIP training frameworks. Key achievements will include the successful automation of Quantum Espresso (QE) calculations, ensuring robust parameterisation and error handling, and the implementation of a Delta Learning training strategy to produce physically realistic potentials. This cycle also involves setting up the modern Python project environment using `pyproject.toml` and the `uv` package manager, which is essential for ensuring reproducibility and efficient dependency management throughout the project's lifecycle.

## 2. System Architecture

The architecture for Cycle 1 is intentionally minimal, focusing only on the two modules that form the computational core of the entire pipeline. These modules are designed to be robust and self-contained, with clear inputs and outputs, allowing them to be orchestrated by the workflow manager in later cycles.

**Module C: Labeling Engine (Automated DFT)**
The Labeling Engine acts as a sophisticated wrapper around the Quantum Espresso (QE) DFT code. Its primary responsibility is to abstract away the complexity of running DFT calculations. It receives a specific atomic structure and returns the same structure annotated with its DFT-calculated energy, forces, and stress tensor. The module will automatically handle the generation of QE input files, applying best-practice settings based on a given configuration. For this cycle, the configuration will be provided manually, but the logic is designed to be driven by the `exec_config_dump.yaml` in the future. A key feature is the adoption of the SSSP (Standard Solid State Pseudopotentials) protocol, which provides standardised, high-precision settings for pseudopotentials and cutoff energies, thereby removing the need for expert-level manual tuning. The engine will also include basic error recovery logic, for instance, attempting to restart a calculation with adjusted parameters if the initial self-consistent field (SCF) calculation fails to converge.

**Module D: Training Engine (Delta Learning)**
The Training Engine is responsible for constructing the MLIP from the labeled data produced by Module C. A central feature of this module is the implementation of **Delta Learning**. Instead of learning the absolute DFT energies and forces directly, the model learns the *difference* (the delta) between the DFT values and the predictions of a simpler, physics-based reference potential (such as a Ziegler-Barlit-Biersack potential for short-range repulsion). This approach is highly beneficial because it embeds fundamental physical constraints directly into the model. The reference potential ensures correct physical behaviour, like the strong repulsive force when atoms are too close, even in regions where training data is sparse. This leads to more data-efficient training and results in MLIPs that are more stable and reliable during simulations. The engine will be designed to work with a modern MLIP framework, such as Atomic Cluster Expansion (ACE), which offers a good balance of accuracy and computational performance.

## 3. Design Architecture

The software design for this cycle focuses on creating two primary classes, one for each module, with clear and simple interfaces.

**Class: `LabelingEngine`**
This class will manage all interactions with the Quantum Espresso executable.
-   **`__init__(self, config: dict)`**: The constructor will take a configuration dictionary containing paths to executables and DFT parameters (e.g., k-point density, smearing settings).
-   **`run_calculation(self, atoms: ase.Atoms) -> ase.Atoms`**: This is the main public method. It takes an `ase.Atoms` object as input. Internally, it will call private helper methods to:
    1.  `_generate_pw_input()`: Create the `pw.in` text file from the `Atoms` object and configuration parameters.
    2.  `_execute_qe()`: Run the `pw.x` executable as a subprocess, capturing its output and errors.
    3.  `_parse_pw_output()`: Read the QE output file to extract the final energy, forces on each atom, and the system's stress tensor.
    The method will return the original `Atoms` object, but with the results added to its `info` (for energy/stress) and `arrays` (for forces) dictionaries.

**Class: `TrainingEngine`**
This class will manage the MLIP training process.
-   **`__init__(self, config: dict)`**: The constructor will take a configuration dictionary specifying the MLIP model type, its hyperparameters (e.g., cutoff radius), and the choice of reference potential for Delta Learning.
-   **`train(self, atoms_list: list[ase.Atoms]) -> Any`**: This is the main public method. It accepts a list of labeled `ase.Atoms` objects. Its internal workflow will be:
    1.  For each `Atoms` object, calculate the energy and forces using the baseline reference potential.
    2.  Subtract these baseline values from the DFT labels to get the "delta" values.
    3.  Prepare the training dataset in the format required by the chosen MLIP framework (e.g., ACE).
    4.  Fit the MLIP model to the delta values.
    5.  The method will return the trained model artifact, which could be a file path or an in-memory object, ready for use in simulations.

**Data Flow:**
The expected data flow for this cycle is straightforward: A user prepares a list of `ase.Atoms` objects. This list is iterated over, and each object is passed to an instance of `LabelingEngine`. The returned, labeled objects are collected into a new list, which is then passed to an instance of `TrainingEngine` to produce the final MLIP.

## 4. Implementation Approach

The implementation will be carried out in a series of logical steps to ensure a solid foundation.

1.  **Project Setup:** Initialise the project using `git`. Create the `pyproject.toml` file to define project metadata and dependencies. Key dependencies will include `ase`, `numpy`, a DFT parsing library (if available), and the chosen MLIP framework. Set up a virtual environment using `uv`.
2.  **`LabelingEngine` Shell:** Create the `c_labeling_engine.py` file and define the `LabelingEngine` class structure with placeholder methods.
3.  **QE Input Generation:** Implement the `_generate_pw_input` method. This function will be purely programmatic, translating the atom positions, cell dimensions, and DFT parameters into the specific NAMELIST format required by Quantum Espresso. This is a critical step and will be thoroughly tested.
4.  **Subprocess Execution:** Implement the `_execute_qe` method to call `pw.x` using Python's `subprocess.run`. The implementation must handle command execution, wait for completion, and check for non-zero exit codes to detect errors.
5.  **QE Output Parsing:** Implement the `_parse_pw_output` method. This will involve robustly parsing the plain-text output file from QE to find and extract the specific lines containing the final energy, forces, and stress. Regular expressions or careful string searching will be used here.
6.  **`TrainingEngine` Shell:** Create the `d_training_engine.py` file and define the `TrainingEngine` class structure.
7.  **Delta Learning Logic:** Implement the core Delta Learning feature. This involves selecting a reference potential (e.g., from within ASE's calculator library) and writing the code to compute its contribution and subtract it from the DFT labels.
8.  **MLIP Framework Integration:** Write the code to interface with the chosen MLIP library. This will involve converting the prepared delta-learning dataset into the library's expected format and calling its fitting/training function.
9.  **Integration Script:** Create a simple Python script (`run_cycle1_manual.py`) that demonstrates the complete workflow. This script will create a sample `Atoms` object, instantiate the engines, and call their methods in the correct sequence, printing the final result.

## 5. Test Strategy

Testing in Cycle 1 is crucial to verify the correctness of our interactions with the external DFT and MLIP libraries.

**Unit Testing Approach (Min 300 words):**
We will use `pytest` for all unit testing. The `LabelingEngine` will be tested extensively without needing to run an actual DFT calculation. We will use `unittest.mock.patch` to mock the `subprocess.run` call. This allows us to isolate the input generation and output parsing logic. One set of tests will focus on `_generate_pw_input`. We will create an `ase.Atoms` object for a simple material like a silicon dimer, and assert that the generated `pw.in` string exactly matches a pre-written, known-good input file. We will test this for various settings (e.g., with and without spin polarisation). Another set of tests will target `_parse_pw_output`. We will have a collection of saved QE output files (including ones with errors) and write tests that assert our parsing function correctly extracts the numerical values or raises the appropriate exceptions.

For the `TrainingEngine`, unit tests will verify the Delta Learning logic. We will create a mock dataset of labeled `Atoms` objects and a mock reference potential. The tests will assert that the data passed to the underlying MLIP framework's fitting function is correctly transformed (i.e., the reference potential's contribution has been subtracted). We will also have a simple test that confirms the `train` method runs without crashing and returns an object of the expected type when given a small, valid dataset.

**Integration Testing Approach (Min 300 words):**
Integration tests will verify that the data flows correctly between the `LabelingEngine` and the `TrainingEngine`. The main integration test will not perform a full DFT calculation but will use a mocked `LabelingEngine`. This mock engine's `run_calculation` method will be programmed to return a pre-computed, correctly labeled `Atoms` object. The test will then pass this object to the `TrainingEngine`. The primary assertion will be to check for data compatibility. We will verify that the data structures and units (e.g., for energy and forces) produced by the `LabelingEngine` are exactly what the `TrainingEngine` expects. For example, we will confirm that the keys used to store energy in `atoms.info` and forces in `atoms.arrays` are consistent between the two modules. This test ensures that the "contract" between the two modules is respected. A second, more comprehensive integration test will be a full end-to-end run on a tiny system (e.g., H2 molecule) that executes a real (but very fast) QE calculation and then trains a simple MLIP. This will be a slow test, marked as such in the test suite, but it will provide the ultimate confirmation that the two modules work together correctly in a real-world scenario.
