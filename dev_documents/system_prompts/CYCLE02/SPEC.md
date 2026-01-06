# Specification: MLIP-AutoPipe - Cycle 2

## 1. Summary

This document provides the detailed technical specification for the second and final implementation cycle of the MLIP-AutoPipe project. This cycle represents the culmination of the project's vision, building upon the robust architectural foundation established in Cycle 1 to implement the advanced, scientifically complex components that constitute the core value of the tool. The primary objective is to transform the application from a basic architectural skeleton into a powerful and sophisticated data generation engine capable of producing high-quality, diverse datasets for training state-of-the-art Machine Learning Interatomic Potentials. Where Cycle 1 was about building a reliable "chassis," Cycle 2 is about designing and installing the high-performance "engine." The key deliverables for this cycle are a feature-complete `LabelingEngine` capable of interfacing with external Density Functional Theory (DFT) codes, the development of the sophisticated hybrid Molecular Dynamics / Monte Carlo (MD/MC) exploration engine, the integration of the Farthest Point Sampling (FPS) algorithm for intelligent and diverse data selection, the expansion of the structure generator library to support a wider range of material systems, and the creation of a user-friendly web-based UI for interactive workflow configuration.

The implementation of the `LabelingEngine` will be a critical task, involving the replacement of the mock engine from Cycle 1 with a production-ready version. This new engine will be designed to robustly manage external simulation processes (like Quantum Espresso or VASP) using Python's `subprocess` module. It will be responsible for the entire lifecycle of an external calculation: dynamically generating formatted input files based on the atomic structure and configuration, executing the external code in a managed subprocess, and, most importantly, parsing the complex, text-based output files to reliably extract the key physical quantities—energy, forces, and stress. This component must be built with resilience in mind, capable of handling common failure modes of DFT calculations. The MD/MC exploration engine will be the most scientifically complex component developed in this cycle. It will integrate a standard ASE-based MD simulator with custom Monte Carlo moves (e.g., atom swaps) to enable efficient exploration of both configurational and compositional phase space. A key feature will be the automatic thermodynamic ensemble switching (NVT for surfaces, NPT for bulk), a critical detail for ensuring physically meaningful simulations. Finally, the sampling module will be upgraded from a simple random selection to include FPS, which leverages SOAP (Smooth Overlap of Atomic Positions) descriptors to select a maximally diverse subset of structures from the potentially vast simulation trajectories. By the successful completion of this cycle, MLIP-AutoPipe will be a feature-complete application, delivering on its promise to automate and elevate the process of MLIP training data generation.

## 2. System Architecture

Cycle 2 focuses on enriching the existing, stable architecture with the advanced scientific capabilities that were previously represented by placeholders. The architectural changes are primarily additive, involving the replacement of mock components from Cycle 1 with functional, feature-rich implementations and the introduction of entirely new modules for advanced functionality. The modular design established in Cycle 1 is what makes this seamless upgrade possible; the core `WorkflowOrchestrator` will be modified to accommodate the new stages, but the fundamental data flow and separation of concerns remain unchanged. The system's robustness is enhanced by ensuring that each new component adheres to the established principles of clear interfaces and schema-driven data contracts. The introduction of a `sampling` module and a `webui` directory represents a significant expansion of the tool's capabilities, moving it from a simple linear pipeline to a more sophisticated workflow with intelligent data reduction steps and an alternative, user-friendly interface. This evolution reflects a mature software architecture that can accommodate significant new functionality without requiring a foundational redesign.

**File Structure (Cycle 2 Focus):**
The files and directories to be created or significantly modified in this cycle are marked in **bold**. This file tree is the exact blueprint for the final state of the application's source code.

```
.
└── src
    └── mlip_autopipec
        ├── **cli.py**              # **Modified**: Add a new command 'webui' to launch the web interface.
        ├── **config.py**           # **Modified**: Add Pydantic models for Exploration, FPS, and DFT settings.
        ├── core
        │   ├── ase_db_wrapper.py # No changes. The interface from Cycle 1 is stable.
        │   ├── **generators.py**     # **Modified**: Add new generator classes (e.g., IonicGenerator).
        │   ├── **labeling.py**       # **Modified**: Replace MockLabelingEngine with a real DFTLabelingEngine.
        │   ├── **exploration.py**    # **NEW**: A new module for the MD/MC exploration engine.
        │   ├── **sampling.py**       # **NEW**: A new module for sampling algorithms (Random and FPS).
        │   └── **workflow.py**       # **Modified**: Integrate the new Exploration and Sampling stages.
        └── **webui**
            ├── **__init__.py**     # **NEW**
            └── **app.py**          # **NEW**: The Streamlit or FastAPI application for the web UI.
```

The architectural evolution is strategic and layered:
-   `config.py`: The application's "public API" will be expanded with new Pydantic models (`ExplorationConfig`, `SamplingConfig`, `DFTConfig`) to provide users with fine-grained control over the new, complex stages. This maintains the principle of a single, validated source of truth for all settings.
-   `generators.py`: The library of available structure generators will be expanded beyond the basic `AlloyGenerator` to support a wider range of material systems (e.g., ionic crystals, surfaces). This demonstrates the extensibility of the factory pattern established in Cycle 1.
-   `labeling.py`: The simple mock engine will be replaced by a `DFTLabelingEngine`. This is a critical upgrade, transforming the pipeline from a demonstration into a scientifically useful tool. It will encapsulate all the complexity of interacting with external, non-Python codes.
-   `exploration.py` & `sampling.py`: These new modules will house the core scientific logic of the application. Their introduction as separate modules maintains the high cohesion of the architecture.
-   `workflow.py`: The `WorkflowOrchestrator` will be the central point of integration. It will be significantly updated to manage the new, full pipeline sequence: Generation -> Exploration -> Sampling -> Labeling -> Storage.
-   `webui/app.py`: This new user-facing component will provide an alternative, interactive entry point to the application's core logic, reusing the same `WorkflowOrchestrator` to ensure consistency.

## 3. Design Architecture

The schema-first design philosophy that was foundational to Cycle 1 will be rigorously extended in Cycle 2 to encompass the new, advanced features. Every new piece of user-configurable functionality will be represented by a corresponding Pydantic model, ensuring that the system remains robust, self-documenting, and easy to validate. This disciplined approach is crucial for managing the increased complexity of this cycle. By defining the data structures for MD, MC, and FPS settings upfront, we establish clear contracts for what the new `MDMCExplorer` and `FPSSampler` components will consume. This prevents ambiguity and reduces the likelihood of integration errors. The use of `Literal` types in the Pydantic models (e.g., `Literal['random', 'fps']`) acts as a powerful form of static validation, ensuring that only valid methods or engines can be specified by the user, which is far more robust than relying on string comparisons deep within the application logic. This extension of the schema-first principle ensures that the application's configuration remains a single, coherent, and centrally-validated object, even as the tool's capabilities grow significantly.

**Pydantic Schema Design (`config.py` extension):**
The `FullConfig` model will be expanded to become the complete representation of any possible workflow. The following code is the precise blueprint for the additions to be made in `src/mlip_autopipec/config.py`.

```python
# config.py (additions and modifications)
from pydantic import BaseModel, Field
from typing import List, Dict, Literal, Optional

# ... existing models from Cycle 1 ...

class MDConfig(BaseModel):
    ensemble: Literal['npt', 'nvt', 'auto'] = Field('auto', description="Thermodynamic ensemble. 'auto' detects if the system is a slab.")
    temperature_K: float = Field(..., gt=0, description="MD simulation temperature in Kelvin.")
    pressure_bar: float = Field(1.0, description="MD simulation pressure in bar (for NPT).")
    timestep_fs: float = Field(1.0, gt=0, description="MD timestep in femtoseconds.")
    steps: int = Field(1000, gt=0, description="Total number of MD steps.")

class MCConfig(BaseModel):
    enabled: bool = Field(False, description="Enable hybrid MD/MC.")
    swap_frequency: int = Field(100, gt=0, description="Attempt an atom swap every N steps.")

class ExplorationConfig(BaseModel):
    engine: Literal['md', 'hybrid_md_mc'] = 'md'
    md: MDConfig
    mc: MCConfig
    mlip_model_path: str = Field(..., description="Path to the trained MLIP model file.")

class SamplingConfig(BaseModel):
    method: Literal['random', 'fps'] = 'random'
    num_samples: int = Field(..., gt=0, description="The final number of structures to select.")
    # Parameters for dscribe.SOAP
    fps_soap_config: Optional[Dict] = Field(None, description="Configuration for SOAP descriptors for FPS.")

# Modify the LabelingConfig from Cycle 1
class LabelingConfig(BaseModel):
    engine: Literal['mock', 'dft'] = 'dft'
    dft_code: Literal['quantum_espresso'] = 'quantum_espresso'
    dft_executable: str = "pw.x"
    # Further DFT parameters like pseudopotentials, k-points, etc. would go here

class FullConfig(BaseModel):
    """ The updated root model for the entire configuration file for Cycle 2. """
    system: SystemConfig
    generation: GenerationConfig
    # NEW sections below, replacing the simple labeling config
    exploration: ExplorationConfig
    sampling: SamplingConfig
    labeling: LabelingConfig
    storage: StorageConfig
```
This expanded schema provides the necessary fine-grained control over the new simulation and sampling capabilities. The architecture ensures that these validated and strongly-typed configuration objects are passed down to the specific components that need them, maintaining a clean and decoupled flow of information. For example, the `WorkflowOrchestrator` will extract the `ExplorationConfig` object and pass it, and only it, to the `MDMCExplorer` instance. This prevents components from having access to configuration settings that are not relevant to their function.

## 4. Implementation Approach

The implementation of Cycle 2 is a sequence of upgrades and additions, transforming the Cycle 1 skeleton into a feature-complete application. The approach is to build and test each new scientific component independently before integrating it into the main workflow.

1.  **Dependencies:** The first step is to update the `pyproject.toml` file to include the new dependencies required for this cycle's advanced features. This will include `dscribe` for calculating the SOAP descriptors needed for FPS, `streamlit` for the Web UI, and potentially `pymatgen` if it is required for the more advanced structure generators.

2.  **`DFTLabelingEngine` (`labeling.py`):** This is a top-priority task.
    *   A new class, `DFTLabelingEngine`, will be implemented, inheriting from a `BaseLabelingEngine`.
    *   Its `run` method will be the core of its functionality. Given an `ase.Atoms` object, it must perform several steps in a temporary directory:
        a. **Input File Generation:** It will use a template or a library like `ase.io.espresso.write_espresso_in` to generate the text input file for the specified DFT code (e.g., Quantum Espresso). This file will contain the atomic positions, cell parameters, and calculation settings derived from the `LabelingConfig`.
        b. **Subprocess Execution:** It will use `subprocess.run([self.executable, '-in', 'input.pwi'], shell=False, capture_output=True, text=True, check=True)`. The use of `shell=False` is critical for security, and `check=True` will automatically raise an exception if the DFT code returns a non-zero exit code, simplifying error handling.
        c. **Output Parsing:** Upon successful execution, it will read the captured `stdout` from the completed process. It will use a series of robust regular expressions to find and parse the final total energy, the block of atomic forces, and the final stress tensor from this text output.
    *   It will construct and return a `DFTResult` Pydantic model with the parsed data. Extensive logging will be added to this component to aid in debugging failed DFT runs.

3.  **Exploration Engine (`exploration.py`):**
    *   A new module, `core/exploration.py`, will be created.
    *   The `MDMCExplorer` class will be implemented. Its main method, `run_exploration`, will accept a list of `ase.Atoms` objects and the `ExplorationConfig`.
    *   Inside a loop over the input structures, it will:
        a. **Load Calculator:** Dynamically load the MLIP model from the path specified in the config and attach it as the ASE calculator.
        b. **Determine Ensemble:** Implement the `auto_ensemble` logic by checking `atoms.pbc.all()` to distinguish between bulk and surface systems and select the appropriate ASE dynamics object (e.g., `NPT` or `Langevin`).
        c. **Run Dynamics:** Attach a `Trajectory` logger to the dynamics object to save the frames. If hybrid mode is enabled, a custom function that performs a Monte Carlo atom swap will be attached as an observer (`dyn.attach(mc_swap_function, interval=mc_config.swap_frequency)`).
        d. The method will return a list of paths to the generated trajectory files.

4.  **Samplers (`sampling.py`):**
    *   The `BaseSampler` abstract class will be defined.
    *   `RandomSampler` will be implemented as a simple baseline.
    *   `FPSSampler` will be the main implementation. Its `sample` method will take a list of trajectory files. It will iterate through all frames in all trajectories, computing the SOAP descriptors for each one using the `dscribe` library and the settings from the `fps_soap_config`. It will then use an iterative FPS algorithm (either a library implementation or a custom one) to select the `num_samples` structures.

5.  **`WorkflowOrchestrator` Integration (`workflow.py`):**
    *   The main `run_workflow` method will be substantially updated to reflect the full pipeline sequence.
    *   The sequence of operations will now be:
        1. Generation -> `List[ase.Atoms]`
        2. Exploration -> `List[TrajectoryFiles]`
        3. Sampling -> `List[ase.Atoms]` (the final, curated set)
        4. Labeling -> Loop over the sampled atoms and get `List[DFTResult]`
        5. Storage -> Save the atoms and their corresponding results to the DB.
    *   The orchestrator will now instantiate and call the `MDMCExplorer` and the chosen `Sampler` from the new factory, passing the correct configuration objects and managing the data flow between these new stages.

6.  **Web UI (`webui/app.py`):**
    *   A new `cli.py` command, `webui`, will be added to launch the Streamlit application using `subprocess.run(['streamlit', 'run', 'app.py'])`.
    *   The `webui/app.py` file will use Streamlit's widgets (`st.slider`, `st.text_input`, `st.selectbox`) to create a form that mirrors the structure of the `FullConfig` Pydantic model.
    *   A "Run Pipeline" button (`st.button`) will, when clicked, gather the state from all the UI widgets, construct a `FullConfig` object, instantiate the main `WorkflowOrchestrator`, and call its `run_workflow` method. The standard output will be redirected to the Streamlit interface to provide real-time feedback to the user.

## 5. Test Strategy

The testing strategy for Cycle 2 is designed to rigorously validate the correctness of the new, complex scientific algorithms and to ensure their seamless and robust integration into the existing application framework. The focus shifts from architectural validation to algorithmic and integration correctness.

**Unit Testing Approach (Min 300 words):**
*   **`test_labeling.py`:** Unit testing the `DFTLabelingEngine` is the highest priority. Given its interaction with external processes, mocking is essential for creating fast, reliable, and portable tests. We will use `unittest.mock.patch` to completely replace the built-in `subprocess.run` function. Our mock will be a callable object that we can configure on a per-test basis. For the "happy path" test, we will create a string variable containing the text of a real, successful Quantum Espresso output file. We will configure our mock to return a `CompletedProcess` object where `returncode=0` and `stdout` is our success string. The test will then call the `labeling_engine.run` method and assert that the returned `DFTResult` object contains the exact energy, forces, and stress values that we know are in the sample output string. We will have another test for failure modes. We will create a sample output string from a QE run that failed to converge. The test will assert that our parsing logic correctly identifies this failure and raises a custom `DFTRuntimeError`. A third test will configure the mock to return `returncode=1`, and we will assert that this also raises our custom exception. This suite of tests provides high confidence in our engine's robustness without ever needing to install or run a real DFT code.

*   **`test_exploration.py`:** While testing the stochastic nature of a full MD run in a unit test is impractical, we can and must test the deterministic decision-making logic within the `MDMCExplorer`. The most important of these is the `auto_ensemble` feature. The test suite will contain a function `test_auto_ensemble_logic`. Inside this function, we will create two `ase.Atoms` objects. The first, `bulk_si`, will be a standard silicon crystal with `pbc=True` (or `pbc=[True, True, True]`). The second, `si_slab`, will be a silicon surface slab with `pbc=[True, True, False]`. We will then pass each of these objects to the internal helper function responsible for choosing the ensemble and assert that it returns the string `'npt'` for `bulk_si` and `'nvt'` for `si_slab`. This directly validates the core logic of this important feature.

*   **`test_sampling.py`:** Testing the `FPSSampler` requires a carefully constructed, non-trivial test case to ensure it is behaving more intelligently than a random sampler. The test function `test_fps_selects_diverse_structures` will create a synthetic trajectory. For example, it could be a 20-frame trajectory where the first 15 frames are created by applying only a very small, random rattle to the same initial structure (making them all very similar). The last 5 frames will be created by applying a large, distinct deformation to the initial structure, making them unique and different from the first 15. We will then configure the `FPSSampler` to select 5 structures from this 20-frame trajectory. A correct implementation of FPS must preferentially select the 5 unique, deformed structures. The test will assert that the indices of the 5 structures returned by the sampler are `[15, 16, 17, 18, 19]`. This provides a clear and unambiguous validation that the FPS algorithm is correctly identifying and prioritizing structural diversity.

**Integration Testing Approach (Min 300 words):**
*   **`test_full_pipeline_integration.py`:** A new, comprehensive integration test will be created in `tests/integration/` to verify the full, end-to-end advanced workflow. This test will serve as the final quality gate, ensuring all the new components from Cycle 2 are correctly integrated into the `WorkflowOrchestrator` and function together as expected. To ensure the test can run quickly and reliably in a CI/CD environment, several strategic choices will be made. The test will be configured to run on a very small atomic system (e.g., an 8-atom SiGe cell). The MD simulation will be set to run for a very short duration (e.g., 10-20 steps). Crucially, the `LabelingEngine` will still be mocked using `unittest.mock.patch`, even in this integration test. This is a deliberate decision to isolate the workflow logic from the unpredictable and time-consuming nature of real DFT calculations, ensuring the test remains fast and focused on the integration of the pipeline stages. The test will invoke the CLI with a YAML file configured to use the `hybrid_md_mc` explorer and the `fps` sampler (requesting, for instance, 5 samples). The assertions will then verify the entire chain of events:
    1.  The command exits with a `0` status code.
    2.  The final database file is created.
    3.  The number of rows in the database is exactly `5` (the number of samples requested by FPS), not `20` (the number of MD steps). This is the most critical assertion, as it proves the Exploration and Sampling stages were executed correctly and their output was passed to the storage stage.
    4.  The data within the database (e.g., the energy) matches the dummy data provided by our mocked labeling engine.

*   **Manual Web UI Testing:** The Web UI, being a highly interactive component, will be validated through a manual test plan. A section in the UAT document will provide a script for a human tester. The script will instruct the tester to launch the UI, use the various widgets to construct a configuration identical to the one used in the integration test, click the "Run" button, and observe the output. The final step of the manual test will be to use a tool (like `sqlite3` or a simple Python script provided in the UAT notebook) to inspect the database generated by the Web UI run and confirm that its contents are identical to the database generated by the CLI integration test. This ensures that both user interfaces are consistent front-ends for the same robust backend logic.
