# CYCLE 05: SPEC.md - Advanced Features and Usability

## 1. Summary

Cycle 05 is the final phase of the initial development plan, focusing on refining the MLIP-AutoPipe system by expanding its capabilities to handle more complex material types and advanced physical phenomena. This cycle is also dedicated to improving overall usability and robustness, transforming the powerful core engine into a more polished and user-friendly tool. The primary objectives are twofold: first, to broaden the scientific scope of the application, and second, to enhance the user experience with better feedback and documentation.

On the scientific front, this cycle will enhance the `StructureGenerator` to support ionic materials (using Ab Initio Random Structure Searching - AIRSS) and covalent/amorphous systems (using a Melt-Quench protocol). This will make the pipeline applicable to a much wider range of materials beyond the alloys and molecules covered in Cycle 02. The `LabelingEngine` will also be upgraded to handle more complex DFT calculations, specifically for magnetic materials, by introducing automated initialization of spin states and logic for handling Spin-Orbit Coupling (SOC). This is a critical feature for accurately modeling many technologically important materials. On the usability front, this cycle will focus on improving the command-line interface with clearer progress indicators and providing better user documentation and tutorials to facilitate adoption.

## 2. System Architecture

The architectural changes in this cycle are primarily enhancements to existing modules rather than the introduction of new ones. The `StructureGenerator` and `LabelingEngine` will see the most significant modifications.

**File Structure (Cycle 05 Focus):**

```
.
├── src/
│   └── mlip_autopipec/
│       ├── **cli.py**                # Enhanced with better progress indicators
│       ├── configs/
│       │   └── **models.py**         # Add config options for new features
│       ├── modules/
│       │   ├── **a_structure_generator.py** # Add Ionic (AIRSS) & Covalent (Melt-Quench)
│       │   └── **c_labeling_engine.py**     # Add magnetism and SOC handling
│       └── ...
└── docs/                             # NEW directory for user documentation
    └── **tutorial.md**
```

The files marked in **bold** are the primary areas of work. The `a_structure_generator.py` and `c_labeling_engine.py` will be extended with new functionality. The `cli.py` will be polished for better user feedback. A new `docs` directory will be created to house user-facing documentation, starting with a comprehensive tutorial.

## 3. Design Architecture

The design for this cycle focuses on adding conditional logic to existing modules to handle the new, advanced cases.

**Pydantic Schema (`configs/models.py`):**

The configuration will be extended to control the new features.

*   `DFTComputeConfig` will be updated with:
    *   `magnetism: Optional[str] = None` (e.g., "ferromagnetic" or "random_afm")
    *   `soc_enabled: bool = False`

*   `StructureGeneratorConfig` (a new model, or added to `MainConfig`):
    *   `melt_quench_temps: Optional[List[float]] = None` (e.g., `[3000, 300]`)
    *   `airss_max_steps: Optional[int] = None`

**Class and Module Design:**

*   **`StructureGenerator` (`modules/a_structure_generator.py`):**
    *   The `run()` method will be extended with new conditional blocks based on the `structure_type` determined by the `ConfigExpander`.
    *   **If 'ionic'**: Implement a simplified AIRSS protocol. This involves generating random atomic positions within a cell, subject to constraints on inter-atomic distances based on ionic radii, followed by a quick relaxation with the surrogate model to generate sensible candidate structures.
    *   **If 'covalent'**: Implement a Melt-Quench protocol. This will use the surrogate model (MACE) to run an MD simulation that first heats the system to a high temperature (melting it) and then rapidly cools it down (quenching) to generate amorphous or defected structures.

*   **`LabelingEngine` (`modules/c_labeling_engine.py`):**
    *   The logic for generating Quantum Espresso input files will be significantly enhanced.
    *   **Magnetism Handling:**
        1.  If the `magnetism` flag is set in the config, the engine will add the necessary QE parameters to the input file (e.g., `nspin=2`, and `starting_magnetization`).
        2.  It will implement logic to automatically initialize the magnetic moments. For "ferromagnetic", all moments will be aligned. For more advanced cases, it might involve randomizing initial spins.
        3.  The error recovery logic will be updated to handle convergence issues common in magnetic calculations (e.g., by attempting different spin initializations).
    *   **SOC Handling:**
        1.  If `soc_enabled` is true, it will add the parameters for a non-collinear calculation (`noncolin = .true.`, `lspinorb = .true.`) to the QE input, which is a more computationally demanding calculation required for SOC.

*   **`cli.py`:**
    *   The CLI will be improved by adding progress bars for long-running tasks, such as the labeling of many structures or the training process. Libraries like `rich` or `tqdm` will be used to provide a better user experience.

## 4. Implementation Approach

1.  **Update Configuration:** Add the new Pydantic model fields for magnetism, SOC, and the new structure generation methods. Update the `ConfigExpander` to recognize ionic and covalent material types and set appropriate default configurations.
2.  **Enhance StructureGenerator:**
    *   Implement the AIRSS logic. This can be done using ASE's built-in tools for generating random structures combined with a surrogate-model relaxation step.
    *   Implement the Melt-Quench protocol. This will involve writing a new MD simulation loop within the generator that uses the surrogate model to control the temperature according to the specified profile.
3.  **Enhance LabelingEngine:**
    *   Modify the QE input writer function to correctly add the `SYSTEM` and `ELECTRONS` card parameters for magnetic and SOC calculations based on the configuration.
    *   Implement the logic for initializing magnetic moments on the `ase.Atoms` object before the input file is written.
4.  **Improve CLI:** Integrate the `rich` library into `cli.py`. Wrap the loops in the `WorkflowOrchestrator` (e.g., the loop over structures in the `LabelingEngine`) with `rich.progress` to display real-time progress to the user.
5.  **Write Documentation:** Create the `docs/` directory. Write a `tutorial.md` file that walks a new user through a complete, end-to-end example, explaining the concepts, the minimal `input.yaml`, how to run the tool, and how to interpret the results.

## 5. Test Strategy

Testing will focus on ensuring the new, complex logic is correctly triggered and executed.

**Unit Testing Approach:**
(Located in `tests/unit/`)
*   **`modules/a_structure_generator.py`:**
    *   Write a test for the AIRSS logic, asserting that it produces randomly generated but physically plausible structures (e.g., atoms are not too close).
    *   Write a test for the Melt-Quench logic. Mock the surrogate MD run and assert that the final structures are disordered (e.g., by checking their radial distribution function).
*   **`modules/c_labeling_engine.py`:**
    *   This is a critical test. Write a new test case for a magnetic material like Iron (Fe).
    *   Provide a config with `magnetism: "ferromagnetic"`.
    *   Run the engine's input generation logic.
    *   Assert that the generated Quantum Espresso input string contains the correct QE flags: `nspin = 2` and `starting_magnetization(1) = ...`.
    *   Create another test case for SOC, asserting that `noncolin = .true.` is present in the input string when `soc_enabled: True`.

**Integration Testing Approach:**
(Located in `tests/e2e/`)
*   **Magnetic Material End-to-End Test:**
    *   **Setup:** Create a minimal `input.yaml` for `elements: ["Fe"]`.
    *   **Execution:** Run the full pipeline using the CLI runner. Mock the `subprocess.run` call for QE. However, the mock will now inspect the input it receives.
    *   **Verification:**
        1.  The test will assert that the `ConfigExpander` correctly identified Fe as a magnetic material and set the magnetism flag.
        2.  The primary assertion will be within the `subprocess.run` mock: it will check that the QE input file it was asked to execute contains the correct magnetism parameters. This verifies that the configuration flowed correctly all the way through the system to the final external process call.
*   **CLI UI Test:**
    *   A simple test to run the CLI on a standard workflow and capture the `stdout`. The test will then assert that the output contains the patterns characteristic of the `rich` progress bar, confirming that the UI enhancements have been applied.
