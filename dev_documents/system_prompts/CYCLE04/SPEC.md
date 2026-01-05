# Cycle 4 Specification: Robustness, Advanced Generators, and Physics-Based Intelligence

## 1. Summary

Cycle 4 transitions the MLIP-AutoPipe from a powerful scientific tool into a robust, production-ready engineering application. The focus of this cycle is on stability, automation, and the incorporation of deeper physical knowledge into the pipeline. While previous cycles built the core functionality, this cycle makes it "smart" and resilient, capable of handling more complex systems and running reliably for extended periods without supervision.

The key deliverables for this cycle are:
1.  **Enhanced Simulation Robustness:** High-temperature simulations, which are crucial for generating diverse data, are prone to failure (e.g., atoms getting too close, leading to numerical instability). This cycle will implement two critical features to mitigate this:
    *   **ZBL Potential Mixing:** The MLIP calculator will be wrapped with a mechanism that blends it with a classical ZBL potential at very short interatomic distances. This provides a strong, physically correct repulsive force that prevents atoms from fusing, dramatically improving simulation stability.
    *   **Robust Error Handling:** The `ExplorationEngine`'s worker processes will be wrapped in a robust `try...except` block. If a simulation fails catastrophically (a "Coulomb explosion"), the error will be caught, the failed structure will be logged and saved for debugging, and the worker will exit gracefully without crashing the entire pipeline.
2.  **Advanced `IonicGenerator`:** A new, sophisticated structure generator for ionic materials (e.g., oxides) will be implemented. This generator will go beyond simple random placement and will enforce charge neutrality, a fundamental physical constraint in ionic systems.
3.  **Physics-Based Intelligence (`Auto Ensemble Switching`):** The pipeline will be imbued with the ability to automatically inspect a system and choose the correct simulation parameters. A `detect_vacuum` function will be implemented to determine if a structure is a bulk material or a slab/surface. Based on this, the `ExplorationEngine` will automatically switch between using an NPT (constant pressure, for relaxing bulk crystals) and an NVT (constant volume, for preserving surfaces) thermodynamic ensemble. This removes a complex choice from the user and prevents common simulation artifacts.

By the end of Cycle 4, MLIP-AutoPipe will be significantly more automated and reliable. It will be able to handle a wider class of materials (ionic crystals), run high-temperature simulations with a much lower failure rate, and make intelligent, physics-based decisions about simulation parameters on behalf of the user, truly advancing the goal of "removing the human expert from the loop."

## 2. System Architecture

This cycle introduces a new `IonicGenerator`, enhances the core MD explorer with major robustness features, and adds new utilities for physics-based analysis. Files to be created or modified are in **bold**.

```
src/mlip_autopipec/
├── __init__.py
├── **config.py**               # Pydantic models updated for new features
├── database/
│   └── ase_db_wrapper.py       # (No changes)
├── domain/
│   └── ...                     # (No changes)
├── engines/
│   ├── generation_engine.py    # Updated to select IonicGenerator
│   └── **exploration_engine.py** # Updated to handle failures and select ensemble
├── explorers/
│   └── **md_mc_explorer.py**     # Enhanced with ZBL mixing
├── generators/
│   ├── __init__.py
│   ├── base_generator.py
│   ├── alloy_generator.py
│   └── **ionic_generator.py**    # **New implementation**
├── samplers/
│   └── ...                     # (No changes)
├── utils/
│   ├── __init__.py
│   └── **physics.py**            # **New file** for vacuum detection
└── workflow_orchestrator.py      # (No changes)
```

**File Blueprints:**

*   **`config.py`**:
    *   A new `IonicGeneratorConfig` Pydantic model will be created, with fields for `cation_elements`, `anion_elements`, their respective `oxidation_states`, etc.
    *   `MainConfig` will be updated to use a discriminated union for the `generator` field, allowing the user to select `"alloy"` or `"ionic"` mode from the YAML file.
    *   A `ZBLConfig` model will be added to `ExplorationConfig` to allow the user to enable and tune the ZBL potential mixing.
    *   A field for `ensemble_mode: "auto" | "nvt" | "npt"` will be added to `ExplorationConfig`.

*   **`generators/ionic_generator.py`**:
    *   A new `IonicGenerator` class will be implemented, inheriting from `BaseGenerator`.
    *   Its `generate` method will be significantly more complex than the alloy generator. It will need to place cations and anions in a way that the total charge of the unit cell sums to zero. It may use heuristics or simple geometric packing algorithms to achieve a reasonable starting structure. It will still be subject to the `_check_overlap` validation from the base class.

*   **`utils/physics.py`**:
    *   A new file for housing physics-related utility functions.
    *   The primary function will be `detect_vacuum(atoms, probe_resolution=2.0)`. This function will implement a grid-based algorithm: it will create a 3D grid of points within the simulation cell, calculate the distance from each grid point to the nearest atom, and if any grid point is further away than a certain threshold, it will classify the system as having a vacuum (i.e., being a slab or surface).

*   **`explorers/md_mc_explorer.py`**:
    *   The `HybridMDMC` class will be modified. The calculator setup logic will be enhanced to wrap the MLIP calculator (e.g., MACE) inside ASE's `mixing.LinearCombination` calculator, combining it with an `ase.calculators.zbl.ZBL` calculator. The ZBL potential will be configured to only be active at very short ranges.
    *   The logic for selecting the ASE dynamics object will be updated. Instead of being hardcoded to `Langevin` (NVT), it will now conditionally choose between `Langevin` (NVT) and `NPT` based on the ensemble mode passed in its configuration.

*   **`engines/exploration_engine.py`**:
    *   The `run()` method will be updated to first call the `detect_vacuum` utility for each structure. This will determine the `ensemble_mode` to be used for that specific simulation.
    *   The `run()` method's `executor.submit` call will be wrapped in a `try...except` block. If the worker process raises an exception (e.g., a `PhysicsViolationError` from the simulation), the main engine process will catch it.
    *   Upon catching an exception, the engine will log the failure, save the problematic input structure to a "quarantine" directory for later inspection, and then continue processing the rest of the structures, ensuring the entire pipeline does not halt due to one failed simulation.

## 3. Design Architecture

The design for Cycle 4 focuses on delegating decisions to autonomous components and building resilient, fault-tolerant systems.

**Pydantic Models (`config.py`):**

```python
from pydantic import BaseModel, Field
from typing import Literal, Union

# ... existing configs ...

class IonicGeneratorConfig(BaseModel):
    # ... fields for ions, oxidation states, etc. ...
    mode: Literal["ionic"] = "ionic"

class AlloyGeneratorConfig(BaseModel):
    # ... existing fields ...
    mode: Literal["alloy"] = "alloy"

class ExplorationConfig(BaseModel):
    # ... existing fields ...
    ensemble: Literal["auto", "nvt", "npt"] = "auto"
    enable_zbl_repulsion: bool = True

class MainConfig(BaseModel):
    generator: Union[AlloyGeneratorConfig, IonicGeneratorConfig] = Field(..., discriminator='mode')
    exploration: Optional[ExplorationConfig] = None
    # ... rest of config ...
```

**Design Principles:**

*   **Fault Tolerance:** The previous `ExplorationEngine` was brittle; a single simulation crash would halt all progress. The new design follows a robust "supervisor" pattern. The engine acts as a supervisor that dispatches jobs to disposable workers. If a worker fails, the supervisor logs the failure, isolates the problematic input, and continues, ensuring the overall process completes. This is a critical feature for long-running, autonomous pipelines.
*   **Delegation of Responsibility:** The choice of thermodynamic ensemble is a complex one that depends on the physical system. Instead of forcing the user to know this, the new design delegates this responsibility to the `physics.py` module. The `ExplorationEngine` *asks* the `detect_vacuum` function what kind of system it is dealing with and then acts on that information. This is a key principle of building "smart" or autonomous systems: encapsulate expert knowledge into specific, queryable components.
*   **Decorator/Wrapper Pattern:** The ZBL potential mixing is a perfect example of the Decorator or Wrapper pattern. The core `Mace` calculator object remains unchanged. We "decorate" it with additional functionality by wrapping it inside another ASE calculator (`LinearCombination`). This wrapper intercepts calls to calculate forces and energies, combines the results from MACE and ZBL, and returns the combined value. The rest of the simulation code is completely unaware of this change; it just interacts with the wrapper as if it were a normal calculator. This is a very clean and extensible way to add functionality.

## 4. Implementation Approach

1.  **Physics Utility (`physics.py`):**
    *   Create the `utils/physics.py` file.
    *   Implement `detect_vacuum`. This will involve getting the cell dimensions, creating a `numpy.meshgrid` of points, getting atom positions, and using efficient NumPy/SciPy methods to calculate the minimum distance from each grid point to any atom.

2.  **Ionic Generator (`ionic_generator.py`):**
    *   Create the `IonicGenerator` class.
    *   Implement the logic to generate a charge-neutral structure. A simple approach is to determine the simplest integer ratio of cations to anions that results in a neutral formula unit (e.g., MgCl2 -> 1 Mg, 2 Cl) and then build a supercell that contains an integer number of these formula units. The atoms can then be placed on a lattice and randomly perturbed.

3.  **ZBL Integration (`md_mc_explorer.py`):**
    *   Modify the calculator setup section in the `HybridMDMC` class.
    *   Add an `if config.enable_zbl_repulsion:` block.
    *   Inside the block, create the `MACE` calculator, create a `ZBL` calculator, and then create the `LinearCombination` calculator that combines them. The combined calculator is then attached to the atoms object.

4.  **Auto-Ensemble and Error Handling (`exploration_engine.py`):**
    *   In the `ExplorationEngine.run` method, before the main loop that submits jobs, add a loop that iterates through the structures and calls `detect_vacuum` for each one, storing the result (e.g., in a dictionary mapping structure ID to "nvt" or "npt").
    *   Modify the call to `executor.submit` to pass the determined ensemble mode to the `HybridMDMC` worker.
    *   Wrap the `future.result()` call (or the whole submission loop) inside a `try...except Exception as e:` block.
    *   Inside the `except` block, add logging to record the error and the ID of the failing structure. Add code to save the input `atoms` object that caused the crash to a file like `quarantine/failed_structure_{id}.xyz`.

## 5. Test Strategy

Testing this cycle's features is focused on verifying the correctness of the new physical models and the resilience of the pipeline.

**Unit Testing Approach (Min 300 words):**

*   **`detect_vacuum`:** This function is critical and must be tested thoroughly.
    1.  Test 1 (Bulk): Create a standard `ase.build.bulk("Cu")` structure. Assert that `detect_vacuum` returns `False`.
    2.  Test 2 (Slab): Create a slab structure using `ase.build.surface`. This structure will have a vacuum layer. Assert that `detect_vacuum` returns `True`.
    3.  Test 3 (Molecule): Create a single molecule in a large box. Assert that `detect_vacuum` returns `True`. This ensures it works for non-periodic systems as well.
*   **`IonicGenerator`:**
    1.  Instantiate the generator for a known system like MgO (Mg: +2, O: -2).
    2.  Generate a few structures.
    3.  For each structure, get the list of chemical symbols, map them to their oxidation states, and assert that the sum of all oxidation states in the cell is exactly zero.
*   **ZBL Integration:** We can't easily test the combined potential directly, but we can test the setup. We will create a `HybridMDMC` instance and inspect its `atoms` object. We will assert that `atoms.calc` is an instance of `ase.calculators.mixing.LinearCombination`, and that this combination calculator contains both a `MACE` instance and a `ZBL` instance. This verifies the wrapper pattern was implemented correctly.
*   **`ExplorationEngine` Error Handling:** We will test this using mocks. We will create a mock worker function that is programmed to raise an exception. We will submit this mock function to the engine. We will then assert that the engine does *not* crash, that the error was logged, and that a mock "save to quarantine" function was called.

**Integration Testing Approach (Min 300 words):**

The main integration test will be a "stress test" designed to trigger the new robustness features.
*   **Test Setup:**
    1.  The test will configure the `IonicGenerator` to create a set of initial structures for a system known to be unstable at high temperatures.
    2.  The `ExplorationConfig` will be set to a very high temperature (e.g., 2000 K).
    3.  Crucially, the test will be run **twice**.
        *   **Run 1:** `enable_zbl_repulsion: false`.
        *   **Run 2:** `enable_zbl_repulsion: true`.
*   **Execution:** The full pipeline will be executed for both runs.
*   **Assertions:**
    1.  **Run 1 (ZBL off):** We expect many simulations to fail. The test will assert that the pipeline *completes* (i.e., does not crash), but that the logs show multiple caught exceptions. It will also assert that a `quarantine/` directory has been created and contains failed structures.
    2.  **Run 2 (ZBL on):** We expect the simulations to be much more stable. The test will assert that the pipeline completes and that the number of caught exceptions is significantly lower (ideally zero) compared to Run 1. This provides a direct, comparative validation of the ZBL mixing's effectiveness.
*   **Auto-Ensemble Test:** A separate integration test will generate both a bulk and a slab structure. It will run the exploration with `ensemble: "auto"`. The test will need to inspect the logs produced by the pipeline (which should be modified to state which ensemble is being used) to confirm that the bulk structure was run with NPT and the slab structure was run with NVT.
