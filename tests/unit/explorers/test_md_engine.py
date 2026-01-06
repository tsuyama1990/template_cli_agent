# tests/unit/explorers/test_md_engine.py

import pytest
from ase import Atoms
from ase.calculators.emt import EMT

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.explorers.md_engine import MDEngine


@pytest.fixture
def md_config() -> FullConfig:
    """Provides a valid config for an MD run."""
    config_dict = {
        "system": {
            "elements": ["H"],
            "composition": {"H": 1.0},
            "supercell_size": [1, 1, 1],
        },
        "exploration": {
            "temperature_k": 300,
            "pressure_gpa": 0,
            "timestep_fs": 1.0,
            "n_steps": 10,  # A short run for testing
        },
        "sampling": {"n_samples": 1},
    }
    return FullConfig.model_validate(config_dict)


def test_md_engine_runs_simulation_and_creates_trajectory(md_config: FullConfig) -> None:
    """
    Test that the MDEngine correctly runs a simulation and returns a trajectory.
    """
    # Use at least two atoms to avoid ZeroDivisionError in Langevin dynamics
    initial_structure = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    initial_structure.set_cell([10, 10, 10])  # type: ignore[no-untyped-call]

    # Use the EMT calculator as it's a simple, fast, and real calculator
    # that can handle dynamic atoms, unlike SinglePointCalculator.
    calculator = EMT()  # type: ignore[no-untyped-call]

    engine = MDEngine(config=md_config, calculator=calculator)

    # Run the exploration
    trajectory = engine.explore(initial_structure)

    # The number of frames depends on the logging frequency of the dynamics object.
    # For this test, we just care that it produced more than the initial frame.
    assert len(trajectory) > 1
    assert isinstance(trajectory[0], Atoms)

    # Check that the final structure has energy and forces attached
    final_atoms = trajectory[-1]
    assert "energy" in final_atoms.info
    assert "forces" in final_atoms.arrays
