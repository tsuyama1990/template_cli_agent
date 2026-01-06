# tests/unit/generators/test_base_generator.py
import pytest
from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.base import BaseStructureGenerator


class ConcreteGenerator(BaseStructureGenerator):
    """A concrete implementation of the abstract base class for testing."""

    def generate(self) -> list[Atoms]:
        # This method is not needed for testing the validation logic.
        return []


@pytest.fixture
def mock_config() -> FullConfig:
    """Provides a mock FullConfig object for testing."""
    config_dict = {
        "system": {
            "elements": ["H"],
            "composition": {"H": 1.0},
            "supercell_size": [1, 1, 1],
        },
        "exploration": {"temperature_k": 300, "pressure_gpa": 0, "timestep_fs": 1.0, "n_steps": 100},
        "sampling": {"n_samples": 1},
        "db_path": "test.db",
    }
    return FullConfig.model_validate(config_dict)


def test_validate_structure_accepts_valid_structure(mock_config: FullConfig) -> None:
    """Test that a structure with no overlaps is considered valid."""
    generator = ConcreteGenerator(config=mock_config)
    # Create two atoms far apart
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 2.0]])
    assert generator._validate_structure(atoms) is True


def test_validate_structure_rejects_overlapping_structure(mock_config: FullConfig) -> None:
    """Test that a structure with overlapping atoms is considered invalid."""
    generator = ConcreteGenerator(config=mock_config)
    # Create two atoms that are too close
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.5]])
    assert generator._validate_structure(atoms) is False
