"""Unit tests for the AlloyGenerator."""

import pytest
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.common.pydantic_models import SystemConfig
from ase import Atoms

@pytest.fixture
def alloy_config() -> SystemConfig:
    """Provides a default SystemConfig for an alloy."""
    return SystemConfig(
        elements=["Fe", "Pt"],
        composition={"Fe": 0.5, "Pt": 0.5},
        supercell_size=[3, 3, 3]
    )

def test_alloy_generator_initialization(alloy_config: SystemConfig) -> None:
    """Tests that the AlloyGenerator can be initialized."""
    generator = AlloyGenerator(config=alloy_config)
    assert generator is not None

def test_alloy_generator_returns_list_of_atoms(alloy_config: SystemConfig) -> None:
    """
    Tests that the generate method returns a list of ASE Atoms objects.
    """
    generator = AlloyGenerator(config=alloy_config)
    structures = generator.generate()
    assert isinstance(structures, list)
    assert len(structures) > 0
    assert isinstance(structures[0], Atoms)
