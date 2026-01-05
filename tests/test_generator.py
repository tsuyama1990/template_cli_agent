"""Unit tests for the AlloyGenerator."""
import pytest
from collections import Counter
from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator

@pytest.fixture
def simple_alloy_config() -> SystemConfig:
    """Provides a simple SystemConfig for a binary alloy."""
    return SystemConfig(
        elements=['Fe', 'Pt'],
        composition={'Fe': 0.5, 'Pt': 0.5},
        lattice='fcc',
        num_structures=5
    )

def test_alloy_generator_correct_count(simple_alloy_config: SystemConfig):
    """
    Test that the AlloyGenerator produces the correct number of structures.
    """
    generator = AlloyGenerator(config=simple_alloy_config)
    structures = generator.generate()
    assert len(structures) == simple_alloy_config.num_structures

def test_alloy_generator_correct_composition(simple_alloy_config: SystemConfig):
    """
    Test that each generated structure has the correct composition.
    """
    generator = AlloyGenerator(config=simple_alloy_config)
    structures = generator.generate()

    assert len(structures) > 0, "Generator should produce at least one structure."

    for atoms in structures:
        # A 3x3x3 supercell of a conventional fcc cell (4 atoms) has 108 atoms.
        expected_total_atoms = 108
        assert len(atoms) == expected_total_atoms

        expected_fe_count = 54 # 108 * 0.5
        expected_pt_count = 54 # 108 * 0.5

        counts = Counter(atoms.get_chemical_symbols())

        assert counts['Fe'] == expected_fe_count
        assert counts['Pt'] == expected_pt_count


def test_alloy_generator_physical_validity(simple_alloy_config: SystemConfig):
    """
    Test that all generated structures are physically plausible.
    """
    generator = AlloyGenerator(config=simple_alloy_config)
    structures = generator.generate()

    assert len(structures) > 0

    for atoms in structures:
        # Use the generator's own validation method to check
        assert generator._is_physically_valid(atoms, min_dist=1.0)
