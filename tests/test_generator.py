"""Unit tests for the structure generators."""
import numpy as np

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator_count_and_composition():
    """Test that the AlloyGenerator produces the correct number and composition of structures."""
    config = SystemConfig(
        elements=["Fe", "Pt"],
        composition={"Fe": 0.5, "Pt": 0.5},
        lattice="fcc",
        num_structures=5,
    )
    generator = AlloyGenerator(config)
    structures = generator.generate()

    assert len(structures) == 5

    for atoms in structures:
        symbols = atoms.get_chemical_symbols()
        num_atoms = len(symbols)
        expected_fe = int(round(0.5 * num_atoms))
        expected_pt = num_atoms - expected_fe
        assert symbols.count("Fe") == expected_fe
        assert symbols.count("Pt") == expected_pt


def test_alloy_generator_physical_validity():
    """Test that the generated structures are physically valid."""
    config = SystemConfig(
        elements=["Fe", "Pt"],
        composition={"Fe": 0.5, "Pt": 0.5},
        lattice="fcc",
        num_structures=1,
    )
    generator = AlloyGenerator(config)
    structures = generator.generate()

    for atoms in structures:
        distances = atoms.get_all_distances(mic=True)
        min_dist = np.min(distances[np.nonzero(distances)])
        assert min_dist > 1.0
