"""Unit tests for the structure generator component."""
import numpy as np
import pytest
from ase.atoms import Atoms

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.generators.base import BaseStructureGenerator


def test_alloy_generator_interface() -> None:
    """Test that AlloyGenerator correctly implements the base interface."""
    assert issubclass(AlloyGenerator, BaseStructureGenerator)


def test_alloy_generator_creation() -> None:
    """Test the creation of alloy structures."""
    config_dict = {
        "elements": ["Cu", "Au"],
        "composition": {"Cu": 0.5, "Au": 0.5},
        "lattice": "fcc",
        "num_structures": 5,
    }
    system_config = SystemConfig(**config_dict)
    generator = AlloyGenerator(system_config)
    structures = generator.generate()

    assert isinstance(structures, list)
    assert len(structures) == 5
    assert all(isinstance(s, Atoms) for s in structures)

    # Check a single structure for correctness
    atoms = structures[0]
    assert len(atoms) > 0  # Should have atoms
    symbols = atoms.get_chemical_symbols()
    assert "Cu" in symbols
    assert "Au" in symbols
    assert len(symbols) == atoms.get_number_of_atoms()
    # Check composition ratio (allowing for slight floating point inaccuracies)
    cu_count = symbols.count("Cu")
    au_count = symbols.count("Au")
    total_atoms = len(symbols)
    assert np.isclose(cu_count / total_atoms, 0.5, atol=0.1)
    assert np.isclose(au_count / total_atoms, 0.5, atol=0.1)


def test_physical_validity_check() -> None:
    """Verify that generated structures are physically plausible."""
    config_dict = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "bcc",
        "num_structures": 1,
    }
    system_config = SystemConfig(**config_dict)
    generator = AlloyGenerator(system_config)
    structures = generator.generate()
    atoms = structures[0]

    distances = atoms.get_all_distances(mic=True)
    min_dist = np.min(distances[np.nonzero(distances)])
    assert min_dist > 1.0, f"Atoms are overlapping with min_dist={min_dist}"
