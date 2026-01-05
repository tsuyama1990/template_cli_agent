"""Unit tests for the structure generators."""

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator():
    """Test the AlloyGenerator for correct count, composition, and validity."""
    config_data = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": 5,
    }
    system_config = SystemConfig(**config_data)
    generator = AlloyGenerator(system_config)
    structures = generator.generate()

    assert len(structures) == 5

    for atoms in structures:
        assert len(atoms) >= 64  # Ensure supercell is large enough
        symbols = atoms.get_chemical_symbols()
        num_fe = symbols.count("Fe")
        num_pt = symbols.count("Pt")
        total_atoms = len(symbols)

        # Check if composition is approximately correct
        assert abs(num_fe / total_atoms - 0.5) < 0.1
        assert abs(num_pt / total_atoms - 0.5) < 0.1

        # Check physical validity
        distances = atoms.get_all_distances(mic=True)
        min_dist = min(d for d in distances.flatten() if d > 0)
        assert min_dist > 1.0
