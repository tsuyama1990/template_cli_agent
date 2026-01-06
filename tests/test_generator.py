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
        num_fe = sum(atoms.numbers == 26)
        num_pt = sum(atoms.numbers == 78)
        total_atoms = len(atoms)
        assert abs(num_fe / total_atoms - 0.5) < 0.1
        assert abs(num_pt / total_atoms - 0.5) < 0.1


def test_physical_validity():
    """Test that the generated structures are physically valid (no overlapping atoms)."""
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
        min_distance = distances[distances > 0].min()
        assert min_distance > 1.0
