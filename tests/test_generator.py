from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator():
    """Tests the AlloyGenerator for correct count and composition."""
    config_dict = {
        "elements": ["Fe", "Pt"],
        "composition": {"Fe": 0.5, "Pt": 0.5},
        "lattice": "fcc",
        "num_structures": 10,
    }
    system_config = SystemConfig(**config_dict)
    generator = AlloyGenerator(system_config)

    structures = generator.generate()

    # Check that the correct number of structures were generated
    assert len(structures) == system_config.num_structures

    for atoms in structures:
        # Check that the composition is correct
        composition = atoms.get_chemical_symbols()
        num_fe = composition.count("Fe")
        num_pt = composition.count("Pt")
        total_atoms = len(atoms)
        assert abs(num_fe / total_atoms - 0.5) < 0.1
        assert abs(num_pt / total_atoms - 0.5) < 0.1

        # Check for physical validity (no overlapping atoms)
        distances = atoms.get_all_distances(mic=True)
        min_distance = distances[distances > 0].min()
        assert min_distance > 1.0
