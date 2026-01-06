"""Unit tests for the structure generators."""

import numpy as np

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator_output() -> None:
    # Corresponds to SPEC.md, Section 2: Component Blueprint (`generators/alloy.py`)
    # Verifies that the AlloyGenerator produces the correct number of structures,
    # with the correct composition, and that they are physically valid.
    # This also supports UAT-C1-005, which checks for non-overlapping atoms.
    """Test the AlloyGenerator to ensure it produces valid structures.

    This test checks for correct composition and count.
    """
    # 1. Define a configuration for a simple binary alloy
    config_dict = {
        "elements": ["Cu", "Au"],
        "composition": {"Cu": 0.75, "Au": 0.25},
        "lattice": "fcc",
        "num_structures": 5,
    }
    system_config = SystemConfig(**config_dict)

    # 2. Instantiate the generator
    generator = AlloyGenerator(config=system_config)

    # 3. Generate the structures
    structures = generator.generate()

    # 4. Assert the correctness of the output
    # Check if the correct number of structures was generated
    assert len(structures) == system_config.num_structures

    for atoms in structures:
        # Check the total number of atoms in the supercell (4x4x4 fcc = 256 atoms)
        assert len(atoms) == 256

        # Check the composition
        symbols = atoms.get_chemical_symbols()
        num_cu = symbols.count("Cu")
        num_au = symbols.count("Au")

        expected_cu = int(0.75 * 256)
        expected_au = int(0.25 * 256)

        # Allow for a small tolerance due to rounding
        assert abs(num_cu - expected_cu) <= 1
        assert abs(num_au - expected_au) <= 1
        assert num_cu + num_au == 256

        # Check for physical validity (no overlapping atoms)
        distances = atoms.get_all_distances(mic=True)
        min_dist = np.min(distances[np.nonzero(distances)])
        assert min_dist > 1.0, "Generated structure has overlapping atoms."
