# tests/mlip_autopipec/unit/test_generators.py
"""Unit tests for structure generators."""
from collections import Counter

from mlip_autopipec.common.pydantic_models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator_creates_correct_structure() -> None:
    """Test that the AlloyGenerator produces a structure with the correct properties."""
    # 1. Define the configuration for a simple binary alloy
    config = SystemConfig(
        elements=["Fe", "Pt"],
        composition={"Fe": 0.75, "Pt": 0.25},
        supercell_size=[2, 2, 2],
    )

    # 2. Initialize the generator and create the structure
    generator = AlloyGenerator(config)
    generated_structures = generator.generate()

    # 3. Perform assertions
    # Ensure it returns a list with a single Atoms object
    assert isinstance(generated_structures, list)
    assert len(generated_structures) == 1
    atoms = generated_structures[0]

    # A 2x2x2 supercell of an FCC primitive cell (4 atoms) should have 32 atoms
    expected_total_atoms = 4 * (2 * 2 * 2)
    assert len(atoms) == expected_total_atoms

    # Verify the composition
    expected_fe_count = int(round(expected_total_atoms * 0.75))
    expected_pt_count = int(round(expected_total_atoms * 0.25))

    actual_counts = Counter(atoms.get_chemical_symbols())

    assert actual_counts["Fe"] == expected_fe_count
    assert actual_counts["Pt"] == expected_pt_count

def test_alloy_generator_handles_composition_rounding() -> None:
    """Test the edge case where composition rounding can affect atom counts."""
    config = SystemConfig(
        elements=["Cu", "Au"],
        composition={"Cu": 0.333, "Au": 0.667},
        supercell_size=[2, 2, 2],
    )

    generator = AlloyGenerator(config)
    generated_structures = generator.generate()
    atoms = generated_structures[0]

    expected_total_atoms = 32
    assert len(atoms) == expected_total_atoms

    actual_counts = Counter(atoms.get_chemical_symbols())
    assert actual_counts["Cu"] + actual_counts["Au"] == expected_total_atoms
