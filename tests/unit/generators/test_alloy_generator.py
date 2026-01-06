# tests/unit/generators/test_alloy_generator.py
from unittest.mock import patch

import pytest

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


@pytest.fixture
def alloy_config() -> FullConfig:
    """Provides a valid config for a binary alloy."""
    config_dict = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {
            "temperature_k": 300,
            "pressure_gpa": 0,
            "timestep_fs": 1.0,
            "n_steps": 100,
        },
        "sampling": {"n_samples": 1},
    }
    return FullConfig.model_validate(config_dict)


def test_alloy_generator_creates_correct_composition(alloy_config: FullConfig) -> None:
    """Test that the generator creates a structure with the correct atom counts."""
    generator = AlloyGenerator(config=alloy_config)

    # We patch the validation to isolate the generation logic for this test
    with patch.object(AlloyGenerator, "_validate_structure", return_value=True):
        structures = generator.generate()

    assert len(structures) > 0
    atoms = structures[0]

    # ASE's `bulk(..., cubic=True)` creates a 4-atom conventional cell.
    # A 2x2x2 supercell of this is 4 * 8 = 32 atoms.
    # Fe: 0.5 * 32 = 16 atoms, Pt: 0.5 * 32 = 16 atoms.
    expected_counts = {"Fe": 16, "Pt": 16}
    actual_counts = atoms.get_atomic_numbers().tolist()  # type: ignore[no-untyped-call]
    fe_count = actual_counts.count(26)  # Atomic number for Fe
    pt_count = actual_counts.count(78)  # Atomic number for Pt

    assert fe_count == expected_counts["Fe"]
    assert pt_count == expected_counts["Pt"]
    assert len(atoms) == sum(expected_counts.values())


def test_alloy_generator_validation_rejects_bad_structure(alloy_config: FullConfig) -> None:
    """Test that if validation fails, the structure is not returned."""
    generator = AlloyGenerator(config=alloy_config)
    with patch.object(AlloyGenerator, "_validate_structure", return_value=False):
        structures = generator.generate()
    assert len(structures) == 0
