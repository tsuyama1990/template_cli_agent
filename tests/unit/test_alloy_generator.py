"""Unit tests for the AlloyGenerator module."""

import pytest

from mlip_autopipec.config import FullConfig, ModelType
from mlip_autopipec.constants import ALLOY_TARGET_ATOMS
from mlip_autopipec.modules.alloy_generator import AlloyGenerator


@pytest.fixture
def alloy_config() -> FullConfig:
    """Provides a valid FullConfig for an alloy system."""
    config_dict = {
        "db_path": "test.db",
        "qe_command": "pw.x",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": "FePt",
            "structure_type": "alloy",
            "melting_point_guess": 1600.0,
        },
        "simulation": {
            "temperature_steps": [300, 600, 900],
            "initial_structures_to_generate": 5,
        },
        "explorer": {"surrogate_model": "mace_mp"},
        "dft_compute": {
            "ecutwfc": 60.0,
            "ecutrho": 240.0,
            "kpoints_density": 5.0,
            "magnetism": "ferromagnetic",
            "control": {"calculation": "scf"},
            "pseudopotentials": {"Fe": "Fe.UPF", "Pt": "Pt.UPF"},
        },
        "mlip_training": {
            "model_type": ModelType.ACE,
            "r_cut": 5.0,
            "loss_weights": {"energy": 1.0, "forces": 100.0},
        },
    }
    return FullConfig.model_validate(config_dict)


def test_alloy_generator_initialization(alloy_config):
    """Test that the AlloyGenerator initializes correctly."""
    generator = AlloyGenerator(alloy_config)
    assert generator.config == alloy_config


def test_alloy_generator_returns_correct_number_of_structures(alloy_config):
    """Test that generate() returns the correct number of atoms objects."""
    generator = AlloyGenerator(alloy_config)
    structures = generator.generate()
    assert len(structures) == alloy_config.simulation.initial_structures_to_generate


def test_alloy_generator_creates_correct_composition(alloy_config):
    """Test that the generated structures have the correct composition."""
    generator = AlloyGenerator(alloy_config)
    structures = generator.generate()
    base_structure = structures[0]

    expected_atoms = ALLOY_TARGET_ATOMS
    assert len(base_structure) == expected_atoms

    symbols = base_structure.get_chemical_symbols()
    assert symbols.count("Fe") == expected_atoms / 2
    assert symbols.count("Pt") == expected_atoms / 2


def test_alloy_generator_applies_strains(alloy_config):
    """Test that the generated structures include strained variations."""
    generator = AlloyGenerator(alloy_config)
    structures = generator.generate()
    base_structure = structures[0]

    assert any(
        abs(s.get_volume() - base_structure.get_volume()) > 1e-6 for s in structures[1:]
    ), "Expected strained structures with different volumes."
