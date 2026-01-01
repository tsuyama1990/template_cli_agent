"""Unit tests for the MoleculeGenerator module."""

import numpy as np
import pytest

from mlip_autopipec.config import FullConfig, ModelType
from mlip_autopipec.modules.molecule_generator import MoleculeGenerator


@pytest.fixture
def molecule_config() -> FullConfig:
    """Provides a valid FullConfig for a molecular system."""
    config_dict = {
        "db_path": "test.db",
        "qe_command": "pw.x",
        "system": {
            "elements": ["H", "O"],
            "composition": "H2O",
            "structure_type": "molecule",
            "melting_point_guess": 273.0,
        },
        "simulation": {
            "temperature_steps": [300],
            "initial_structures_to_generate": 8,
        },
        "dft_compute": {
            "ecutwfc": 80.0,
            "ecutrho": 320.0,
            "kpoints_density": 1.0, # Not relevant for molecules
            "control": {"calculation": "scf"},
            "pseudopotentials": {"H": "H.UPF", "O": "O.UPF"},
        },
        "mlip_training": {
            "model_type": ModelType.ACE,
            "r_cut": 4.0,
            "loss_weights": {"energy": 1.0, "forces": 100.0},
        },
    }
    return FullConfig.model_validate(config_dict)


def test_molecule_generator_initialization(molecule_config):
    """Test that the MoleculeGenerator initializes correctly."""
    generator = MoleculeGenerator(molecule_config)
    assert generator.config == molecule_config


def test_molecule_generator_returns_correct_number_of_structures(molecule_config):
    """Test that generate() returns the correct number of atoms objects."""
    generator = MoleculeGenerator(molecule_config)
    structures = generator.generate()
    assert len(structures) == molecule_config.simulation.initial_structures_to_generate


def test_molecule_generator_creates_correct_base_molecule(molecule_config):
    """Test that the first generated structure is the correct equilibrium molecule."""
    generator = MoleculeGenerator(molecule_config)
    structures = generator.generate()
    base_structure = structures[0]

    assert base_structure.get_chemical_formula() == "H2O"
    assert len(base_structure) == 3


def test_molecule_generator_applies_displacements(molecule_config):
    """Test that subsequent structures are distorted."""
    generator = MoleculeGenerator(molecule_config)
    structures = generator.generate()
    base_structure = structures[0]

    # Check that other structures are indeed different
    assert any(
        np.any(np.abs(s.positions - base_structure.positions) > 1e-6) for s in structures[1:]
    ), "Expected distorted structures with different positions."

def test_molecule_generator_raises_for_invalid_formula(molecule_config):
    """Test that an error is raised for an unknown molecular formula."""
    molecule_config.system.composition = "InvalidFormula123"
    generator = MoleculeGenerator(molecule_config)
    with pytest.raises(ValueError, match="Could not create molecule"):
        generator.generate()
