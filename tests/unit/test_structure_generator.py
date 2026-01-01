"""Unit tests for the StructureGenerator module."""

from unittest.mock import MagicMock, patch

import pytest
from ase import Atoms

from mlip_autopipec.config import FullConfig, ModelType
from mlip_autopipec.modules.structure_generator import StructureGenerator


# A minimal valid FullConfig for testing
@pytest.fixture
def mock_config() -> FullConfig:
    """Provides a mock FullConfig object that satisfies the schema."""
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


@pytest.fixture
def mock_db_wrapper() -> MagicMock:
    """Provides a mock AseDBWrapper."""
    return MagicMock()


def test_structure_generator_initialization(mock_config, mock_db_wrapper):
    """Test that the StructureGenerator initializes correctly."""
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    assert generator.config == mock_config
    assert generator.db == mock_db_wrapper


def test_generate_skips_if_db_not_empty(mock_config, mock_db_wrapper):
    """Test that structure generation is skipped if the database is not empty."""
    mock_db_wrapper.is_empty.return_value = False
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    generator.generate()
    mock_db_wrapper.add_atoms.assert_not_called()


@patch("mlip_autopipec.modules.structure_generator.StructureGenerator._dispatcher")
def test_generate_calls_dispatcher_if_db_empty(mock_dispatcher, mock_config, mock_db_wrapper):
    """Test that the dispatcher is called if the database is empty."""
    mock_db_wrapper.is_empty.return_value = True
    mock_dispatcher.return_value = [Atoms("H")]  # Return a dummy structure
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    generator.generate()
    mock_dispatcher.assert_called_once()
    mock_db_wrapper.add_atoms.assert_called_once()


def test_dispatcher_calls_sqs_for_alloy(mock_config, mock_db_wrapper):
    """Test the dispatcher calls the correct method for alloys."""
    mock_config.system.structure_type = "alloy"
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    with patch.object(generator, "_generate_sqs_for_alloy", return_value=[]) as mock_method:
        generator._dispatcher()
        mock_method.assert_called_once()


def test_dispatcher_calls_nms_for_molecule(mock_config, mock_db_wrapper):
    """Test the dispatcher calls the correct method for molecules."""
    mock_config.system.structure_type = "molecule"
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    with patch.object(generator, "_generate_nms_for_molecule", return_value=[]) as mock_method:
        generator._dispatcher()
        mock_method.assert_called_once()


def test_dispatcher_raises_for_unknown_type(mock_config, mock_db_wrapper):
    """Test the dispatcher raises a ValueError for an unknown structure type."""
    mock_config.system.structure_type = "unknown_type"
    generator = StructureGenerator(mock_config, mock_db_wrapper)
    with pytest.raises(ValueError, match="Unknown structure type: unknown_type"):
        generator._dispatcher()


def test_generate_sqs_for_alloy_creates_structures_ase(mock_config, mock_db_wrapper):
    """Test the ASE-based random alloy generation logic."""
    mock_config.system.elements = ["Fe", "Pt"]
    mock_config.system.composition = "FePt"
    mock_config.system.structure_type = "alloy"

    generator = StructureGenerator(mock_config, mock_db_wrapper)
    structures = generator._generate_sqs_for_alloy()

    # Check the generated structures
    assert len(structures) == mock_config.simulation.initial_structures_to_generate
    base_structure = structures[0]

    # Verify composition
    composition = base_structure.get_chemical_symbols()
    assert composition.count("Fe") == 32
    assert composition.count("Pt") == 32
    assert "Fe32Pt32" in base_structure.get_chemical_formula()


    # Check that strained structures were also created
    assert any(
        abs(s.get_volume() - base_structure.get_volume()) > 1e-6 for s in structures[1:]
    ), "Expected strained structures with different volumes."

# This is a placeholder for a more complex NMS test.
# A full test would require mocking geometry optimization and Hessian calculations.
def test_generate_nms_for_molecule_placeholder(mock_config, mock_db_wrapper):
    """Placeholder test for NMS generation."""
    mock_config.system.structure_type = "molecule"
    mock_config.system.composition = "H2O"
    mock_config.system.elements = ["H", "O"]

    generator = StructureGenerator(mock_config, mock_db_wrapper)
    # The actual implementation is not done, so we expect a warning and empty list
    with pytest.warns(UserWarning, match="NMS generation is not yet implemented"):
        structures = generator._generate_nms_for_molecule()

    assert structures == []
