from unittest.mock import MagicMock, call

import pytest
from ase import Atoms

from mlip_autopipec.data.models import StructureGeneration
from mlip_autopipec.modules.a_structure_generator import StructureGenerator


@pytest.fixture
def mock_db_wrapper():
    """Provides a mock AseDBWrapper with a callable `add_structures` method."""
    mock = MagicMock()
    mock.add_structures = MagicMock()
    return mock


@pytest.fixture
def structure_generation_config() -> StructureGeneration:
    """Provides a sample StructureGeneration config model."""
    return StructureGeneration(
        generation_strategy="sqs",
        supercell_size=8,
        strains=[-0.01, 0.0, 0.01],
    )


def test_structure_generator_execute(
    mock_db_wrapper: MagicMock, structure_generation_config: StructureGeneration
):
    """
    Tests that the StructureGenerator's execute method generates and saves
    the correct number and type of structures.
    """
    # 1. Arrange: Instantiate the generator with config and the mock DB
    generator = StructureGenerator(
        config=structure_generation_config, db_wrapper=mock_db_wrapper
    )

    # 2. Act: Run the execute method
    generator.execute()

    # 3. Assert:
    # Check that the database's add_structures method was called exactly once
    mock_db_wrapper.add_structures.assert_called_once()

    # Check the content of the call
    # Get the list of Atoms objects that were passed to add_structures
    args, _ = mock_db_wrapper.add_structures.call_args
    added_structures = args[0]

    # Assert that the number of structures is correct (should equal the number of strains)
    assert len(added_structures) == len(structure_generation_config.strains)

    # Assert that each structure is an ASE Atoms object
    for structure in added_structures:
        assert isinstance(structure, Atoms)

    # Assert that the strains were applied correctly by checking cell volumes
    # (assuming a simple 1x1x1 cell for the base atom)
    base_volume = 1.0
    expected_volumes = [
        pytest.approx(base_volume * (1 + strain) ** 3)
        for strain in structure_generation_config.strains
    ]
    actual_volumes = sorted([s.get_volume() for s in added_structures])
    expected_volumes.sort(key=lambda x: x.expected)

    # Compare the sorted lists
    assert actual_volumes == expected_volumes
