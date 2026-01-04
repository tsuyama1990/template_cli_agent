from unittest.mock import patch

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator_returns_correct_number_of_structures(
    mock_config: FullConfig,
) -> None:
    """
    Tests that the AlloyGenerator produces the expected number of Atoms objects.
    """
    # Arrange
    generator = AlloyGenerator(
        system_config=mock_config.system, generation_config=mock_config.generation
    )

    # Act
    structures = generator.generate()

    # Assert
    assert len(structures) == mock_config.system.num_initial_structures


def test_alloy_generator_validates_structure(mock_config: FullConfig) -> None:
    """
    Tests that the _validate_structure method is called for each generated structure.
    """
    # Arrange
    generator = AlloyGenerator(
        system_config=mock_config.system, generation_config=mock_config.generation
    )

    # Act
    with patch.object(
        generator, "_validate_structure", return_value=True
    ) as mock_validate:
        generator.generate()

        # Assert
        assert mock_validate.call_count == mock_config.system.num_initial_structures
