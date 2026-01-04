from pathlib import Path
from unittest.mock import MagicMock

import pytest
from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.interfaces import IExplorer, ISampler, IStructureGenerator
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator
from mlip_autopipec.storage.database_manager import DatabaseManager


@pytest.fixture
def mock_generator() -> MagicMock:
    """Provides a mock IStructureGenerator that returns a single Atoms object."""
    mock = MagicMock()
    mock.generate.return_value = [Atoms("H")]
    return mock


@pytest.fixture
def mock_explorer() -> MagicMock:
    """Provides a mock IExplorer that returns a single Atoms object."""
    mock = MagicMock()
    mock.run.return_value = [Atoms("He")]
    return mock


@pytest.fixture
def mock_sampler() -> MagicMock:
    """Provides a mock ISampler that returns a single Atoms object."""
    mock = MagicMock()
    mock.sample.return_value = [Atoms("Li")]
    return mock


@pytest.fixture
def mock_db_manager() -> MagicMock:
    """Provides a mock DatabaseManager."""
    return MagicMock()


def test_orchestrator_run_sequence(
    mock_config: FullConfig,
    mock_generator: IStructureGenerator,
    mock_explorer: IExplorer,
    mock_sampler: ISampler,
    mock_db_manager: DatabaseManager,
    tmp_path: Path,
) -> None:
    """
    Tests that the orchestrator calls each component in the correct sequence.
    """
    # Arrange
    orchestrator = PipelineOrchestrator(
        config=mock_config,
        generator=mock_generator,
        explorer=mock_explorer,
        sampler=mock_sampler,
        db_manager=mock_db_manager,
        output_dir=tmp_path,
    )

    # Act
    orchestrator.run()

    # Assert
    mock_generator.generate.assert_called_once()
    mock_explorer.run.assert_called_once()
    mock_sampler.sample.assert_called_once()
    mock_db_manager.write_atoms.assert_called_once()
