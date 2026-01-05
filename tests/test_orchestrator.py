"""Unit tests for the PipelineRunner orchestrator."""
from unittest.mock import patch, MagicMock
import pytest
import yaml

from mlip_autopipec.core.orchestrator import PipelineRunner
from mlip_autopipec.config.models import FullConfig

# A valid config to instantiate the PipelineRunner
VALID_CONFIG_YAML = """
project_name: test_orchestrator
system:
  elements: ['Fe', 'Pt']
  composition: {'Fe': 0.5, 'Pt': 0.5}
  lattice: 'fcc'
  num_structures: 10
exploration:
  temperature: 300.0
sampling:
  method: 'random'
  fraction: 0.8
"""

@pytest.fixture
def valid_config() -> FullConfig:
    """Provides a valid FullConfig object for testing."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    return FullConfig(**config_dict)

@patch('mlip_autopipec.core.orchestrator.AseDBWrapper')
@patch('mlip_autopipec.core.orchestrator.RandomSampler')
@patch('mlip_autopipec.core.orchestrator.MDExplorer')
@patch('mlip_autopipec.core.orchestrator.AlloyGenerator')
def test_pipeline_runner_sequence(
    mock_alloy_gen: MagicMock,
    mock_md_explorer: MagicMock,
    mock_sampler: MagicMock,
    mock_db_wrapper: MagicMock,
    valid_config: FullConfig
):
    """
    Verify that the PipelineRunner calls each component in the correct order
    and passes data between them correctly.
    """
    # Arrange: Set up mock return values to trace the data flow
    mock_initial_structs = [MagicMock(name="initial_struct")]
    mock_explored_structs = [MagicMock(name="explored_struct")]
    mock_sampled_structs = [MagicMock(name="sampled_struct")]

    # Configure the return values of the mocked instances' methods
    mock_alloy_gen.return_value.generate.return_value = mock_initial_structs
    mock_md_explorer.return_value.run_md.return_value = mock_explored_structs
    mock_sampler.return_value.sample.return_value = mock_sampled_structs

    # The AseDBWrapper is a context manager, so we need to mock that behavior
    mock_db_instance = mock_db_wrapper.return_value.__enter__.return_value

    # Act: Instantiate and run the pipeline
    runner = PipelineRunner(config=valid_config)
    runner.run()

    # Assert: Verify that each component was instantiated and called correctly

    # 1. Generator
    mock_alloy_gen.assert_called_once_with(valid_config.system)
    mock_alloy_gen.return_value.generate.assert_called_once()

    # 2. Explorer
    mock_md_explorer.assert_called_once()
    mock_md_explorer.return_value.run_md.assert_called_once_with(mock_initial_structs)

    # 3. Sampler
    mock_sampler.assert_called_once_with(valid_config.sampling)
    mock_sampler.return_value.sample.assert_called_once_with(mock_explored_structs)

    # 4. Storage
    mock_db_wrapper.assert_called_once_with(runner.db_path)
    mock_db_instance.write_structures.assert_called_once_with(mock_sampled_structs)
