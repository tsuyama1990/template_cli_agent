"""Unit tests for the pipeline orchestrator."""
from unittest.mock import MagicMock, patch

import pytest

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


@pytest.fixture
def mock_full_config() -> FullConfig:
    """Return a mock FullConfig object."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
    }
    return FullConfig(**config_dict)


@patch("mlip_autopipec.core.orchestrator.AseDBWrapper")
@patch("mlip_autopipec.core.orchestrator.RandomSampler")
@patch("mlip_autopipec.core.orchestrator.MDExplorer")
@patch("mlip_autopipec.core.orchestrator.AlloyGenerator")
def test_pipeline_runner_execution_flow(
    mock_alloy_generator: MagicMock,
    mock_md_explorer: MagicMock,
    mock_random_sampler: MagicMock,
    mock_db_wrapper: MagicMock,
    mock_full_config: FullConfig,
) -> None:
    """Test that the pipeline runner calls components in the correct sequence."""
    # Mock the instances and their methods' return values
    mock_gen_instance = mock_alloy_generator.return_value
    mock_gen_instance.generate.return_value = [MagicMock()] * 10

    mock_exp_instance = mock_md_explorer.return_value
    mock_exp_instance.run_md.return_value = [MagicMock()] * 10

    mock_samp_instance = mock_random_sampler.return_value
    mock_samp_instance.sample.return_value = [MagicMock()] * 8

    mock_db_instance = mock_db_wrapper.return_value
    mock_db_instance.__enter__.return_value = mock_db_instance

    # Instantiate and run the pipeline
    runner = PipelineRunner(mock_full_config)
    runner.run()

    # Assert that components were instantiated correctly
    mock_alloy_generator.assert_called_once_with(mock_full_config.system)
    mock_md_explorer.assert_called_once()
    mock_random_sampler.assert_called_once_with(mock_full_config.sampling)
    mock_db_wrapper.assert_called_once()

    # Assert that the methods were called in the correct order
    mock_gen_instance.generate.assert_called_once()
    mock_exp_instance.run_md.assert_called_once_with(
        mock_gen_instance.generate.return_value
    )
    mock_samp_instance.sample.assert_called_once_with(
        mock_exp_instance.run_md.return_value
    )
    mock_db_instance.write_structures.assert_called_once_with(
        mock_samp_instance.sample.return_value
    )
