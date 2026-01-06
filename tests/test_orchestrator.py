# -*- coding: utf-8 -*-
"""Unit tests for the PipelineRunner orchestrator."""
from unittest.mock import MagicMock, patch

import pytest

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


@patch("mlip_autopipec.core.orchestrator.AseDBWrapper")
@patch("mlip_autopipec.core.orchestrator.RandomSampler")
@patch("mlip_autopipec.core.orchestrator.MDExplorer")
@patch("mlip_autopipec.core.orchestrator.AlloyGenerator")
def test_pipeline_runner_orchestration(
    mock_alloy_generator, mock_md_explorer, mock_random_sampler, mock_ase_db_wrapper,
    mocker
):
    """Verify that the PipelineRunner calls all components in the correct order."""
    # 1. Setup Mocks
    mock_generator_instance = MagicMock()
    mock_explorer_instance = MagicMock()
    mock_sampler_instance = MagicMock()
    mock_db_instance = MagicMock()

    mock_alloy_generator.return_value = mock_generator_instance
    mock_md_explorer.return_value = mock_explorer_instance
    mock_random_sampler.return_value = mock_sampler_instance
    mock_ase_db_wrapper.return_value.__enter__.return_value = mock_db_instance

    mock_generated_structures = [MagicMock()]
    mock_generator_instance.generate.return_value = mock_generated_structures
    mock_explorer_instance.run_md.return_value = mock_generated_structures
    mock_sampler_instance.sample.return_value = mock_generated_structures

    # 2. Prepare Config
    config = FullConfig.model_validate({
        "project_name": "test", "system": {"elements": ["Cu"], "composition": {"Cu": 1.0}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300}, "sampling": {"method": "random", "fraction": 1.0}
    })

    # 3. Run Pipeline
    runner = PipelineRunner(config)
    runner.run()

    # 4. Assert Calls
    mock_alloy_generator.assert_called_once_with(config=config.system)
    mock_generator_instance.generate.assert_called_once()
    mock_md_explorer.assert_called_once()
    mock_explorer_instance.run_md.assert_called_once_with(mock_generated_structures)
    mock_random_sampler.assert_called_once_with(config=config.sampling)
    mock_sampler_instance.sample.assert_called_once_with(mock_generated_structures)
    mock_ase_db_wrapper.assert_called_once()
    mock_db_instance.write_structures.assert_called_once_with(mock_generated_structures)


@patch("mlip_autopipec.core.orchestrator.AlloyGenerator")
def test_pipeline_runner_error_handling(mock_alloy_generator, caplog):
    """Test that exceptions during a pipeline stage are caught and logged."""
    # 1. Setup Mock to raise an exception
    mock_generator_instance = MagicMock()
    mock_generator_instance.generate.side_effect = ValueError("Generation failed")
    mock_alloy_generator.return_value = mock_generator_instance

    # 2. Prepare Config
    config = FullConfig.model_validate({
        "project_name": "test_error", "system": {"elements": ["Cu"], "composition": {"Cu": 1.0}, "lattice": "fcc", "num_structures": 1},
        "exploration": {"temperature": 300}, "sampling": {"method": "random", "fraction": 1.0}
    })

    # 3. Run and Assert
    runner = PipelineRunner(config)
    with pytest.raises(ValueError, match="Generation failed"):
        runner.run()

    # Check that the error was logged
    assert "An unexpected error occurred" in caplog.text
    assert "Generation failed" in caplog.text
