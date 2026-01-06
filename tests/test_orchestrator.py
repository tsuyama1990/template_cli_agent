"""Unit tests for the pipeline orchestrator."""

from unittest.mock import MagicMock, patch

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


@patch("mlip_autopipec.core.orchestrator.AseDBWrapper")
@patch("mlip_autopipec.core.orchestrator.RandomSampler")
@patch("mlip_autopipec.core.orchestrator.MDExplorer")
@patch("mlip_autopipec.core.orchestrator.AlloyGenerator")
def test_pipeline_runner_orchestration(
    mock_alloy_generator, mock_md_explorer, mock_random_sampler, mock_ase_db_wrapper
) -> None:
    # Corresponds to SPEC.md, Section 2: Component Blueprint (`core/orchestrator.py`)
    # This test verifies that the `PipelineRunner` correctly instantiates and
    # calls each of the four pipeline stages in the exact required sequence.
    """Test that the PipelineRunner calls each component in the correct sequence."""
    # 1. Setup mock instances and return values
    mock_generator_instance = MagicMock()
    mock_generator_instance.generate.return_value = [MagicMock()]
    mock_alloy_generator.return_value = mock_generator_instance

    mock_explorer_instance = MagicMock()
    mock_explorer_instance.run_md.return_value = [MagicMock(), MagicMock()]
    mock_md_explorer.return_value = mock_explorer_instance

    mock_sampler_instance = MagicMock()
    mock_sampler_instance.sample.return_value = [MagicMock()]
    mock_random_sampler.return_value = mock_sampler_instance

    mock_db_instance = MagicMock()
    mock_db_context_manager = MagicMock()
    mock_db_context_manager.__enter__.return_value = mock_db_instance
    mock_ase_db_wrapper.return_value = mock_db_context_manager

    # 2. Create a valid configuration
    config_dict = {
        "project_name": "test_orchestration",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 1,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 1.0},
    }
    full_config = FullConfig(**config_dict)

    # 3. Instantiate and run the pipeline
    runner = PipelineRunner(config=full_config)
    runner.run()

    # 4. Assert calls and data flow
    mock_alloy_generator.assert_called_once_with(config=full_config.system)
    mock_generator_instance.generate.assert_called_once()

    mock_md_explorer.assert_called_once_with(config=full_config.exploration)
    mock_explorer_instance.run_md.assert_called_once_with(
        structures=mock_generator_instance.generate.return_value
    )

    mock_random_sampler.assert_called_once_with(config=full_config.sampling)
    mock_sampler_instance.sample.assert_called_once_with(
        structures=mock_explorer_instance.run_md.return_value
    )

    mock_ase_db_wrapper.assert_called_once_with(db_path="test_orchestration.db")
    mock_db_instance.write_structures.assert_called_once_with(
        structures=mock_sampler_instance.sample.return_value
    )
