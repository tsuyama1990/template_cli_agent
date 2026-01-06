from unittest.mock import MagicMock, patch

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


@patch("mlip_autopipec.core.orchestrator.AlloyGenerator")
@patch("mlip_autopipec.core.orchestrator.MDExplorer")
@patch("mlip_autopipec.core.orchestrator.RandomSampler")
@patch("mlip_autopipec.core.orchestrator.AseDBWrapper")
def test_pipeline_runner(
    mock_storage: MagicMock,
    mock_sampler: MagicMock,
    mock_explorer: MagicMock,
    mock_generator: MagicMock,
):
    """Tests that the PipelineRunner calls components in the correct order."""
    config_dict = {
        "project_name": "test_project",
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.5},
    }
    config = FullConfig(**config_dict)

    # Mock the return values of the components
    mock_generator.return_value.generate.return_value = [MagicMock()] * 10
    mock_explorer.return_value.run_md.return_value = [MagicMock()] * 10
    mock_sampler.return_value.sample.return_value = [MagicMock()] * 5

    runner = PipelineRunner(config)
    runner.run()

    # Assert that each component was called once
    mock_generator.assert_called_once()
    mock_explorer.assert_called_once()
    mock_sampler.assert_called_once()
    mock_storage.assert_called_once()

    # Assert that the methods were called
    mock_generator.return_value.generate.assert_called_once()
    mock_explorer.return_value.run_md.assert_called_once()
    mock_sampler.return_value.sample.assert_called_once()
    mock_storage.return_value.write_structures.assert_called_once()
