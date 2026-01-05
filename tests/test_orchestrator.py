"""Unit tests for the pipeline orchestrator."""

from unittest.mock import MagicMock

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


def test_pipeline_runner(mocker):
    """
    Test that the PipelineRunner calls all components in the correct order.
    """
    # Mock all the components
    mocker.patch(
        "mlip_autopipec.core.orchestrator.AlloyGenerator",
        return_value=MagicMock(),
    )
    mocker.patch(
        "mlip_autopipec.core.orchestrator.MDExplorer",
        return_value=MagicMock(),
    )
    mocker.patch(
        "mlip_autopipec.core.orchestrator.RandomSampler",
        return_value=MagicMock(),
    )
    mock_write_structures = mocker.patch("mlip_autopipec.core.orchestrator.write_structures")

    # Create a valid config
    valid_data = {
        "system": {
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        "exploration": {"temperature": 300.0},
        "sampling": {"method": "random", "fraction": 0.8},
        "project_name": "test_project",
    }
    config = FullConfig(**valid_data)

    # Run the pipeline
    runner = PipelineRunner(config)
    runner.run()

    # Assert that the components were called
    assert PipelineRunner.__module__ + ".AlloyGenerator.generate.called_once"
    assert PipelineRunner.__module__ + ".MDExplorer.run_md.called_once"
    assert PipelineRunner.__module__ + ".RandomSampler.sample.called_once"
    mock_write_structures.assert_called_once()
