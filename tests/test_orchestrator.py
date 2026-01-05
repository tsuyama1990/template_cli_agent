"""Unit tests for the PipelineRunner."""
from unittest.mock import MagicMock, patch

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


def test_pipeline_runner_orchestration():
    """Test that the PipelineRunner calls components in the correct order."""
    config = FullConfig(
        project_name="test_project",
        system={
            "elements": ["Fe", "Pt"],
            "composition": {"Fe": 0.5, "Pt": 0.5},
            "lattice": "fcc",
            "num_structures": 10,
        },
        exploration={"temperature": 300.0},
        sampling={"method": "random", "fraction": 0.8},
    )

    with patch("mlip_autopipec.core.orchestrator.AlloyGenerator") as mock_generator, patch(
        "mlip_autopipec.core.orchestrator.MDExplorer"
    ) as mock_explorer, patch(
        "mlip_autopipec.core.orchestrator.RandomSampler"
    ) as mock_sampler, patch(
        "mlip_autopipec.core.orchestrator.AseDBWrapper"
    ) as mock_db_wrapper:
        mock_gen_instance = MagicMock()
        mock_generator.return_value = mock_gen_instance
        mock_gen_instance.generate.return_value = [MagicMock()] * 10

        mock_exp_instance = MagicMock()
        mock_explorer.return_value = mock_exp_instance
        mock_exp_instance.run_md.return_value = [MagicMock()] * 10

        mock_samp_instance = MagicMock()
        mock_sampler.return_value = mock_samp_instance
        mock_samp_instance.sample.return_value = [MagicMock()] * 8

        mock_db_instance = MagicMock()
        mock_db_wrapper.return_value = mock_db_instance

        runner = PipelineRunner(config)
        runner.run()

        mock_generator.assert_called_once_with(config.system)
        mock_gen_instance.generate.assert_called_once()

        mock_explorer.assert_called_once()
        mock_exp_instance.run_md.assert_called_once()

        mock_sampler.assert_called_once_with(config.sampling)
        mock_samp_instance.sample.assert_called_once()

        mock_db_wrapper.assert_called_once_with("test_project.db")
        mock_db_instance.write_structures.assert_called_once()
