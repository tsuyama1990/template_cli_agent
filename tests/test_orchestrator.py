from unittest.mock import MagicMock, patch

from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner


def test_pipeline_runner_orchestration() -> None:
    """Test that the PipelineRunner calls components in the correct order."""
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
    config = FullConfig.model_validate(config_dict)

    with patch("mlip_autopipec.core.orchestrator.AlloyGenerator") as mock_generator, \
         patch("mlip_autopipec.core.orchestrator.MDExplorer") as mock_explorer, \
         patch("mlip_autopipec.core.orchestrator.RandomSampler") as mock_sampler, \
         patch("mlip_autopipec.core.orchestrator.AseDBWrapper") as mock_db_wrapper:

        mock_gen_instance = MagicMock()
        mock_exp_instance = MagicMock()
        mock_sam_instance = MagicMock()
        mock_db_instance = MagicMock()

        mock_generator.return_value = mock_gen_instance
        mock_explorer.return_value = mock_exp_instance
        mock_sampler.return_value = mock_sam_instance
        mock_db_wrapper.return_value.__enter__.return_value = mock_db_instance

        runner = PipelineRunner(config)
        runner.run()

        mock_gen_instance.generate.assert_called_once()
        mock_exp_instance.explore.assert_called_once()
        mock_sam_instance.sample.assert_called_once()
        mock_db_instance.write_structures.assert_called_once()
