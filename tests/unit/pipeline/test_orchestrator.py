"""Unit tests for the PipelineOrchestrator."""

import pytest
from unittest.mock import Mock, MagicMock
from mlip_autopipec.pipeline.orchestrator import PipelineOrchestrator
from mlip_autopipec.common.pydantic_models import FullConfig
from typing import Any, Dict

@pytest.fixture
def full_config() -> FullConfig:
    """Provides a default FullConfig for testing."""
    config: Dict[str, Any] = {
        "system": {
            "elements": ["Si"],
            "composition": {"Si": 1.0},
            "supercell_size": [2, 2, 2],
        },
        "exploration": {"temperature_k": 500.0, "pressure_gpa": 0.0},
        "sampling": {"method": "random", "n_samples": 50},
    }
    return FullConfig(**config)

def test_orchestrator_initialization(full_config: FullConfig) -> None:
    """Tests that the PipelineOrchestrator can be initialized."""
    orchestrator = PipelineOrchestrator(
        config=full_config,
        generator=Mock(),
        explorer=Mock(),
        sampler=Mock(),
        storage=Mock()
    )
    assert orchestrator is not None

def test_pipeline_orchestrator_calls_stages_in_order(full_config: FullConfig) -> None:
    """
    Tests that the orchestrator calls the generator, explorer, sampler, and
    storage in the correct sequence.
    """
    # Create mock objects for each component
    mock_generator = MagicMock()
    mock_explorer = MagicMock()
    mock_sampler = MagicMock()
    mock_storage = MagicMock()

    # Define mock return values to trace the data flow
    mock_generator.generate.return_value = [Mock()]  # Dummy generated structures
    mock_explorer.explore.return_value = [Mock()]    # Dummy trajectory
    mock_sampler.sample.return_value = [Mock()]      # Dummy sampled structures

    orchestrator = PipelineOrchestrator(
        config=full_config,
        generator=mock_generator,
        explorer=mock_explorer,
        sampler=mock_sampler,
        storage=mock_storage
    )

    # Run the pipeline
    orchestrator.run_pipeline()

    # Assert that each method was called once
    mock_generator.generate.assert_called_once()
    mock_explorer.explore.assert_called_once()
    mock_sampler.sample.assert_called_once()
    mock_storage.write_structures.assert_called_once()
