from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.domain_models import CycleManifest
from ac_cdd_core.graph_nodes import CycleNodes
from ac_cdd_core.services.coder_service import CoderService


@pytest.fixture
def mock_coder_service() -> CoderService:
    return MagicMock(spec=CoderService)


@pytest.fixture
def cycle_nodes(mock_coder_service: CoderService) -> CycleNodes:
    return CycleNodes(MagicMock(), MagicMock(), mock_coder_service)


@pytest.mark.asyncio
async def test_coder_session_node_delegates_to_service(
    cycle_nodes: CycleNodes, mock_coder_service: CoderService
) -> None:
    """Test that the coder_session_node correctly delegates to the CoderService."""
    # Setup
    mock_coder_service.run_coder_session = AsyncMock(
        return_value={"status": "ready_for_audit"}
    )
    state = {"cycle_id": "01", "iteration_count": 1, "resume_mode": True}

    # Execute
    result = await cycle_nodes.coder_session_node(state)

    # Verify
    mock_coder_service.run_coder_session.assert_awaited_once_with("01", True)
    assert result["status"] == "ready_for_audit"
