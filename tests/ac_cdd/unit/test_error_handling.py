from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.graph_nodes import CycleNodes


@pytest.fixture
def mock_sandbox() -> MagicMock:
    return MagicMock()


@pytest.fixture
def mock_jules_client() -> AsyncMock:
    return AsyncMock()


@pytest.fixture
def mock_session_manager() -> AsyncMock:
    return AsyncMock()


@pytest.fixture
def mock_git_manager() -> AsyncMock:
    return AsyncMock()


@pytest.fixture
def cycle_nodes(
    mock_sandbox: MagicMock,
    mock_jules_client: AsyncMock,
    mock_session_manager: AsyncMock,
    mock_git_manager: AsyncMock,
) -> CycleNodes:
    return CycleNodes(
        sandbox=mock_sandbox,
        jules=mock_jules_client,
        session_manager=mock_session_manager,
        git_manager=mock_git_manager,
    )


@pytest.mark.asyncio
async def test_architect_node_handles_jules_exception(
    cycle_nodes: CycleNodes,
    mock_jules_client: AsyncMock,
    mock_session_manager: AsyncMock,
) -> None:
    """
    Verify architect_session_node handles exceptions from JulesClient,
    updates state to failed, and returns a failed status.
    """
    # Arrange
    error_message = "Jules API is down"
    mock_jules_client.run_session.side_effect = Exception(error_message)
    state = {"cycle_id": "01"}

    # Act
    result = await cycle_nodes.architect_session_node(state)

    # Assert
    assert result["status"] == "architect_failed"
    mock_session_manager.update_cycle_state.assert_awaited_once_with(
        "01", status="architect_failed", error=error_message
    )


@pytest.mark.asyncio
async def test_coder_node_handles_jules_exception(
    cycle_nodes: CycleNodes,
    mock_jules_client: AsyncMock,
    mock_session_manager: AsyncMock,
) -> None:
    """
    Verify coder_session_node handles exceptions from JulesClient,
    updates state to failed, and returns a failed status.
    """
    # Arrange
    error_message = "Jules API timed out"
    mock_jules_client.run_session.side_effect = Exception(error_message)
    # Simulate a non-resume scenario
    state = {"cycle_id": "02", "iteration_count": 1, "resume_mode": False}
    mock_session_manager.get_cycle.return_value.jules_session_id = None

    # Act
    result = await cycle_nodes.coder_session_node(state)

    # Assert
    assert result["status"] == "failed"
    assert result["error"] == error_message
    mock_session_manager.update_cycle_state.assert_awaited_once_with(
        "02", status="failed", error=error_message
    )
