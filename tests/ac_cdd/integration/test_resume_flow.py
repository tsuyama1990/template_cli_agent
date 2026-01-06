from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.domain_models import CycleManifest, JulesSessionResult
from ac_cdd_core.graph_nodes import CycleNodes
from ac_cdd_core.services.jules_client import JulesClient


@pytest.fixture
def mock_sandbox_runner() -> MagicMock:
    """Provides a mock SandboxRunner."""
    return MagicMock()


@pytest.fixture
def mock_jules_client() -> MagicMock:
    """Provides a mock JulesClient."""
    mock = MagicMock(spec=JulesClient)
    # Ensure async methods are awaitable
    mock.continue_session = AsyncMock()
    mock.run_session = AsyncMock()
    mock.wait_for_completion = AsyncMock()
    return mock


@pytest.fixture
def mock_session_manager() -> AsyncMock:
    """Provides a mock SessionManager."""
    return AsyncMock()


@pytest.fixture
def mock_git_manager() -> AsyncMock:
    """Provides a mock GitManager."""
    return AsyncMock()


@pytest.fixture
def cycle_nodes(
    mock_sandbox_runner: MagicMock,
    mock_jules_client: MagicMock,
    mock_session_manager: AsyncMock,
    mock_git_manager: AsyncMock,
) -> CycleNodes:
    """Fixture to create CycleNodes with mocked dependencies."""
    return CycleNodes(
        sandbox=mock_sandbox_runner,
        jules=mock_jules_client,
        session_manager=mock_session_manager,
        git_manager=mock_git_manager,
    )


@pytest.mark.asyncio
async def test_resume_logic_hot_resume(
    cycle_nodes: CycleNodes,
    mock_jules_client: MagicMock,
    mock_session_manager: AsyncMock,
) -> None:
    """Test hot resume flow where a session ID already exists."""
    # Arrange
    cycle_manifest = CycleManifest(id="01", jules_session_id="existing-session-123")
    mock_session_manager.get_cycle.return_value = cycle_manifest
    mock_jules_client.continue_session.return_value = JulesSessionResult(
        status="success",
        session_name="existing-session-123",
        pr_url="http://github.com/pr/hot-resume",
    )
    state = {"cycle_id": "01", "resume_mode": True}

    # Act
    result = await cycle_nodes.coder_session_node(state)

    # Assert
    mock_jules_client.continue_session.assert_awaited_once_with("existing-session-123")
    mock_jules_client.run_session.assert_not_called()
    mock_session_manager.update_cycle_state.assert_awaited_with(
        "01", pr_url="http://github.com/pr/hot-resume", status="completed"
    )
    assert result["status"] == "ready_for_audit"


@pytest.mark.asyncio
async def test_resume_logic_cold_start_persists_id(
    cycle_nodes: CycleNodes,
    mock_jules_client: MagicMock,
    mock_session_manager: AsyncMock,
) -> None:
    """Test cold start flow where no session ID exists, and it gets persisted."""
    # Arrange
    cycle_manifest = CycleManifest(id="01", jules_session_id=None)
    mock_session_manager.get_cycle.return_value = cycle_manifest
    mock_jules_client.run_session.return_value = JulesSessionResult(
        status="running", session_name="new-session-456", pr_url=None
    )
    mock_jules_client.wait_for_completion.return_value = JulesSessionResult(
        status="success",
        session_name="new-session-456",
        pr_url="http://github.com/pr/cold-start",
    )
    state = {"cycle_id": "01", "resume_mode": False}

    # Act
    result = await cycle_nodes.coder_session_node(state)

    # Assert
    mock_jules_client.run_session.assert_awaited_once()
    # Check that the new session ID was persisted
    mock_session_manager.update_cycle_state.assert_any_await(
        "01", jules_session_id="new-session-456", status="in_progress"
    )
    # Check that the final state was persisted after completion
    mock_session_manager.update_cycle_state.assert_awaited_with(
        "01", pr_url="http://github.com/pr/cold-start", status="completed"
    )
    assert result["status"] == "ready_for_audit"
