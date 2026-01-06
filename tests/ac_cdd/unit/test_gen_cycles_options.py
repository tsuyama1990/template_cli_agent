from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.domain_models import JulesSessionResult
from ac_cdd_core.graph_nodes import CycleNodes
from ac_cdd_core.services.jules_client import JulesClient
from ac_cdd_core.state import CycleState
from rich.console import Console


@pytest.fixture
def mock_jules_client() -> AsyncMock:
    return AsyncMock(spec=JulesClient)


@pytest.fixture
def mock_session_manager() -> AsyncMock:
    return AsyncMock()


@pytest.fixture
def mock_git_manager() -> AsyncMock:
    return AsyncMock()


@pytest.fixture
def mock_sandbox() -> MagicMock:
    return MagicMock()


@pytest.fixture
def mock_console() -> MagicMock:
    return MagicMock(spec=Console)


@pytest.mark.asyncio
class TestGenCyclesCountOption:
    async def test_prompt_injection_with_count(
        self,
        mock_sandbox: MagicMock,
        mock_jules_client: AsyncMock,
        mock_session_manager: AsyncMock,
        mock_git_manager: AsyncMock,
        mock_console: MagicMock,
    ) -> None:
        """Verify max_cycles is injected into the prompt when --count is used."""
        # Arrange
        nodes = CycleNodes(mock_sandbox, mock_jules_client, mock_session_manager, mock_git_manager)
        state = CycleState(cycle_id="01", requested_cycle_count=5)
        mock_jules_client.run_session.return_value = JulesSessionResult(
            status="success", session_name="test", pr_url="http://pr"
        )

        # Act
        await nodes.architect_session_node(state.model_dump())

        # Assert
        mock_jules_client.run_session.assert_awaited_once()
        call_kwargs = mock_jules_client.run_session.await_args.kwargs
        assert "max_cycles=5" in call_kwargs["user_prompt"]

    async def test_prompt_no_injection_without_count(
        self,
        mock_sandbox: MagicMock,
        mock_jules_client: AsyncMock,
        mock_session_manager: AsyncMock,
        mock_git_manager: AsyncMock,
        mock_console: MagicMock,
    ) -> None:
        """Verify max_cycles is NOT in the prompt when --count is omitted."""
        # Arrange
        nodes = CycleNodes(mock_sandbox, mock_jules_client, mock_session_manager, mock_git_manager)
        state = CycleState(cycle_id="01", requested_cycle_count=None)  # Explicitly None
        mock_jules_client.run_session.return_value = JulesSessionResult(
            status="success", session_name="test", pr_url="http://pr"
        )

        # Act
        await nodes.architect_session_node(state.model_dump())

        # Assert
        mock_jules_client.run_session.assert_awaited_once()
        call_kwargs = mock_jules_client.run_session.await_args.kwargs
        assert "max_cycles" not in call_kwargs["user_prompt"]

    @pytest.mark.parametrize("count", [1, 2, 3, 5, 10])
    async def test_prompt_injection_various_counts(
        self,
        count: int,
        mock_sandbox: MagicMock,
        mock_jules_client: AsyncMock,
        mock_session_manager: AsyncMock,
        mock_git_manager: AsyncMock,
        mock_console: MagicMock,
    ) -> None:
        """Test with various valid cycle counts."""
        # Arrange
        nodes = CycleNodes(mock_sandbox, mock_jules_client, mock_session_manager, mock_git_manager)
        state = CycleState(cycle_id="01", requested_cycle_count=count)
        mock_jules_client.run_session.return_value = JulesSessionResult(
            status="success", session_name="test", pr_url="http://pr"
        )

        # Act
        await nodes.architect_session_node(state.model_dump())

        # Assert
        mock_jules_client.run_session.assert_awaited_once()
        call_kwargs = mock_jules_client.run_session.await_args.kwargs
        assert f"max_cycles={count}" in call_kwargs["user_prompt"]
