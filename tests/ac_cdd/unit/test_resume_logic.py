from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.domain_models import CycleManifest, JulesSessionResult
from ac_cdd_core.graph_nodes import CycleNodes


@pytest.mark.asyncio
class TestResumeLogic:
    @pytest.fixture
    def mock_dependencies(self) -> tuple[MagicMock, AsyncMock]:
        """Provides mocked SandboxRunner and JulesClient."""
        sandbox = MagicMock()
        jules = AsyncMock()
        return sandbox, jules

    async def test_resume_from_existing_session_and_persist(
        self, mock_dependencies: tuple[MagicMock, AsyncMock]
    ) -> None:
        """Test that an existing session ID is used and the result is persisted."""
        # Arrange
        sandbox, jules = mock_dependencies
        mock_session_manager = AsyncMock()
        mock_git_manager = AsyncMock()

        # Arrange: Setup manifest with an existing session ID
        cycle = CycleManifest(id="01", jules_session_id="jules-existing-123")
        mock_session_manager.get_cycle.return_value = cycle
        jules.continue_session.return_value = JulesSessionResult(
            status="success",
            session_name="jules-existing-123",
            pr_url="http://pr-existing",
        )

        nodes = CycleNodes(sandbox, jules, mock_session_manager, mock_git_manager)
        state = {"cycle_id": "01", "iteration_count": 1, "resume_mode": True}

        # Act
        result = await nodes.coder_session_node(state)

        # Assert
        jules.continue_session.assert_awaited_once_with("jules-existing-123")
        mock_session_manager.update_cycle_state.assert_awaited_with(
            "01", pr_url="http://pr-existing", status="completed"
        )
        assert result["pr_url"] == "http://pr-existing"
        assert result["status"] == "ready_for_audit"

    async def test_fallback_to_new_session_and_persist(
        self, mock_dependencies: tuple[MagicMock, AsyncMock]
    ) -> None:
        """Test that if no session exists, a new one is started and its ID is persisted."""
        # Arrange
        sandbox, jules = mock_dependencies
        mock_session_manager = AsyncMock()
        mock_git_manager = AsyncMock()

        # Arrange: Setup manifest with NO existing session
        cycle = CycleManifest(id="01", jules_session_id=None)
        mock_session_manager.get_cycle.return_value = cycle

        # Arrange: Setup Jules to return a new session that requires waiting
        jules.run_session.return_value = JulesSessionResult(
            session_name="jules-new-456", status="running", pr_url=None
        )
        jules.wait_for_completion.return_value = JulesSessionResult(
            status="success",
            session_name="jules-new-456",
            pr_url="http://pr-new",
        )

        nodes = CycleNodes(sandbox, jules, mock_session_manager, mock_git_manager)
        state = {"cycle_id": "01", "iteration_count": 1, "resume_mode": True}

        # Act
        result = await nodes.coder_session_node(state)

        # Assert
        jules.run_session.assert_awaited_once()
        # The node should request the session ID early for persistence
        assert jules.run_session.await_args.kwargs["require_plan_approval"] is True

        # Verify the new session ID was immediately persisted
        mock_session_manager.update_cycle_state.assert_any_await(
            "01", jules_session_id="jules-new-456", status="in_progress"
        )
        # Verify the node then waits for the new session to complete
        jules.wait_for_completion.assert_awaited_once_with("jules-new-456")

        assert result["status"] == "ready_for_audit"
        assert result["pr_url"] == "http://pr-new"
        # Verify the final state is persisted
        mock_session_manager.update_cycle_state.assert_awaited_with(
            "01", pr_url="http://pr-new", status="completed"
        )
