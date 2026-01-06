from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.domain_models import JulesSessionResult
from ac_cdd_core.graph import GraphBuilder
from ac_cdd_core.service_container import ServiceContainer
from ac_cdd_core.services.artifacts import ArtifactManager
from ac_cdd_core.services.contracts import ContractManager
from ac_cdd_core.services.file_ops import FilePatcher
from ac_cdd_core.state import CycleState


@pytest.fixture
def mock_sandbox_runner() -> AsyncMock:
    return AsyncMock()


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
def mock_file_patcher() -> MagicMock:
    return MagicMock(spec=FilePatcher)


@pytest.fixture
def mock_contract_manager() -> MagicMock:
    return MagicMock(spec=ContractManager)


@pytest.fixture
def mock_artifact_manager() -> MagicMock:
    return MagicMock(spec=ArtifactManager)


@pytest.fixture
def mock_services(
    mock_jules_client: AsyncMock,
    mock_file_patcher: MagicMock,
    mock_contract_manager: MagicMock,
    mock_artifact_manager: MagicMock,
) -> ServiceContainer:
    return ServiceContainer(
        jules=mock_jules_client,
        file_patcher=mock_file_patcher,
        contract_manager=mock_contract_manager,
        artifact_manager=mock_artifact_manager,
    )


@pytest.mark.asyncio
@patch("ac_cdd_core.graph.SandboxRunner")
@patch("ac_cdd_core.graph.GitManager")
@patch("ac_cdd_core.graph.SessionManager")
@patch("ac_cdd_core.graph.JulesClient")
async def test_graph_builder_initialization(
    mock_jules_cls: MagicMock,
    mock_session_manager_cls: MagicMock,
    mock_git_manager_cls: MagicMock,
    mock_sandbox_cls: MagicMock,
) -> None:
    """Tests that GraphBuilder initializes its components correctly."""
    mock_jules_instance = mock_jules_cls.return_value
    mock_session_manager_instance = mock_session_manager_cls.return_value
    mock_git_manager_instance = mock_git_manager_cls.return_value
    # Use AsyncMock for the instance so its methods are awaitable
    mock_sandbox_instance = AsyncMock()
    mock_sandbox_cls.return_value = mock_sandbox_instance

    services = ServiceContainer(
        jules=mock_jules_instance,
        file_patcher=MagicMock(spec=FilePatcher),
        contract_manager=MagicMock(spec=ContractManager),
        artifact_manager=MagicMock(spec=ArtifactManager),
    )
    builder = GraphBuilder(services)

    assert builder.sandbox == mock_sandbox_instance
    assert builder.jules == mock_jules_instance
    assert builder.session_manager == mock_session_manager_instance
    assert builder.git_manager == mock_git_manager_instance

    # Test cleanup
    await builder.cleanup()
    mock_sandbox_instance.cleanup.assert_awaited_once()


@pytest.mark.asyncio
async def test_architect_graph_execution(
    mock_services: ServiceContainer, mock_jules_client: AsyncMock
) -> None:
    """Tests that the architect graph executes the main node."""
    # Arrange
    mock_jules_client.run_session.return_value = JulesSessionResult(
        status="success",
        session_name="test-session",
        pr_url="http://github.com/pr/1",
    )
    initial_state = CycleState(cycle_id="01")
    config = {"configurable": {"thread_id": "test-thread"}}

    with patch("ac_cdd_core.graph.CycleNodes", spec=True) as mock_cycle_nodes_cls:
        mock_nodes_instance = mock_cycle_nodes_cls.return_value
        mock_nodes_instance.architect_session_node.return_value = {
            "status": "architect_completed",
            "pr_url": "http://github.com/pr/1",
        }

        builder = GraphBuilder(mock_services)
        graph = builder.build_architect_graph()

        # Act
        result_state = await graph.ainvoke(initial_state.model_dump(), config=config)

        # Assert
        mock_nodes_instance.architect_session_node.assert_awaited_once()
        # Verify that the node was called with the correct state object
        call_args, _ = mock_nodes_instance.architect_session_node.await_args
        assert isinstance(call_args[0], CycleState)
        assert call_args[0].cycle_id == "01"
        assert result_state["status"] == "architect_completed"


@pytest.mark.asyncio
async def test_coder_graph_full_flow(
    mock_services: ServiceContainer, mock_jules_client: AsyncMock
) -> None:
    """Tests a simplified happy path of the coder graph."""
    initial_state = CycleState(cycle_id="01", iteration_count=1)
    config = {"configurable": {"thread_id": "test-thread"}}

    with patch("ac_cdd_core.graph.CycleNodes", spec=True) as mock_cycle_nodes_cls:
        mock_nodes = mock_cycle_nodes_cls.return_value
        # Mock node behaviors
        mock_nodes.coder_session_node.return_value = {"status": "ready_for_audit"}
        mock_nodes.check_coder_outcome.return_value = "ready_for_audit"
        mock_nodes.auditor_node.return_value = {"status": "approved"}
        mock_nodes.committee_manager_node.return_value = {"status": "cycle_approved"}
        mock_nodes.route_committee.return_value = "uat_evaluate"
        mock_nodes.uat_evaluate_node.return_value = {"status": "cycle_completed"}

        builder = GraphBuilder(mock_services)
        graph = builder.build_coder_graph()

        # Act
        result_state = await graph.ainvoke(initial_state.model_dump(), config=config)

        # Assert
        mock_nodes.coder_session_node.assert_awaited_once()
        mock_nodes.auditor_node.assert_awaited_once()
        mock_nodes.committee_manager_node.assert_called_once()
        mock_nodes.uat_evaluate_node.assert_awaited_once()
        assert result_state["status"] == "cycle_completed"


@pytest.mark.asyncio
async def test_coder_graph_retry_loop(
    mock_services: ServiceContainer, mock_jules_client: AsyncMock
) -> None:
    """Tests the retry loop between coder and auditor."""
    initial_state = CycleState(cycle_id="01", iteration_count=1)
    config = {"configurable": {"thread_id": "test-thread"}}

    with patch("ac_cdd_core.graph.CycleNodes", spec=True) as mock_cycle_nodes_cls:
        mock_nodes = mock_cycle_nodes_cls.return_value
        # Mock node behaviors to simulate one rejection, then approval
        mock_nodes.coder_session_node.side_effect = [
            {"status": "ready_for_audit"},  # First call
            {"status": "ready_for_audit"},  # Second call after retry
        ]
        mock_nodes.check_coder_outcome.return_value = "ready_for_audit"
        mock_nodes.auditor_node.side_effect = [
            {"status": "rejected"},
            {"status": "approved"},
        ]
        mock_nodes.committee_manager_node.side_effect = [
            {"status": "retry_fix"},
            {"status": "cycle_approved"},
        ]
        mock_nodes.route_committee.side_effect = ["coder_session", "uat_evaluate"]
        mock_nodes.uat_evaluate_node.return_value = {"status": "cycle_completed"}

        builder = GraphBuilder(mock_services)
        graph = builder.build_coder_graph()

        # Act
        await graph.ainvoke(initial_state.model_dump(), config=config)

        # Assert
        assert mock_nodes.coder_session_node.await_count == 2
        assert mock_nodes.auditor_node.await_count == 2
        assert mock_nodes.committee_manager_node.call_count == 2
        mock_nodes.uat_evaluate_node.assert_awaited_once()
