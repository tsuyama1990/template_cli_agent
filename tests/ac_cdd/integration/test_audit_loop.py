from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.domain_models import JulesSessionResult
from ac_cdd_core.graph import GraphBuilder
from ac_cdd_core.service_container import ServiceContainer
from ac_cdd_core.services.artifacts import ArtifactManager
from ac_cdd_core.services.contracts import ContractManager
from ac_cdd_core.services.file_ops import FilePatcher


@pytest.fixture
def mock_services(
    mock_jules_client: AsyncMock,
    mock_file_patcher: MagicMock,
    mock_contract_manager: MagicMock,
    mock_artifact_manager: MagicMock,
) -> ServiceContainer:
    """Provides a mock ServiceContainer."""
    return ServiceContainer(
        jules=mock_jules_client,
        file_patcher=mock_file_patcher,
        contract_manager=mock_contract_manager,
        artifact_manager=mock_artifact_manager,
    )


@pytest.mark.asyncio
@patch("ac_cdd_core.graph.SandboxRunner", return_value=AsyncMock())
@patch("ac_cdd_core.graph.GitManager", return_value=AsyncMock())
@patch("ac_cdd_core.graph.SessionManager", return_value=AsyncMock())
@patch("ac_cdd_core.graph.JulesClient")
async def test_architect_session_node_creates_pr(
    mock_jules_client_cls: MagicMock,
    _mock_session_manager_cls: MagicMock,
    _mock_git_manager_cls: MagicMock,
    _mock_sandbox_cls: MagicMock,
) -> None:
    """
    Tests that the architect_session_node correctly calls the Jules client
    and returns the expected state upon successful PR creation.
    """
    # Arrange: Configure a mock instance for JulesClient
    mock_jules_instance = AsyncMock()
    mock_jules_instance.run_session = AsyncMock(
        return_value=JulesSessionResult(
            status="success",
            pr_url="https://github.com/test/repo/pull/1",
            session_name="sessions/architect-123",
        )
    )
    mock_jules_client_cls.return_value = mock_jules_instance

    # Arrange: Create a ServiceContainer with the mocked client
    services = ServiceContainer(
        jules=mock_jules_instance,
        file_patcher=MagicMock(spec=FilePatcher),
        contract_manager=MagicMock(spec=ContractManager),
        artifact_manager=MagicMock(spec=ArtifactManager),
    )

    # Arrange: Instantiate the GraphBuilder. It will use the patched dependencies.
    builder = GraphBuilder(services)
    state = {"cycle_id": "01", "project_session_id": "test-session"}

    # Act
    result = await builder.nodes.architect_session_node(state)

    # Assert
    assert result["status"] == "architect_completed"
    assert result["pr_url"] == "https://github.com/test/repo/pull/1"
