from unittest.mock import AsyncMock, patch

import pytest
from ac_cdd_core.graph import GraphBuilder
from ac_cdd_core.state import CycleState


@pytest.mark.asyncio
async def test_init_branch_node(mock_services):
    """Test init_branch node."""
    builder = GraphBuilder(mock_services)
    # Mock git manager
    builder.git = AsyncMock()
    builder.git.create_working_branch.return_value = "feature/cycle01"

    state = CycleState(cycle_id="01")
    result = await builder.init_branch_node(state)
    assert result["current_phase"] == "branch_ready"
    assert result["active_branch"] == "feature/cycle01"


@pytest.mark.asyncio
async def test_architect_session_failed(mock_services):
    """Test architect session failure when SPEC is missing."""
    builder = GraphBuilder(mock_services)

    with patch("pathlib.Path.exists", return_value=False):
        state = CycleState(cycle_id="01")
        result = await builder.architect_session_node(state)
        assert result["current_phase"] == "architect_failed"
        assert "ALL_SPEC.md not found" in result["error"]


@pytest.mark.asyncio
async def test_architect_graph_structure(mock_services):
    """Test architect graph structure and compilation."""
    builder = GraphBuilder(mock_services)
    graph = builder.build_architect_graph()

    # Verify nodes
    assert "init_branch" in graph.nodes
    assert "architect_session" in graph.nodes
    assert "commit" in graph.nodes


@pytest.mark.asyncio
async def test_coder_graph_structure(mock_services):
    """Test coder graph structure and compilation."""
    builder = GraphBuilder(mock_services)
    graph = builder.build_coder_graph()

    # Verify nodes
    assert "checkout_branch" in graph.nodes
    assert "coder_session" in graph.nodes
    assert "run_tests" in graph.nodes
    assert "uat_evaluate" in graph.nodes
    assert "auditor" in graph.nodes
    assert "commit" in graph.nodes
