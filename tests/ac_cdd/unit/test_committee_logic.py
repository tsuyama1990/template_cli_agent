from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.domain_models import AuditResult, CommitteeDecision
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


def test_committee_logic_flow(cycle_nodes: CycleNodes) -> None:
    # Test case 1: Approved by the first auditor, moves to the next
    state_approved_first = {
        "cycle_id": "01",
        "iteration_count": 1,
        "num_auditors": 3,
        "current_auditor_index": 0,
        "audit_result": AuditResult(
            status="approved", feedback="LGTM", is_approved=True
        ).model_dump(),
    }
    result = cycle_nodes.committee_manager_node(state_approved_first)
    decision = CommitteeDecision(**result["committee_decision"])
    assert decision.status == "next_auditor"
    assert result["current_auditor_index"] == 1

    # Test case 2: Approved by the final auditor, cycle is approved
    state_approved_final = {
        "cycle_id": "01",
        "iteration_count": 1,
        "num_auditors": 3,
        "current_auditor_index": 2,
        "audit_result": AuditResult(
            status="approved", feedback="All good", is_approved=True
        ).model_dump(),
    }
    result = cycle_nodes.committee_manager_node(state_approved_final)
    decision = CommitteeDecision(**result["committee_decision"])
    assert decision.status == "cycle_approved"

    # Test case 3: Rejected by an auditor, requires retry
    state_rejected = {
        "cycle_id": "01",
        "iteration_count": 1,
        "num_auditors": 3,
        "current_auditor_index": 1,
        "audit_result": AuditResult(
            status="rejected", feedback="Needs changes", is_approved=False
        ).model_dump(),
    }
    result = cycle_nodes.committee_manager_node(state_rejected)
    decision = CommitteeDecision(**result["committee_decision"])
    assert decision.status == "retry_fix"
    assert result["iteration_count"] == 2  # Iteration is incremented

    # Test case 4: Audit process failed
    state_failed = {
        "cycle_id": "01",
        "iteration_count": 1,
        "num_auditors": 3,
        "current_auditor_index": 0,
        "audit_result": AuditResult(
            status="failed", feedback="Crashed", is_approved=False
        ).model_dump(),
    }
    result = cycle_nodes.committee_manager_node(state_failed)
    decision = CommitteeDecision(**result["committee_decision"])
    assert decision.status == "cycle_failed"


def test_committee_pipeline_handover(cycle_nodes: CycleNodes) -> None:
    # Test routing for 'next_auditor'
    state_next = {
        "committee_decision": CommitteeDecision(
            status="next_auditor", reason="Approved"
        ).model_dump()
    }
    assert cycle_nodes.route_committee(state_next) == "auditor"

    # Test routing for 'retry_fix'
    state_retry = {
        "committee_decision": CommitteeDecision(status="retry_fix", reason="Rejected").model_dump()
    }
    assert cycle_nodes.route_committee(state_retry) == "coder_session"

    # Test routing for 'cycle_approved'
    state_approved = {
        "committee_decision": CommitteeDecision(
            status="cycle_approved", reason="All auditors approved"
        ).model_dump()
    }
    assert cycle_nodes.route_committee(state_approved) == "uat_evaluate"

    # Test routing for 'cycle_failed'
    state_failed = {
        "committee_decision": CommitteeDecision(status="cycle_failed", reason="Crashed").model_dump()
    }
    assert cycle_nodes.route_committee(state_failed) == "failed"
