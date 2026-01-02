from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.graph import GraphBuilder
from ac_cdd_core.state import CycleState
from ac_cdd_core.domain_models import AuditResult
from ac_cdd_core.service_container import ServiceContainer


@pytest.mark.asyncio
async def test_audit_rejection_loop():
    """
    Test that the audit loop functions correctly when changes are requested.
    Verifies that the graph iterates through retries.
    """
    # Mock Services
    mock_services = MagicMock(spec=ServiceContainer)
    mock_services.git = AsyncMock()
    mock_services.jules = MagicMock()
    mock_services.sandbox = MagicMock()
    mock_services.reviewer = MagicMock()

    # Mock Jules run session
    mock_services.jules.run_session = AsyncMock(return_value={"pr_url": "http://pr", "status": "success"})

    # Mock Reviewer
    mock_services.reviewer.review_code = AsyncMock(return_value="CHANGES_REQUESTED: Please fix X.")

    # Build Graph
    # Important: GraphBuilder uses DI from ServiceContainer in refactored code.
    # We pass ServiceContainer instance directly.
    # However, GraphBuilder expects services.jules and services.sandbox.

    # Also, we need to mock PlanAuditor inside AuditOrchestrator because it uses PydanticAI
    # which fails if API key is missing (even in test if validation runs).
    # We patch `ac_cdd_core.services.audit_orchestrator.PlanAuditor` to avoid PydanticAI instantiation.

    with patch("ac_cdd_core.services.audit_orchestrator.PlanAuditor") as MockAuditorClass:
        mock_plan_auditor = MockAuditorClass.return_value
        # If AuditOrchestrator is used, it might call `auditor.audit_plan`.
        # But here we are testing Coder Graph which uses LLMReviewer (auditor_node).
        # CycleNodes init: `self.audit_orchestrator = AuditOrchestrator(...)`.
        # This triggers PlanAuditor init. So patching class works.

        builder = GraphBuilder(mock_services)
        graph = builder.build_coder_graph()

        initial_state = CycleState(
            cycle_id="01",
            session_id="test_session",
            integration_branch="dev/main",
        )

        with patch("ac_cdd_core.config.settings.max_audit_retries", 2):

            async def mock_coder_node(state: CycleState):
                return {
                    "status": "ready_for_audit",
                    "pr_url": "http://pr",
                    "iteration_count": state.iteration_count + 1
                }

            with patch.object(builder.nodes, "coder_session_node", side_effect=mock_coder_node) as mock_coder:

                # Fix for Checkpointer error: pass a thread_id configuration
                config = {"configurable": {"thread_id": "test_audit_loop"}}

                # Run the graph.
                # Recursion limit increased to handle the loop steps.
                # If it fails with GraphRecursionError, it means the termination condition (count >= max) wasn't met.

                final_state = await graph.ainvoke(initial_state, {**config, "recursion_limit": 25})

                assert mock_coder.call_count == 2
                assert mock_services.reviewer.review_code.call_count == 2

                # Verify final state has rejected status (since we forced loop to exhaustion)
                # AuditResult is nested in state
                assert final_state["audit_result"].status == "REJECTED"
                assert final_state["audit_result"].is_approved is False
