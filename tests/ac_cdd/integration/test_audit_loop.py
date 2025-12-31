from unittest.mock import AsyncMock, MagicMock

import pytest
from ac_cdd_core.graph import GraphBuilder
from ac_cdd_core.state import CycleState


@pytest.mark.asyncio
async def test_audit_rejection_loop():
    """
    Test that the audit loop functions correctly when changes are requested.
    Verifies that the graph iterates through 3 auditors * 2 reviews each = 6 cycles.
    """
    # Mock Services
    mock_services = MagicMock()
    mock_services.git = AsyncMock()
    mock_services.jules = MagicMock()
    mock_services.sandbox = MagicMock()
    mock_services.reviewer = MagicMock()

    # Mock Jules run session
    mock_services.jules.run_session = AsyncMock(return_value={"pr_url": "http://pr"})
    mock_services.jules.continue_session = AsyncMock(
        return_value={"pr_url": "http://pr"}
    )

    # Mock Sandbox
    mock_services.sandbox.run_lint_check = AsyncMock(return_value=(True, "OK"))

    # Mock Reviewer
    mock_services.reviewer.review_code = AsyncMock(
        return_value="CHANGES_REQUESTED: Please fix X."
    )

    # Build Graph
    builder = GraphBuilder(mock_services)

    # Mock implementation nodes to isolate loop logic
    # We patch the instance methods BEFORE building the graph
    builder.checkout_branch_node = AsyncMock(return_value={"current_phase": "coding"})

    # Coder Node needs to return correct state keys to satisfy future nodes AND increment counters
    def coder_side_effect(state: CycleState):
        # Determine if we are in initial run or fix loop
        # Initial run: review_count=0 (default) or 1?
        # State defaults: review_count=1, auditor_idx=1

        # If we are called, it means we are either implementing (Iter 1) or Fixing.
        # If it's a FIX, we must increment counters for the NEXT review.
        # But wait, initially review_count=1. Auditor reviews. If fail -> Coder(Fix).
        # Coder(Fix) must increment so next Auditor call sees review_count=2.

        # But for Iter 1, we shouldn't increment?
        # Graph flow: default state (rev=1) -> Coder -> Auditor (reads rev 1)
        # -> Loops back if needed.
        # If loops back, Coder must increment?
        # Let's see current logic in graph.py:
        # "next_review_count = state.current_auditor_review_count + 1"
        # "if next_review_count > REVIEWS_PER_AUDITOR: ..."

        # So YES, the Coder node is responsible for advancing the state.

        # Logic:
        # If it's the very first run (iteration=0 or 1), maybe we don't increment?
        # Actually logic is: Coder runs -> then Auditor runs.
        # Auditor checks "Review {count}".
        # So for Review 1, state must be 1.
        # After Auditor 1 rejects, we go back to Coder. Coder fixes.
        # NOW request is for Review 2. So Coder must increment state to 2.

        # BUT... the first time Coder runs (State=1). It returns. Auditor runs (State=1).
        # If Mock increments blindly, first run becomes State=2? No.
        # Detailed logic:
        # We need to know if this is a "fix" run.
        # Since we are mocking, let's keep it simple: Just always increment?
        # No, initial run is special.

        # Let's inspect state.
        # If we just depend on call count?

        current_rev = state.current_auditor_review_count
        current_aud = state.current_auditor_index

        # Current logic in graph.py increments *after* fix.
        # So input state has old count. Coder returns NEW count.

        # Example trace:
        # Start: Rev=1, Aud=1.
        # 1. Coder (Input Rev=1). Returns ? Should return Rev=1?
        #    Auditor runs (Rev=1).
        #    Auditor returns "coder_session" (loop).
        # 2. Coder (Input Rev=1). This is the FIX. Returns Rev=2.
        #    Auditor runs (Rev=2).
        #    Auditor returns "coder_session" (next auditor).

        # So we increment ONLY if we are fixing?
        # How does real code know? "if iteration_count > 1 and state.jules_session_name:"
        # We can use iteration_count from state.

        # BUT for this test we want to force the loop forward.
        # Let's just implement a simple state advancer.

        # Logic matches graph.py:86: "next_review_count = state.current_auditor_review_count + 1"
        # But only if it's a fix?
        # Creating a helper to simulate the logic.

        # settings.REVIEWS_PER_AUDITOR = 2 (default)
        REVIEWS_PER_AUDITOR = 2

        updated_rev = current_rev
        updated_aud = current_aud

        if (
            state.iteration_count > 0
        ):  # graph increments iteration count in coder node too
            # It's a fix or next cycle
            updated_rev += 1
            if updated_rev > REVIEWS_PER_AUDITOR:
                updated_rev = 1
                updated_aud += 1

        return {
            "current_phase": "coded",
            "pr_url": "http://pr",
            "current_auditor_index": updated_aud,
            "current_auditor_review_count": updated_rev,
            "iteration_count": state.iteration_count + 1,
        }

    builder.coder_session_node = AsyncMock(side_effect=coder_side_effect)
    builder.syntax_check_node = AsyncMock(return_value={"current_phase": "syntax_ok"})

    # We use REAL auditor_node logic to test the loop decision
    # But current auditor_node calls self.git.get_changed_files().
    # We must mock that service response.
    mock_services.git.get_changed_files = AsyncMock(return_value=["file.py"])
    # And file reads? auditor_node calls Path.read_text.
    # We should mock that or ensuring the node doesn't crash.
    # The file read catch block logs check error, so it should proceed even if files missing.
    # The reviewer call uses the content.

    graph = builder.build_coder_graph()

    # Initialize State
    initial_state = CycleState(
        cycle_id="01",
        session_id="test_session",
        integration_branch="dev/main",
        cycle_branch="dev/cycle01",
    )

    # Run Graph
    final_state = await graph.ainvoke(initial_state)

    # Verification

    # 1. Total Auditor Calls: 3 auditors * 2 reviews = 6 calls
    assert mock_services.reviewer.review_code.call_count == 6

    # 2. Total Jules Calls
    # Note: Since we mocked `coder_session_node`, we can't assert on `jules_client` calls inside it!
    # We must assert on `coder_session_node` calls instead.
    # Initial Run (1) + Fixes (0) ?
    # Wait, if we mock `coder_session_node`, we must effectively simulate "fix" logic
    # if we want to confirm params.
    # But here we just want to verify the LOOP invokes the coder node 6 times + 1 initial time?
    # No, iteration 1 runs coder. Auditor triggers re-run.
    # Total coder node calls:
    # 1 (Iter 1) -> Audit 1.1 -> Loop
    # 2 (Iter 2) -> Audit 1.2 -> Loop
    # ...
    # 6 (Iter 6) -> Audit 3.2 -> Commit
    # So Coder Node should be called 6 times.
    assert builder.coder_session_node.call_count == 6

    # 3. Final State
    assert "completion_message" not in final_state

    # Check if we reached commit
    # Since we didn't mock commit_coder_node, it might fail?
    # commit_coder_node calls self.git.merge_to_integration.
    assert mock_services.git.merge_to_integration.called


@pytest.mark.asyncio
async def test_audit_loop_early_exit_scenario():
    """Hypothetical future test."""
    pass
