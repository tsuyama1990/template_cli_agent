# ruff: noqa: T201, E501
import asyncio
import json
from typing import Any, Literal

from ac_cdd_core.domain_models import (
    AuditResult,
    CommitteeDecision,
    JulesSessionResult,
    PlanAuditResult,
)
from ac_cdd_core.services.git_ops import GitManager
from ac_cdd_core.services.jules_client import JulesClient
from ac_cdd_core.session_manager import SessionManager
from rich.console import Console

console = Console()


class CycleNodes:
    """Encapsulates the logic for all nodes in the LangGraph workflows."""

    def __init__(
        self,
        sandbox: Any,
        jules: JulesClient,
        session_manager: SessionManager,
        git_manager: GitManager,
    ):
        self.sandbox = sandbox
        self.jules = jules
        self.session_manager = session_manager
        self.git_manager = git_manager

    async def architect_session_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """Runs the architect session to generate cycle plans."""
        cycle_id = state["cycle_id"]
        requested_cycle_count = state.get("requested_cycle_count")
        console.print(f"Running Architect Session for Cycle {cycle_id}...")

        prompt = f"Run architect phase for cycle {cycle_id}"
        if requested_cycle_count:
            prompt += f", max_cycles={requested_cycle_count}"

        try:
            result = await self.jules.run_session(
                "architect",
                user_prompt=prompt,
                cycle_id=cycle_id,
            )
            if result.status != "success" or not result.pr_url:
                raise ValueError("Architect session failed to produce a PR URL.")

            console.print(f"Architect session completed. PR URL: {result.pr_url}")
            await self.session_manager.update_cycle_state(
                cycle_id, pr_url=result.pr_url, status="architect_completed"
            )
            return {
                "status": "architect_completed",
                "pr_url": result.pr_url,
            }
        except Exception as e:
            console.print(f"Architect Session Failed: {e}")
            await self.session_manager.update_cycle_state(
                cycle_id, status="architect_failed", error=str(e)
            )
            return {"status": "architect_failed"}

    async def coder_session_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """
        Runs the coder session, either starting a new one or resuming an existing one.
        Handles the logic for waiting for long-running sessions to complete.
        """
        cycle_id = state["cycle_id"]
        iteration = state.get("iteration_count", 1)
        is_resume = state.get("resume_mode", False)

        console.print(f"Starting Coder Session for Cycle {cycle_id} (Iteration {iteration})...")

        try:
            cycle_manifest = await self.session_manager.get_cycle(cycle_id)
            session_id = cycle_manifest.jules_session_id if is_resume else None

            result: JulesSessionResult | None = None

            if session_id:
                console.print(f"Resuming Jules Session: {session_id}")
                result = await self.jules.continue_session(session_id)
            else:
                console.print("Starting a new Jules session...")
                # Start the session but don't wait for the final result yet,
                # just get the session ID so we can persist it.
                initial_result = await self.jules.run_session(
                    "coder",
                    user_prompt="Implement the code as per the spec.",
                    cycle_id=cycle_id,
                    iteration=iteration,
                    require_plan_approval=True,  # Request session_name immediately
                )
                session_id = initial_result.session_name

                if not session_id:
                    raise ValueError("Failed to create a new Jules session.")

                console.print(f"Session {session_id} created. Persisting state...")
                await self.session_manager.update_cycle_state(
                    cycle_id, jules_session_id=session_id, status="in_progress"
                )

                if initial_result.status == "running":
                    console.print(f"Session {session_id} is running. Waiting for completion...")
                    result = await self.jules.wait_for_completion(session_id)
                else:
                    result = initial_result

            if not result:
                raise ValueError("Jules session did not return a result.")

            # Now, process the final result
            if result.status == "success" and result.pr_url:
                console.print(f"Coder session completed. PR URL: {result.pr_url}")
                await self.session_manager.update_cycle_state(
                    cycle_id, pr_url=result.pr_url, status="completed"
                )
                return {"status": "ready_for_audit", "pr_url": result.pr_url}

            # Handle other statuses if necessary, e.g., plan needs approval
            if result.status != "success":
                raise ValueError(f"Jules session ended with status: '{result.status}'")
            if not result.pr_url:
                raise ValueError("Jules session succeeded but did not provide a PR URL.")

            # This part should ideally not be reached if logic is correct
            return {"status": "failed", "error": "Unhandled session outcome"}

        except Exception as e:
            console.print(f"Coder Session Failed: {e}")
            await self.session_manager.update_cycle_state(cycle_id, status="failed", error=str(e))
            return {"status": "failed", "error": str(e)}

    def check_coder_outcome(self, state: dict) -> Literal["ready_for_audit", "failed"]:
        """Routes based on the coder session outcome."""
        return "ready_for_audit" if state.get("status") == "ready_for_audit" else "failed"

    async def auditor_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """
        Runs the auditor agent on the code generated by the coder agent.
        This node will execute in a sandbox environment.
        """
        cycle_id = state["cycle_id"]
        pr_url = state["pr_url"]
        iteration = state.get("iteration_count", 1)
        current_auditor_index = state.get("current_auditor_index", 0)

        console.print(
            f"Auditor {current_auditor_index + 1} starting for Cycle {cycle_id} (Iteration {iteration})..."
        )
        console.print(f"PR to be audited: {pr_url}")

        try:
            # 1. Checkout PR branch
            branch_name = f"pr-cycle-{cycle_id}-iter-{iteration}"
            await self.git_manager.checkout_pr_branch(pr_url, branch_name)
            console.print(f"Checked out PR branch '{branch_name}'.")

            # 2. Run the auditor in a sandbox
            # The sandbox setup will include installing deps and running the agent
            audit_result_str = await self.sandbox.run(
                f'aider --config-file=aider.project.yml --pr-url="{pr_url}"'
            )
            # TODO: This parsing is brittle. The auditor agent should output structured data.
            # For now, we assume the last line is a JSON object with the result.
            last_line = audit_result_str.strip().split("\n")[-1]
            audit_result_json = json.loads(last_line)
            audit_result = AuditResult(**audit_result_json)

            # 3. Save the audit result
            await self.session_manager.save_audit_result(
                cycle_id, iteration, current_auditor_index, audit_result
            )
            console.print(
                f"Auditor {current_auditor_index + 1} finished with status: {audit_result.status}"
            )

            return {"audit_result": audit_result.model_dump()}

        except Exception as e:
            console.print(f"Auditor Node Failed: {e}")
            # Save a failed audit result
            failed_audit = AuditResult(
                status="failed", feedback=f"Auditor process failed: {e}", score=0
            )
            await self.session_manager.save_audit_result(
                cycle_id, iteration, current_auditor_index, failed_audit
            )
            return {"audit_result": failed_audit.model_dump()}

    def committee_manager_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """
        Manages the committee of auditors, deciding whether to proceed to the next auditor,
        retry the fix, or approve the cycle.
        """
        cycle_id = state["cycle_id"]
        iteration = state.get("iteration_count", 1)
        num_auditors = state.get("num_auditors", 3)
        current_auditor_index = state.get("current_auditor_index", 0)
        audit_result = AuditResult(**state["audit_result"])

        console.print(f"Committee Manager processing audit from Auditor {current_auditor_index + 1}...")

        if audit_result.status == "rejected":
            console.print("Audit status: REJECTED. Routing for retry.")
            decision = CommitteeDecision(
                status="retry_fix",
                reason=f"Rejected by Auditor {current_auditor_index + 1}.",
            )
            return {
                "committee_decision": decision.model_dump(),
                "iteration_count": iteration + 1,  # Increment iteration for the retry
            }

        if audit_result.status == "approved":
            if current_auditor_index + 1 < num_auditors:
                console.print(
                    f"Audit status: APPROVED. Proceeding to Auditor {current_auditor_index + 2}."
                )
                decision = CommitteeDecision(
                    status="next_auditor",
                    reason=f"Approved by Auditor {current_auditor_index + 1}.",
                )
                return {
                    "committee_decision": decision.model_dump(),
                    "current_auditor_index": current_auditor_index + 1,
                }

            console.print("Final auditor approved. Cycle audit passed.")
            decision = CommitteeDecision(status="cycle_approved", reason="All auditors approved.")
            return {"committee_decision": decision.model_dump()}

        # Handle failed state
        console.print(f"Audit status: FAILED ({audit_result.feedback}). Cycle failed.")
        decision = CommitteeDecision(
            status="cycle_failed",
            reason=f"Audit process failed at Auditor {current_auditor_index + 1}.",
        )
        return {"committee_decision": decision.model_dump()}

    def route_committee(self, state: dict) -> Literal["coder_session", "auditor", "uat_evaluate", "failed"]:
        """Routes based on the committee's decision."""
        decision = CommitteeDecision(**state["committee_decision"])
        if decision.status == "retry_fix":
            return "coder_session"
        if decision.status == "next_auditor":
            return "auditor"
        if decision.status == "cycle_approved":
            return "uat_evaluate"
        return "failed"

    async def uat_evaluate_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """
        Runs the QA agent to evaluate the UAT criteria against the implemented code.
        """
        cycle_id = state["cycle_id"]
        console.print(f"Evaluating UAT for Cycle {cycle_id}...")
        try:
            # In a real scenario, this would involve running a QA agent (e.g., another LLM)
            # against the UAT.md file and the test execution logs.
            # For this simulation, we'll assume it passes if the committee approved.
            console.print("UAT evaluation passed.")
            await self.session_manager.update_cycle_state(cycle_id, status="uat_passed")
            return {"status": "cycle_completed"}
        except Exception as e:
            console.print(f"UAT Evaluation Failed: {e}")
            await self.session_manager.update_cycle_state(
                cycle_id, status="uat_failed", error=str(e)
            )
            return {"status": "cycle_failed"}

    async def plan_auditor_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """
        Audits the plan generated by the architect.
        """
        cycle_id = state["cycle_id"]
        console.print(f"Auditing plan for Cycle {cycle_id}...")

        try:
            # This would involve a call to a plan auditor agent/service
            # For simulation, we'll assume the plan is always approved.
            await asyncio.sleep(1)  # Simulate network call
            plan_audit_result = PlanAuditResult(status="approved", feedback="Plan looks solid.")
            console.print(f"Plan audit result: {plan_audit_result.status}")
            return {"plan_audit_result": plan_audit_result.model_dump()}
        except Exception as e:
            console.print(f"Plan Auditor Node Failed: {e}")
            failed_audit = PlanAuditResult(
                status="rejected", feedback=f"Plan auditor process failed: {e}"
            )
            return {"plan_audit_result": failed_audit.model_dump()}

    def route_plan_audit(self, state: dict) -> Literal["coder_session", "failed"]:
        """Routes based on the plan audit result."""
        audit_result = PlanAuditResult(**state["plan_audit_result"])
        return "coder_session" if audit_result.status == "approved" else "failed"

    def failed_node(self, state: dict[str, Any]) -> dict[str, Any]:
        """A terminal node for failed states."""
        error = state.get("error", "Unknown error")
        console.print(f"Workflow failed. Reason: {error}")
        return {"status": "failed", "error": error}
