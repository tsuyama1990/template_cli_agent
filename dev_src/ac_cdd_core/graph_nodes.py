import asyncio
from typing import Any

from rich.console import Console

from .config import settings
from .sandbox import SandboxRunner
from .services.audit_orchestrator import AuditOrchestrator
from .services.jules_client import JulesClient
from .services.llm_reviewer import LLMReviewer
from .state import CycleState

console = Console()


class CycleNodes:
    """
    Encapsulates the logic for each node in the AC-CDD workflow graph.
    Decoupled from the graph topology definition in GraphBuilder.
    """

    def __init__(self, sandbox_runner: SandboxRunner, jules_client: JulesClient):
        self.sandbox = sandbox_runner
        self.jules = jules_client
        self.audit_orchestrator = AuditOrchestrator(jules_client, sandbox_runner)
        self.llm_reviewer = LLMReviewer(sandbox_runner)

    async def architect_session_node(self, state: CycleState) -> dict[str, Any]:
        """Node for Architect Agent (Jules)."""
        console.print("[bold blue]Starting Architect Session...[/bold blue]")

        # Start session with Jules
        # In architect mode, we expect Jules to generate SYSTEM_ARCHITECTURE.md, SPEC.md, UAT.md

        # Context is ALL_SPEC.md (loaded from settings)
        context_files = [settings.filename_spec, "ARCHITECT_INSTRUCTION.md"]

        # Resolve path for instructions to ensure they are found.
        # JulesClient will need to load it.

        result = await self.jules.start_architect_session(
            files=context_files, instruction_template="ARCHITECT_INSTRUCTION.md"
        )

        if result.get("status") == "success":
            return {"status": "architect_completed"}
        else:
            return {"status": "architect_failed", "error": result.get("error")}

    async def coder_session_node(self, state: CycleState) -> dict[str, Any]:
        """Node for Coder Agent (Jules or Aider)."""
        cycle_id = state["cycle_id"]
        iteration = state["iteration_count"]

        console.print(
            f"[bold green]Starting Coder Session for Cycle {cycle_id} "
            f"(Iteration {iteration})...[/bold green]"
        )

        # Determine if we use Jules (Iter 0) or Aider (Iter > 0)
        # For now, let's assume Jules for Iter 0 and Aider for fixing

        if iteration == 0:
            # Initial implementation by Jules
            # Context: SYSTEM_ARCHITECTURE.md, SPEC.md (for this cycle)
            # spec_file = f"CYCLE{cycle_id}/SPEC.md"
            # files = [settings.filename_arch, spec_file]

            # Start/Resume session
            # We need to wait for a PLAN from Jules (Architecture -> Plan -> Code)
            # This might be handled inside jules_client or audit_orchestrator

            # Actually, for Iteration 0, we might just want to trigger Jules
            # to code based on the spec
            pass

            # Placeholder for actual logic
            await asyncio.sleep(1)
            return {"status": "ready_for_audit"}

        else:
            # Fix implementation by Aider (or Jules Fixer)
            # We have feedback from Auditor
            console.print("[yellow]Starting Fixer Agent...[/yellow]")

            # Run Aider
            # ...

            state["iteration_count"] += 1
            return {"status": "ready_for_audit"}

    async def auditor_node(self, state: CycleState) -> dict[str, Any]:
        """Node for Auditor Agent (Aider/LLM)."""
        console.print("[bold magenta]Starting Auditor...[/bold magenta]")

        result = await self.audit_orchestrator.run_audit(state)

        return {"audit_result": result, "status": result.status}

    async def uat_evaluate_node(self, state: CycleState) -> dict[str, Any]:
        """Node for UAT Evaluation."""
        console.print("[bold cyan]Running UAT Evaluation...[/bold cyan]")

        # Run UAT tests via Sandbox
        # Check UAT.md requirements

        # ...

        return {"status": "cycle_completed"}

    def check_coder_outcome(self, state: CycleState) -> str:
        status = state.get("status")
        if status == "ready_for_audit":
            return "ready_for_audit"
        elif (
            status == "architect_failed"
        ):  # Should not happen in coder graph but handling for safety
            return "failed"
        return "completed"  # Default fallback

    def check_audit_outcome(self, state: CycleState) -> str:
        audit_res = state.get("audit_result")
        if not audit_res:
            return "rejected_retry"  # Fallback

        if audit_res.status == "approved":
            return "approved"

        # Check max retries using settings
        if state["iteration_count"] >= settings.max_audit_retries:
            console.print("[bold red]Max audit retries reached. Stopping.[/bold red]")
            return "rejected_max_retries"

        return "rejected_retry"
