from ac_cdd_core.services.jules_client import JulesClient
from ac_cdd_core.services.session_manager import SessionManager


class CoderService:
    def __init__(self, jules_client: JulesClient, session_manager: SessionManager):
        self.jules_client = jules_client
        self.session_manager = session_manager

    async def run_coder_session(self, cycle_id: str, resume_mode: bool) -> dict:
        cycle = await self.session_manager.get_cycle(cycle_id)

        if resume_mode and cycle.jules_session_id:
            result = await self.jules_client.wait_for_completion(cycle.jules_session_id)
        else:
            result = await self.jules_client.run_session()
            await self.session_manager.update_cycle_state(
                cycle_id, jules_session_id=result["session_name"], status="in_progress"
            )

        if result.get("status") == "success" and result.get("pr_url"):
            return {"status": "ready_for_audit"}
        else:
            return {"status": "failed"}
