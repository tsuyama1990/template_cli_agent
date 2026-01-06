import unittest
from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

from ac_cdd_core.services.jules_client import JulesClient


class TestJulesClientLogic(unittest.IsolatedAsyncioTestCase):
    def setUp(self) -> None:
        # Patch dependencies to avoid real API calls or Auth
        self.auth_patcher = patch("google.auth.default", return_value=(MagicMock(), "test-project"))
        self.auth_patcher.start()

        # Initialize client
        with patch.object(JulesClient, "__init__", lambda x: None):  # Skip init
            self.client = JulesClient()
            self.client.base_url = "https://mock.api"
            self.client.timeout = 5
            self.client.poll_interval = 0.1
            self.client.console = MagicMock()
            self.client.manager_agent = AsyncMock()
            self.client.manager_agent.run.return_value = MagicMock(output="Manager Reply")
            self.client.credentials = MagicMock()
            self.client._get_headers = MagicMock(return_value={})
            self.client.credentials.token = "mock_token"  # noqa: S105
        self.client.git = AsyncMock()

        # Add api_client mock which is now used by wait_for_completion
        self.client.api_client = MagicMock()
        self.client.api_client.api_key = "mock_key"

    def tearDown(self) -> None:
        self.auth_patcher.stop()

    @patch("asyncio.sleep", return_value=None)
    @patch("httpx.AsyncClient")
    async def test_prioritize_inquiry_over_completed_state(
        self, mock_httpx_cls: Any, _mock_sleep: Any
    ) -> None:
        """
        Verify that if state is COMPLETED but there is a NEW inquiry,
        we prioritize answering the inquiry over returning "Success/No PR".
        """
        mock_client = AsyncMock()
        mock_httpx_cls.return_value.__aenter__.return_value = mock_client

        session_id = "sessions/123"
        activity_id = "sessions/123/activities/456"

        self.client.list_activities = MagicMock(return_value=[])
        self.client._send_message = AsyncMock()

        # --- Mock HTTP Responses ---
        r_session_completed = MagicMock()
        r_session_completed.status_code = 200
        r_session_completed.json.return_value = {"state": "COMPLETED", "outputs": []}

        r_acts_question = MagicMock()
        r_acts_question.status_code = 200
        r_acts_question.json.return_value = {
            "activities": [{"name": activity_id, "agentMessaged": {"agentMessage": "Question?"}}]
        }

        r_session_success = MagicMock()
        r_session_success.status_code = 200
        r_session_success.json.return_value = {
            "state": "SUCCEEDED",
            "outputs": [{"pullRequest": {"url": "http://github.com/pr/1"}}],
        }

        r_acts_empty = MagicMock()
        r_acts_empty.status_code = 200
        r_acts_empty.json.return_value = {"activities": []}

        # --- Execution Trace (Number of GET calls) ---
        # Loop 1:
        # 1. GET session -> COMPLETED
        # 2. GET activities (inquiry check) -> Question found, inquiry_processed=True
        # 3. GET activities (success check fallback) -> No PR found, loop continues
        # 4. GET activities (log count)
        # Loop 2:
        # 5. GET session -> SUCCEEDED with PR in outputs. Returns immediately.
        mock_client.get.side_effect = [
            r_session_completed,  # 1
            r_acts_question,  # 2
            r_acts_empty,  # 3
            r_acts_empty,  # 4
            r_session_success,  # 5
        ]

        result = await self.client.wait_for_completion(session_id)

        self.client._send_message.assert_called_once()
        assert "pr_url" in result, "Result dictionary should contain 'pr_url'"
        assert result["pr_url"] == "http://github.com/pr/1"

    @patch("asyncio.sleep", return_value=None)
    @patch("httpx.AsyncClient")
    async def test_deduplication_of_existing_activities(
        self, mock_httpx_cls: Any, _mock_sleep: Any
    ) -> None:
        """
        Verify that existing activities are IGNORED and do not trigger a reply.
        The loop should terminate early on COMPLETED state because no new inquiry is processed.
        """
        mock_client = AsyncMock()
        mock_httpx_cls.return_value.__aenter__.return_value = mock_client

        session_id = "sessions/123"
        old_activity_id = "sessions/123/activities/old"

        self.client.list_activities = MagicMock(
            return_value=[
                {"name": old_activity_id, "agentMessaged": {"agentMessage": "Old Question"}}
            ]
        )
        self.client._send_message = AsyncMock()

        # --- Mock HTTP Responses ---
        r_session_completed = MagicMock()
        r_session_completed.status_code = 200
        r_session_completed.json.return_value = {"state": "COMPLETED", "outputs": []}

        r_acts_old = MagicMock()
        r_acts_old.status_code = 200
        r_acts_old.json.return_value = {
            "activities": [
                {"name": old_activity_id, "agentMessaged": {"agentMessage": "Old Question"}}
            ]
        }

        r_acts_empty = MagicMock()
        r_acts_empty.status_code = 200
        r_acts_empty.json.return_value = {"activities": []}

        # --- Execution Trace (Number of GET calls) ---
        # Loop 1:
        # 1. GET session -> COMPLETED
        # 2. GET activities (inquiry check) -> Old activity found, inquiry_processed=False
        # 3. GET activities (success check fallback) -> No PR found.
        #    -> Loop terminates because (state==COMPLETED and not inquiry_processed) is true.
        mock_client.get.side_effect = [
            r_session_completed,  # 1
            r_acts_old,  # 2
            r_acts_empty,  # 3
        ]

        result = await self.client.wait_for_completion(session_id)

        self.client._send_message.assert_not_called()
        assert "pr_url" not in result
        assert result["status"] == "success"


if __name__ == "__main__":
    unittest.main()
