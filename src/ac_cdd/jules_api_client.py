import json
from pathlib import Path
from typing import Any

import httpx
from rich.console import Console

from ac_cdd.agent_interface import AgentInterface
from ac_cdd.utils import logger

# Constants for cost estimation
# Placeholder rates: Input $0.10/1M, Output $0.40/1M
COST_PER_1M_INPUT_TOKENS = 0.10
COST_PER_1M_OUTPUT_TOKENS = 0.40

console = Console()

class JulesApiClient(AgentInterface):
    """
    Agent interface implementation for the Jules REST API.
    Manages session state and communicates with the Jules API.
    """
    BASE_URL = "https://jules.googleapis.com/v1alpha"
    STATE_FILE = Path(".jules/session_state.json")

    def __init__(self, api_key: str, base_url: str | None = None) -> None:
        self.api_key = api_key
        if base_url:
            self.BASE_URL = base_url

        self.client = httpx.Client(
            headers={
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            },
            timeout=60.0
        )
        self.STATE_FILE.parent.mkdir(exist_ok=True)

    def _get_session_id(self) -> str | None:
        if not self.STATE_FILE.exists():
            return None
        try:
            data = json.loads(self.STATE_FILE.read_text())
            return data.get("last_session_id")
        except (OSError, json.JSONDecodeError):
            return None

    def _save_session_id(self, session_id: str) -> None:
        self.STATE_FILE.write_text(json.dumps({"last_session_id": session_id}))

    def _calculate_cost(self, usage_metadata: dict[str, Any]) -> float:
        """
        Calculates the estimated cost based on token usage.
        """
        if not usage_metadata:
            return 0.0

        prompt_tokens = usage_metadata.get("promptTokenCount", 0)
        candidates_tokens = usage_metadata.get("candidatesTokenCount", 0)

        input_cost = (prompt_tokens / 1_000_000) * COST_PER_1M_INPUT_TOKENS
        output_cost = (candidates_tokens / 1_000_000) * COST_PER_1M_OUTPUT_TOKENS

        return input_cost + output_cost

    def _log_cost(self, usage_metadata: dict[str, Any]) -> None:
        """
        Logs the estimated cost using Rich.
        """
        cost = self._calculate_cost(usage_metadata)
        console.print(f"[yellow]💰 Est. Cost: ${cost:.6f}[/yellow]")
        p_tokens = usage_metadata.get('promptTokenCount', 0)
        c_tokens = usage_metadata.get('candidatesTokenCount', 0)
        logger.info(f"Token Usage - Prompt: {p_tokens}, Candidates: {c_tokens}")


    def start_task(self, prompt: str, **kwargs: Any) -> str:
        """
        Starts a new task by creating a new session via the API.
        """
        session_name = kwargs.get("session_name")
        payload = {"prompt": prompt}
        if session_name:
            payload["session_name"] = session_name

        try:
            response = self.client.post(
                f"{self.BASE_URL}/sessions",
                json=payload
            )
            response.raise_for_status()
            data = response.json()
            session_id = data.get("name")

            if "usageMetadata" in data:
                self._log_cost(data["usageMetadata"])

            if session_id:
                self._save_session_id(session_id)
                # If there's a text response in the creation, return it,
                # otherwise return success msg
                msg = data.get("response", f"API session '{session_id}' created successfully.")
                return str(msg)

            raise ConnectionError("Failed to create API session: No session ID in response.")
        except httpx.HTTPStatusError as e:
            raise ConnectionError(f"API Error: {e.response.status_code} - {e.response.text}") from e
        except httpx.RequestError as e:
            raise ConnectionError(f"Request Error: {e}") from e

    def send_message(self, prompt: str, **kwargs: Any) -> str:
        """
        Sends a message to the currently active session.
        """
        session_id = self._get_session_id()
        if not session_id:
            return self.start_task(prompt, **kwargs)

        try:
            response = self.client.post(
                f"{self.BASE_URL}/{session_id}:sendMessage",
                json={"prompt": prompt}
            )
            response.raise_for_status()
            data = response.json()

            if "usageMetadata" in data:
                self._log_cost(data["usageMetadata"])

            return str(data.get("response", "Message sent successfully."))
        except httpx.HTTPStatusError as e:
            # Check if session not found, maybe restart?
            if e.response.status_code == 404:
                logger.warning("Session not found, starting new session.")
                return self.start_task(prompt, **kwargs)
            raise ConnectionError(f"API Error: {e.response.status_code} - {e.response.text}") from e
        except httpx.RequestError as e:
            raise ConnectionError(f"Request Error: {e}") from e

    def get_status(self) -> str:
        """
        Returns the ID of the last active session.
        """
        session_id = self._get_session_id()
        if session_id:
            return f"Current active session ID: {session_id}"
        return "No active session found."
