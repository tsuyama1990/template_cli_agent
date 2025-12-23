import asyncio
import httpx
import json
from typing import Any, Dict, Optional
from pathlib import Path

from ac_cdd_core.config import settings
from ac_cdd_core.utils import logger
from ac_cdd_core.services.git_ops import GitManager
from rich.console import Console

class JulesSessionError(Exception):
    pass

class JulesTimeoutError(JulesSessionError):
    pass

class JulesClient:
    """
    Client for interacting with the Google Cloud Code Agents API (Jules API).
    Refactored for the Alpha API specification (Global Endpoint).
    """

    def __init__(self) -> None:
        self.base_url = "https://jules.googleapis.com/v1alpha"
        self.timeout = settings.jules.timeout_seconds
        self.poll_interval = settings.jules.polling_interval_seconds
        self.console = Console()
        self.git = GitManager()

        # API Key is mandatory for the Alpha API
        self.api_key = settings.JULES_API_KEY
        if not self.api_key:
            logger.warning("JULES_API_KEY is not set. JulesClient will likely fail.")

    def _get_headers(self) -> Dict[str, str]:
        headers = {
            "Content-Type": "application/json",
        }
        if self.api_key:
            headers["x-goog-api-key"] = self.api_key
        return headers

    async def run_session(
        self,
        session_id: str,
        prompt: str,
        files: list[str],
        completion_signal_file: Path, # Kept for signature compatibility
        runner: Any = None, # Deprecated
    ) -> Dict[str, Any]:
        """
        Orchestrates the Jules session:
        1. Creates a session with 'AUTO_CREATE_PR' mode.
        2. Polls for completion and handles user interaction.
        3. Returns the PR URL.
        """
        if not self.api_key:
             raise JulesSessionError("Missing JULES_API_KEY configuration.")

        # 1. Prepare Source Context
        try:
            repo_url = await self.git.get_remote_url()
            # Parse owner/repo from URL (e.g., https://github.com/owner/repo.git or git@github.com:owner/repo.git)
            if "github.com" in repo_url:
                parts = repo_url.replace(".git", "").split("/")
                repo_name = parts[-1]
                owner = parts[-2].split(":")[-1] # Handle git@github.com:owner
            else:
                raise JulesSessionError(f"Unsupported repository URL format: {repo_url}. Only GitHub is supported.")

            branch = await self.git.get_current_branch()
        except Exception as e:
            raise JulesSessionError(f"Failed to determine git context: {e}") from e

        # 2. Create Session
        logger.info(f"Creating Jules Session {session_id} on branch {branch}...")

        # Global endpoint for sessions
        url = f"{self.base_url}/sessions"

        # Add file list to prompt if provided, to give context
        full_prompt = prompt
        if files:
            file_list_str = "\n".join(files)
            full_prompt += f"\n\nPlease focus on the following files:\n{file_list_str}"

        source_resource_name = f"sources/github/{owner}/{repo_name}"

        payload = {
            "prompt": full_prompt,
            "sourceContext": {
                "source": source_resource_name,
                "githubRepoContext": {
                    "startingBranch": branch
                }
            },
            "automationMode": "AUTO_CREATE_PR"
        }

        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(url, json=payload, headers=self._get_headers(), timeout=30.0)
                if response.status_code != 200:
                    error_msg = response.text
                    try:
                        err_json = response.json()
                        error_msg = err_json.get("error", {}).get("message", response.text)
                    except:
                        pass
                    raise JulesSessionError(f"Failed to create session: {response.status_code} - {error_msg}")

                resp_data = response.json()
                session_name = resp_data.get("name") # e.g., sessions/12345...
                if not session_name:
                    raise JulesSessionError("API did not return a session name.")

            except httpx.RequestError as e:
                raise JulesSessionError(f"Network error creating session: {e}") from e

        # 3. Poll for Completion & Interaction
        logger.info(f"Session created: {session_name}. Waiting for PR creation...")
        return await self.wait_for_completion(session_name)

    async def wait_for_completion(self, session_name: str) -> Dict[str, Any]:
        """
        Polls the session until it succeeds (PR created) or fails.
        Handles user interaction if the agent sends a message.
        """
        url = f"{self.base_url}/{session_name}"
        activities_url = f"{self.base_url}/{session_name}/activities"

        start_time = asyncio.get_event_loop().time()
        last_activity_count = 0
        processed_activity_ids = set()

        status_context = self.console.status("[bold green]Jules is working on your PR...", spinner="dots")
        status_context.start()

        async with httpx.AsyncClient() as client:
            try:
                while True:
                    if asyncio.get_event_loop().time() - start_time > self.timeout:
                        raise JulesTimeoutError("Timed out waiting for Jules to complete.")

                    # --- 1. Check Session Status for PR ---
                    try:
                        response = await client.get(url, headers=self._get_headers(), timeout=10.0)
                        if response.status_code != 200:
                            logger.warning(f"Polling error: {response.status_code} - {response.text}")
                            await asyncio.sleep(self.poll_interval)
                            continue

                        data = response.json()

                        # Check "outputs" for Pull Request
                        outputs = data.get("outputs", [])
                        if outputs:
                            # Usually the first output contains the PR
                            for output in outputs:
                                pr_info = output.get("pullRequest")
                                if pr_info and pr_info.get("url"):
                                    status_context.stop()
                                    pr_url = pr_info.get("url")
                                    logger.info(f"Jules Task Completed. PR Created: {pr_url}")
                                    return {"pr_url": pr_url, "status": "success", "raw": data}

                        # Check overall state if needed (though PR check is primary for AUTO_CREATE_PR)
                        # The API might not return a definitive "FAILED" state easily in Alpha,
                        # but we can check if it says "SUCCEEDED" without PR (which is an issue we handle in graph.py)
                        state = data.get("state") # e.g. SUCCEEDED, FAILED
                        if state == "FAILED":
                             status_context.stop()
                             error_msg = data.get("error", {}).get("message", "Unknown error")
                             raise JulesSessionError(f"Jules Session Failed: {error_msg}")

                        if state == "SUCCEEDED":
                             # If succeeded but no PR found in outputs yet, give it one more quick check or return
                             # We return raw data so graph.py can decide if missing PR is fatal
                             status_context.stop()
                             return {"status": "success", "raw": data} # pr_url might be None

                    except httpx.RequestError as e:
                        logger.warning(f"Network error polling status: {e}")

                    # --- 2. Check Activities for Interaction ---
                    try:
                        act_resp = await client.get(activities_url, headers=self._get_headers(), timeout=10.0)
                        if act_resp.status_code == 200:
                            act_data = act_resp.json()
                            activities = act_data.get("activities", [])

                            # API usually returns activities in some order (chronological or reverse).
                            # We just want to find *new* ones.
                            # Assuming standard list behavior, but let's be robust using IDs if available,
                            # or just count since we only need to show new stuff.

                            if len(activities) > last_activity_count:
                                # New activities found
                                # Identify the new ones. Assuming activities are appended (latest at end or start? Quickstart doesn't say).
                                # We'll just print the ones we haven't 'processed' if we can track IDs.
                                # If no IDs, we rely on count/index.

                                # Let's assume list is chronological for now, or just look at the 'diff'
                                # If we can't be sure, just look at the last one.

                                # Using the count difference to get the 'new' items
                                num_new = len(activities) - last_activity_count
                                # The API might return newest first or last.
                                # Without docs, let's assume standard Google API: often newest first in 'list' endpoints?
                                # Quickstart curl example doesn't show order.
                                # Let's try to handle both or just print the latest 'unseen'.

                                # Safe bet: Just print the activity that was added.
                                # Let's iterate all and print if not seen (by ID or signature)

                                for activity in activities:
                                    # Try to find a unique ID
                                    act_id = activity.get("name") or str(activity) # 'name' is usually the resource ID

                                    if act_id not in processed_activity_ids:
                                        processed_activity_ids.add(act_id)

                                        # It's new.
                                        # Check if it's from the AGENT (heuristic)
                                        # Or just print everything nicely

                                        # Extract message content
                                        # Structure varies: might be 'message' -> 'content' or 'reasoning' etc.
                                        # Quickstart example: response empty, "agent will send its response in the next activity"

                                        # We'll inspect 'message' or 'debugMessage'
                                        content = activity.get("message", {}).get("content") or \
                                                  activity.get("debugMessage") or \
                                                  json.dumps(activity) # Fallback

                                        if content:
                                             status_context.console.print(f"\n[dim]Activity:[/dim] {content}")

                                             # --- INTERACTION HEURISTIC ---
                                             # If it looks like a question from the agent
                                             content_clean = str(content).strip()
                                             if content_clean.endswith("?") or \
                                                "please provide" in content_clean.lower() or \
                                                "confirm" in content_clean.lower():

                                                 # Pause spinner to ask for input
                                                 status_context.stop()

                                                 self.console.print(f"\n[bold magenta]Jules Request:[/bold magenta] {content_clean}")
                                                 user_input = await asyncio.to_thread(self.console.input, "[bold green]Your Reply > [/bold green]")

                                                 # Send Reply
                                                 reply_url = f"{self.base_url}/{session_name}:sendMessage"
                                                 reply_payload = {"prompt": user_input}
                                                 await client.post(reply_url, json=reply_payload, headers=self._get_headers())

                                                 self.console.print("[dim]Reply sent.[/dim]")
                                                 status_context.start()

                                last_activity_count = len(activities)

                    except httpx.RequestError as e:
                        # Ignore activity fetch errors to keep polling status
                        pass

                    await asyncio.sleep(self.poll_interval)

            finally:
                status_context.stop()
