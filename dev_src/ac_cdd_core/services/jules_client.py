# dev_src/ac_cdd_core/services/jules_client.py
import asyncio
import json
import logging
import time
from datetime import UTC, datetime, timedelta
from typing import Any

import google.auth
import google.auth.transport.requests
from rich.console import Console
from rich.panel import Panel

from ac_cdd_core.config import settings
from ac_cdd_core.interfaces import IJulesClient

# Eagerly get logger to avoid issues with multiprocessing
logger = logging.getLogger("AC-CDD")


class JulesTimeoutError(Exception):
    pass


console = Console()


class JulesClient(IJulesClient):
    """
    Client for interacting with the Jules API.
    Handles session creation, messaging, and polling for results.
    """

    def __init__(self, sandbox: Any) -> None:
        self.sandbox = sandbox

    def _get_headers(self) -> dict[str, str]:
        if settings.jules.api_key:
            return {
                "Content-Type": "application/json",
                "x-goog-api-key": settings.jules.api_key.get_secret_value(),
            }

        # Fallback to Application Default Credentials
        creds, project_id = google.auth.default()
        auth_req = google.auth.transport.requests.Request()
        creds.refresh(auth_req)
        return {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {creds.token}",
        }

    def _get_api_base_url(self) -> str:
        return f"https://{settings.jules.api_endpoint}"

    def _get_session_url(self, session_id: str) -> str:
        base_url = self._get_api_base_url()
        return f"{base_url}/v1alpha/sessions/{session_id}"

    async def _send_message(self, url: str, message: str) -> None:
        import httpx

        headers = self._get_headers()
        payload = {"message": {"text": message}}

        async with httpx.AsyncClient() as client:
            try:
                # Add :sendMessage to the session URL
                # Ensure no double slashes
                send_url = url.rstrip("/") + ":sendMessage"
                logger.debug(f"Sending message to {send_url}...")
                response = await client.post(send_url, headers=headers, json=payload)
                response.raise_for_status()

            except httpx.HTTPStatusError as e:
                logger.error(f"Error sending message to Jules: {e.response.text}")
                raise

    async def send_message_to_session(self, session_id: str, message: str) -> None:
        """Sends a message to an existing Jules session."""
        session_url = self._get_session_url(session_id)
        await self._send_message(session_url, message)

    async def create_session(self) -> str:
        import httpx

        headers = self._get_headers()
        base_url = self._get_api_base_url()
        session_url = f"{base_url}/v1alpha/sessions"

        # Dynamically get the remote URL from git
        from ac_cdd_core.services.git_ops import GitManager

        git_manager = GitManager()
        remote_url = await git_manager.get_remote_url()
        if not remote_url:
            raise ValueError("Could not determine git remote URL for source context.")

        payload = {
            "session": {
                "automationMode": "AUTO_CREATE_PR",
                "sourceContext": {"gitSource": {"remoteUri": remote_url}},
            }
        }

        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(session_url, headers=headers, json=payload)
                response.raise_for_status()
                session_data = response.json()
                session_name = session_data.get("name", "").split("/")[-1]
                logger.info(f"Jules session created: {session_name}")
                return session_name
            except httpx.HTTPStatusError as e:
                logger.error(f"Error creating Jules session: {e.response.text}")
                raise

    async def _get_session_activities(self, session_id: str) -> list[dict[str, Any]]:
        import httpx

        headers = self._get_headers()
        session_url = self._get_session_url(session_id)
        activities_url = f"{session_url}/activities?page_size=1000"

        async with httpx.AsyncClient() as client:
            try:
                response = await client.get(activities_url, headers=headers)
                response.raise_for_status()
                return response.json().get("activities", [])
            except httpx.HTTPStatusError as e:
                # It's possible for a session to be created but not yet have activities
                if e.response.status_code == 404:
                    return []
                logger.error(f"Error fetching session activities: {e.response.text}")
                raise

    async def wait_for_activity_type(
        self, session_id: str, activity_type: str, timeout: int
    ) -> dict[str, Any] | None:
        start_time = time.time()
        existing_activity_ids = set()

        while time.time() - start_time < timeout:
            activities = await self._get_session_activities(session_id)

            # Iterate in reverse to get the latest first
            for activity in reversed(activities):
                activity_id = activity.get("name")
                if activity_id and activity_id not in existing_activity_ids:
                    # Log the new activity for debugging
                    # Check for code_review_request and extract state
                    if activity.get("codeReviewRequest"):
                        review_state = activity["codeReviewRequest"].get("state")
                        logger.debug(f"Found new activity: {activity_id} ({review_state})")
                    else:
                        activity_summary = {
                            k: v for k, v in activity.items() if k != "codeReviewRequest"
                        }
                        logger.debug(f"Found new activity: {activity_id} {activity_summary}")

                    existing_activity_ids.add(activity_id)

                    # Check for the target activity type
                    if activity.get(activity_type):
                        return activity[activity_type]

            await asyncio.sleep(settings.jules.polling_interval)

        raise JulesTimeoutError(f"Timed out waiting for activity: {activity_type}")

    async def wait_for_completion(self, session_id: str) -> dict[str, Any]:
        """Polls a Jules session until it completes or fails."""
        start_time = datetime.now(UTC)
        timeout = timedelta(seconds=settings.jules.session_timeout)

        # Container for the most recent PR URL found
        latest_pr_url = None

        console.print(f"Waiting for Jules session '{session_id}' to complete...")

        processed_activity_ids = set()
        last_human_inquiry_ts = None

        while datetime.now(UTC) - start_time < timeout:
            try:
                activities = await self._get_session_activities(session_id)

                for activity in reversed(activities):
                    activity_id = activity.get("name")
                    if not activity_id or activity_id in processed_activity_ids:
                        continue

                    processed_activity_ids.add(activity_id)

                    # 1. Handle Human Inquiry
                    if "humanInquiry" in activity:
                        inquiry_time = activity.get("createTime")
                        if not last_human_inquiry_ts or inquiry_time > last_human_inquiry_ts:
                            last_human_inquiry_ts = inquiry_time
                            inquiry_text = activity["humanInquiry"].get("text", "No text provided.")
                            console.print(
                                Panel(
                                    f"[bold yellow]Jules needs input:[/bold yellow]\n{inquiry_text}",
                                    title="Human-in-the-Loop Inquiry",
                                    border_style="yellow",
                                )
                            )
                            user_response = await asyncio.to_thread(input, "> ")
                            await self.send_message_to_session(session_id, user_response)
                            # After sending a response, restart the loop to get fresh activities
                            continue

                    # 2. Check for PR and update latest_pr_url
                    if "codeReviewRequest" in activity:
                        pr_url = activity["codeReviewRequest"].get("reviewUrl")
                        if pr_url:
                            latest_pr_url = pr_url
                            console.print(f"PR URL found/updated: {latest_pr_url}")

                    # 3. Check for terminal states
                    if "sessionCompletion" in activity:
                        completion_data = activity["sessionCompletion"]
                        final_state = completion_data.get("finalState")

                        # If we have seen a PR, prioritize returning that even on failure
                        if latest_pr_url:
                            console.print(
                                f"Session completed with state '{final_state}'. "
                                f"Returning latest PR URL: {latest_pr_url}"
                            )
                            return {"status": "success", "pr_url": latest_pr_url}

                        # If no PR, handle based on final state
                        if final_state == "SUCCEEDED":
                            # This case is unexpected if no PR was created, but handle it.
                            console.print(
                                "[green]Session completed successfully, but no PR URL was found.[/green]"
                            )
                            return {"status": "success", "pr_url": None}

                        console.print(
                            f"[red]Session failed with final state: {final_state}[/red]"
                        )
                        return {
                            "status": "failed",
                            "error": f"Session completed with state: {final_state}",
                        }

                await asyncio.sleep(settings.jules.polling_interval)

            except Exception as e:
                logger.warning(f"Polling loop unexpected error: {e}")
                await asyncio.sleep(settings.jules.polling_interval * 2)  # Backoff

        raise JulesTimeoutError("Timed out waiting for Jules to complete.")

    async def run_session(
        self,
        session_id: str,
        prompt: str,
        target_files: list[str],
        context_files: list[str],
        require_plan_approval: bool = True,
    ) -> dict[str, Any]:
        """
        Main entry point to create and run a Jules session.
        - Creates a new session ID if one isn't provided.
        - Uploads files to the sandbox.
        - Sends the initial prompt.
        - Waits for a plan, optionally requiring user approval.
        - Waits for the session to complete and returns the PR URL.
        """
        try:
            # Step 1: Create a sandbox run
            console.print("Starting sandbox...")
            await self.sandbox.start()

            # Step 2: Upload all files (targets and context)
            all_files = list(set(target_files + context_files))
            console.print(f"Uploading {len(all_files)} files to the sandbox...")
            await self.sandbox.upload_files(all_files)

            # Step 3: Construct the prompt with file paths
            # Use sandbox-relative paths for the prompt
            sandbox_target_paths = [self.sandbox.to_sandbox_path(p) for p in target_files]
            sandbox_context_paths = [self.sandbox.to_sandbox_path(p) for p in context_files]

            prompt_parts = [
                prompt,
                "\n\n---",
                "You have been provided with the following files:",
            ]
            if sandbox_target_paths:
                prompt_parts.append("\n**Editable Target Files:**")
                prompt_parts.extend(f"- `{p}`" for p in sandbox_target_paths)
            if sandbox_context_paths:
                prompt_parts.append("\n**Read-only Context Files:**")
                prompt_parts.extend(f"- `{p}`" for p in sandbox_context_paths)

            full_prompt = "\n".join(prompt_parts)

            # Step 4: Create Jules session and send prompt
            console.print("Creating Jules session...")
            jules_session_id = await self.create_session()

            console.print(f"Sending prompt to Jules session: {jules_session_id}")
            await self.send_message_to_session(jules_session_id, full_prompt)

            # Step 5 & 6: Wait for plan and optionally approve
            if require_plan_approval:
                console.print("Waiting for plan from Jules...")
                plan_activity = await self.wait_for_activity_type(
                    jules_session_id, "solutionPlan", settings.jules.plan_timeout
                )

                if not plan_activity:
                    return {"status": "failed", "error": "Did not receive a plan from Jules."}

                console.print(
                    Panel(
                        plan_activity.get("text", "No plan text provided."),
                        title="Jules' Proposed Plan",
                        border_style="blue",
                    )
                )

                approval = await asyncio.to_thread(
                    input, "Do you approve this plan? (y/n) "
                )
                if approval.lower() != "y":
                    return {"status": "failed", "error": "Plan was not approved by the user."}

                console.print("Plan approved. Continuing session...")
                await self.send_message_to_session(jules_session_id, "That plan is approved. Please proceed.")

            # Step 7: Wait for the session to complete (PR creation)
            result = await self.wait_for_completion(jules_session_id)
            result["session_name"] = jules_session_id  # Tack on the session name for persistence
            return result

        finally:
            # Ensure sandbox is always stopped
            console.print("Stopping sandbox...")
            if self.sandbox:
                await self.sandbox.close()
