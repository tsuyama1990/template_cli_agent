"""Session management utilities for AC-CDD."""

import json
import subprocess
from datetime import datetime
from pathlib import Path

from ac_cdd_core.utils import logger


class SessionValidationError(Exception):
    """Raised when session validation fails."""

    pass


class SessionManager:
    """Manages session persistence for AC-CDD development sessions."""

    SESSION_FILE = Path(".ac_cdd_session.json")

    @classmethod
    def save_session(cls, session_id: str, integration_branch: str) -> None:
        """Save session information to file."""
        session_data = {
            "session_id": session_id,
            "integration_branch": integration_branch,
            "created_at": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
        }
        cls.SESSION_FILE.write_text(json.dumps(session_data, indent=2))
        logger.info(f"Session saved: {session_id}")

    @classmethod
    def load_session(cls) -> dict | None:
        """Load session information from file."""
        if cls.SESSION_FILE.exists():
            try:
                data = json.loads(cls.SESSION_FILE.read_text())
                logger.info(f"Session loaded: {data.get('session_id')}")
                return data
            except Exception as e:
                logger.warning(f"Failed to load session file: {e}")
                return None
        return None

    @classmethod
    def clear_session(cls) -> None:
        """Clear session file."""
        if cls.SESSION_FILE.exists():
            cls.SESSION_FILE.unlink()
            logger.info("Session file cleared")

    @classmethod
    def validate_session(cls, session_id: str, integration_branch: str) -> tuple[bool, str]:
        """
        Validate that session state is consistent with Git state.

        Returns:
            (is_valid, error_message)
        """
        # Check if integration branch exists locally
        result = subprocess.run(
            ["git", "rev-parse", "--verify", integration_branch],
            capture_output=True,
            check=False,
        )

        if result.returncode != 0:
            from ac_cdd_core.messages import RecoveryMessages

            error_msg = f"Session validation failed: {RecoveryMessages.branch_not_found(integration_branch, str(cls.SESSION_FILE))}"
            return False, error_msg

        # Check if branch exists on remote
        result = subprocess.run(
            ["git", "ls-remote", "--heads", "origin", integration_branch],
            capture_output=True,
            check=False,
        )

        if result.returncode != 0 or not result.stdout.strip():
            from ac_cdd_core.messages import RecoveryMessages

            error_msg = f"Session validation failed: {RecoveryMessages.branch_not_found(integration_branch, str(cls.SESSION_FILE))}"
            logger.warning(error_msg)
            return False, error_msg

        return True, None

    @classmethod
    def reconcile_session(cls) -> dict | None:
        """
        Attempt to reconcile session state from Git branches if session file is missing.

        Returns:
            Reconciled session data or None if reconciliation fails
        """
        logger.info("Attempting to reconcile session from Git state...")

        # List all dev/* branches
        result = subprocess.run(
            ["git", "branch", "--list", "dev/session-*"],
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.warning("Failed to list Git branches for reconciliation")
            return None

        branches = [b.strip().lstrip("* ") for b in result.stdout.splitlines() if b.strip()]

        if not branches:
            logger.info("No session branches found for reconciliation")
            return None

        # Use the most recent branch (lexicographically last, which works for our timestamp format)
        latest_branch = sorted(branches)[-1]
        session_id = latest_branch.replace("dev/", "")

        logger.info(f"Reconciled session from branch: {latest_branch}")

        session_data = {
            "session_id": session_id,
            "integration_branch": latest_branch,
            "created_at": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "reconciled": True,
        }

        # Save the reconciled session
        cls.save_session(session_id, latest_branch)

        return session_data

    @classmethod
    def get_integration_branch(cls, session_id: str) -> str:
        """Get integration branch name for session.

        Args:
            session_id: Session identifier

        Returns:
            Integration branch name (e.g., 'dev/session-20251230-120000')
        """
        from ac_cdd_core.config import settings

        return f"{settings.session.integration_branch_prefix}/{session_id}"

    @classmethod
    def load_or_reconcile_session(
        cls,
        session_id: str | None = None,
        auto_reconcile: bool = True,
    ) -> dict[str, str]:
        """Load session from parameter, file, config, or Git reconciliation.

        Args:
            session_id: Optional explicit session ID
            auto_reconcile: If True, attempt Git reconciliation if session not found

        Returns:
            Dict with 'session_id' and 'integration_branch'

        Raises:
            SessionValidationError: If no session can be found
        """
        from ac_cdd_core.config import settings

        # 1. Use explicit session ID if provided
        if session_id:
            integration_branch = cls.get_integration_branch(session_id)
            return {"session_id": session_id, "integration_branch": integration_branch}

        # 2. Try loading from session file
        saved_session = cls.load_session()
        if saved_session:
            return {
                "session_id": saved_session["session_id"],
                "integration_branch": saved_session["integration_branch"],
            }

        # 3. Try config
        if settings.session.session_id:
            session_id_from_config = settings.session.session_id
            integration_branch = cls.get_integration_branch(session_id_from_config)
            return {
                "session_id": session_id_from_config,
                "integration_branch": integration_branch,
            }

        # 4. Try Git reconciliation
        if auto_reconcile:
            reconciled = cls.reconcile_session()
            if reconciled:
                return {
                    "session_id": reconciled["session_id"],
                    "integration_branch": reconciled["integration_branch"],
                }

        # No session found
        raise SessionValidationError(
            "No session found.\n\n"
            "Recovery options:\n"
            "1. Start a new session:\n"
            "   uv run manage.py gen-cycles\n"
            "2. If you have an existing session, specify it:\n"
            "   uv run manage.py run-cycle --session <session-id>"
        )
