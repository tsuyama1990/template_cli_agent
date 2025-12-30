
class RecoveryMessages:
    """Centralized messages for error recovery instructions."""

    @staticmethod
    def branch_not_found(branch_name: str, session_file: str) -> str:
        return (
            f"Branch '{branch_name}' not found locally.\n"
            "This usually happens if the session file is stale or the branch was deleted manually.\n\n"
            "To fix this:\n"
            f"1. Delete the session file: rm {session_file}\n"
            "2. Re-run your command (it will try to reconcile or start fresh)."
        )

    @staticmethod
    def remote_branch_missing(branch_name: str) -> str:
        return (
            f"Warning: Branch '{branch_name}' does not exist on origin.\n"
            "You might be working in a detached state or offline.\n"
            "Ensure you push your branch: git push -u origin {branch_name}"
        )
