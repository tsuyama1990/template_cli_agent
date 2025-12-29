# ruff: noqa: D101, D102, D103, D104, D105, D107
import subprocess
from pathlib import Path


class ProcessRunner:
    """A wrapper for running external processes."""

    def run(
        self, command: list[str], cwd: Path
    ) -> subprocess.CompletedProcess[str]:
        """Runs a command in a specified directory.

        Args:
            command: The command to run as a list of strings.
            cwd: The working directory to run the command in.

        Returns:
            The CompletedProcess object from subprocess.run.
        """
        return subprocess.run(  # nosec
            command,
            capture_output=True,
            text=True,
            cwd=cwd,
        )
