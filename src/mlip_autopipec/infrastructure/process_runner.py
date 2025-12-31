"""Concrete implementation of the ProcessRunnerPort using subprocess."""

from typing import List
import subprocess
from mlip_autopipec.domain.ports import ProcessRunnerPort


class SubprocessRunner(ProcessRunnerPort):
    """A concrete implementation of ProcessRunner using subprocess."""

    def run(self, command: List[str], input_str: str) -> str:
        """
        Runs a command using subprocess.run and returns its stdout.
        """
        result = subprocess.run(  # noqa: S603
            command,
            input=input_str,
            capture_output=True,
            text=True,
            check=False,  # We will check for errors by parsing the output
        )
        return result.stdout
