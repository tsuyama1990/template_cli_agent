"""Defines the protocol for running external processes."""

import subprocess
from typing import Protocol


class ProcessRunner(Protocol):
    """A protocol for running external processes."""

    def run(self, command: list[str], input_str: str) -> str:
        """
        Runs a command and returns its standard output.

        Args:
            command: The command to run as a list of strings.
            input_str: A string to pass to the process's standard input.

        Returns:
            The standard output of the process as a string.
        """
        ...


class SubprocessRunner:
    """A concrete implementation of ProcessRunner using subprocess."""

    def run(self, command: list[str], input_str: str) -> str:
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
