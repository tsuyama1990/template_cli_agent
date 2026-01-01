import subprocess

from mlip_autopipec.interfaces import IProcessRunner


class SubprocessRunner(IProcessRunner):
    """
    A concrete implementation of IProcessRunner that uses the subprocess module.

    This class is responsible for executing external commands and handling the
    redirection of their output.
    """

    def run(self, command: list[str], stdout_path: str) -> None:
        """
        Runs an external command using subprocess.run and redirects stdout to a file.

        Args:
            command: The command to run as a list of strings.
            stdout_path: The path to the file to redirect stdout to.

        Raises:
            subprocess.CalledProcessError: If the command returns a non-zero exit code.
        """
        with open(stdout_path, "w") as f_out:
            subprocess.run(  # noqa: S603
                command, stdout=f_out, shell=False, check=True
            )
