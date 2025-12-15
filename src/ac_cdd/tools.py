import shutil
import subprocess

from .utils import logger


class ToolNotFoundError(Exception):
    pass

class ToolWrapper:
    def __init__(self, command: str):
        self.command = command
        if not shutil.which(command):
            raise ToolNotFoundError(f"Command '{command}' not found in PATH.")

    def run(
        self,
        args: list[str],
        capture_output: bool = False,
        check: bool = True,
        text: bool = True
    ) -> subprocess.CompletedProcess:
        full_cmd = [self.command] + args
        logger.debug(f"Running command: {' '.join(full_cmd)}")
        try:
            return subprocess.run(  # noqa: S603
                full_cmd, capture_output=capture_output, check=check, text=text
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed: {e}")
            raise

class GitClient:
    """Git操作を行うクライアント"""
    def __init__(self) -> None:
        self.git = ToolWrapper("git")

    def get_diff(self, target: str = "HEAD") -> str:
        """Git差分を取得"""
        try:
            # capture_output=True to get stdout
            result = self.git.run(["diff", target], capture_output=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            # git diffのエラー（差分なし等）は空文字として扱う場合
            return ""
