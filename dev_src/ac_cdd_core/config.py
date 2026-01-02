import os
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

from dotenv import load_dotenv
from pydantic import Field
from pydantic_settings import BaseSettings

# Load environment variables from .env file
load_dotenv()


class Settings(BaseSettings):
    """
    Application settings.
    """

    # Constants for file names
    filename_spec: str = "ALL_SPEC.md"
    filename_arch: str = "SYSTEM_ARCHITECTURE.md"
    max_audit_retries: int = 2

    # Path resolution
    internal_template_path: Path = Field(
        default=Path("/opt/ac-cdd/templates"), description="Path to system prompts inside container"
    )

    model_config = {"extra": "allow"}

    def get_template(self, name: str) -> Path:
        """
        Resolve a template path.
        Tries to find the template in internal_template_path.
        Fallbacks to local dev path if not found (for development).
        """
        # 1. Check configured internal path (Env var or default)
        path = self.internal_template_path / name
        if path.exists():
            return path

        # 2. Fallback: Check env var AC_CDD_TEMPLATE_PATH explicitly if set
        env_path_str = os.environ.get("AC_CDD_TEMPLATE_PATH")
        if env_path_str:
            env_path = Path(env_path_str) / name
            if env_path.exists():
                return env_path

        # 3. Fallback: Local dev structure (relative to this file)
        # dev_src/ac_cdd_core/config.py -> .../dev_documents/system_prompts
        local_dev_path = (
            Path(__file__).parent.parent.parent / "dev_documents" / "system_prompts" / name
        )
        if local_dev_path.exists():
            return local_dev_path

        return path  # Return the default path even if it doesn't exist, to let caller handle error


def load_settings(config_path: str | None = None) -> Settings:
    """
    Load settings from ac_cdd_config.py.

    Tries to load configuration in the following order:
    1. Path specified in config_path argument
    2. Path specified in AC_CDD_CONFIG_PATH environment variable
    3. ac_cdd_config.py in the current working directory (Container /app)
    """

    # Determine config file path
    if config_path:
        file_path = Path(config_path)
    elif os.environ.get("AC_CDD_CONFIG_PATH"):
        file_path = Path(os.environ["AC_CDD_CONFIG_PATH"])
    else:
        # Always try current directory (container /app)
        file_path = Path.cwd() / "ac_cdd_config.py"

    if not file_path or not file_path.exists():
        # If no user config, return default Settings
        # This allows the tool to run even if the user hasn't initialized config yet
        # or if we are just running help commands
        return Settings()

    # Import the module
    spec = spec_from_file_location("ac_cdd_config", file_path)
    if not spec or not spec.loader:
        raise ImportError(f"Could not load configuration from {file_path}")

    module = module_from_spec(spec)
    sys.modules["ac_cdd_config"] = module
    spec.loader.exec_module(module)

    if not hasattr(module, "settings"):
        # If user config doesn't have settings, maybe it's just a definition file?
        # Or we should fallback to default.
        # raising error enforces correct config structure.
        raise AttributeError(f"Configuration file {file_path} must define a 'settings' object.")

    return module.settings


# Global settings object
# Access this via: from dev_src.ac_cdd_core.config import settings
try:
    settings = load_settings()
except Exception:
    # If loading fails (e.g. during installation), we might still need to import this module
    # Print warning but don't crash unless settings are actually accessed
    # print(f"Warning: Failed to load settings: {e}")
    settings = None
