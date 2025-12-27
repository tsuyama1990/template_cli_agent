from pathlib import Path

import yaml
from pydantic import ValidationError

from mlip_autopipec.domain_models import MainConfig


def load_config(config_path: Path) -> MainConfig:
    """
    Loads, validates, and expands the main configuration file.
    """
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file not found at: {config_path}")

    with open(config_path) as f:
        raw_config = yaml.safe_load(f)

    try:
        config = MainConfig(**raw_config)
    except ValidationError as e:
        raise ValueError(f"Configuration validation failed: {e}") from e

    return config
