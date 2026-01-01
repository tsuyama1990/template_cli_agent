import os

import yaml
from pydantic import ValidationError

from mlip_autopipec.data.models import Cycle01Config


def load_config(config_path: str) -> Cycle01Config:
    """
    Loads, validates, and returns the configuration from a YAML file.

    Args:
        config_path: The path to the YAML configuration file.

    Returns:
        A Cycle01Config object.

    Raises:
        FileNotFoundError: If the config file is not found.
        yaml.YAMLError: If the config file is not valid YAML.
        ValidationError: If the configuration fails Pydantic validation.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at '{config_path}'")

    with open(config_path) as f:
        config_dict = yaml.safe_load(f)

    # Resolve the database path relative to the config file location
    config_dir = os.path.dirname(os.path.abspath(config_path))
    db_path = os.path.join(config_dir, config_dict.get('database_path', ''))
    config_dict['database_path'] = db_path

    try:
        return Cycle01Config(**config_dict)
    except ValidationError as e:
        # Re-raise with a more informative message if needed, or handle here
        raise e
