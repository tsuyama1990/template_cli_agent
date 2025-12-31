"""Centralized configuration for the MLIP AutoPipe."""

from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field

class Settings(BaseSettings):
    """Defines the application's configuration settings."""
    model_config = SettingsConfigDict(
        env_file=".env",
        env_nested_delimiter='__',
    )

    db_path: str = Field("mlip.db", description="Path to the ASE database.")
    dft_command: str = Field("pw.x", description="The command to execute Quantum Espresso.")
    pseudo_dir: str = Field(".", description="Directory for pseudopotentials.")
    ecutwfc: float = Field(50.0, description="Plane-wave cutoff energy.")
    model_path: str = Field("model.pt", description="Path to save the trained model.")
    cutoff: float = Field(5.0, description="Cutoff radius for the MACE model.")
