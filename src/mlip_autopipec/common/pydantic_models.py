"""
Core Pydantic models for configuration and data structures.

This module defines the schema for the entire configuration of the MLIP-AutoPipe
application. It uses Pydantic for data validation and settings management.
"""


from pydantic import BaseModel, ConfigDict, Field


class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(..., description="List of element symbols (e.g., ['Fe', 'Pt']).")
    composition: dict[str, float] = Field(
        ..., description="Dictionary of element compositions (e.g., {'Fe': 0.75, 'Pt': 0.25})."
    )
    supercell_size: list[int] = Field(
        ...,
        min_length=3,
        max_length=3,
        description="Dimensions of the supercell (e.g., [3, 3, 3]).",
    )


class MDConfig(BaseModel):
    """Configuration for the Molecular Dynamics exploration stage."""

    model_config = ConfigDict(extra="forbid")

    temperature_k: float = Field(..., gt=0, description="MD simulation temperature in Kelvin.")
    pressure_gpa: float = Field(..., description="MD simulation pressure in GPa.")
    # Add other MD parameters as needed


class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""

    model_config = ConfigDict(extra="forbid")

    method: str = Field("random", description="Sampling method to use ('random' or 'fps').")
    n_samples: int = Field(..., gt=0, description="Number of samples to select.")
    # Add other sampling parameters as needed


class FullConfig(BaseModel):
    """Top-level configuration model for the entire pipeline."""

    model_config = ConfigDict(extra="forbid")

    system: SystemConfig
    exploration: MDConfig
    sampling: SamplingConfig
    # Add other top-level config sections as needed
