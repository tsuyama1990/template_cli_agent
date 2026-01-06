# src/mlip_autopipec/common/pydantic_models.py
"""Core Pydantic models for configuration and data structures."""

from pydantic import BaseModel, ConfigDict, Field


class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(
        ..., min_length=1, description="List of element symbols (e.g., ['Fe', 'Pt'])."
    )
    composition: dict[str, float] = Field(
        ..., description="Composition of the elements (e.g., {'Fe': 0.75, 'Pt': 0.25})."
    )
    supercell_size: list[int] = Field(
        ..., min_length=3, max_length=3, description="Dimensions of the supercell [x, y, z]."
    )


class MDConfig(BaseModel):
    """Configuration for the Molecular Dynamics exploration."""

    model_config = ConfigDict(extra="forbid")

    temperature_k: float = Field(
        ..., gt=0, description="Temperature for the MD simulation in Kelvin."
    )
    pressure_gpa: float | None = Field(
        None, description="Pressure for the MD simulation in GPa. If None, NVT ensemble is used."
    )
    n_steps: int = Field(..., gt=0, description="Number of MD steps to perform.")


class SamplingConfig(BaseModel):
    """Configuration for the sampling process."""

    model_config = ConfigDict(extra="forbid")

    method: str = Field("Random", description="Sampling method to use ('Random' or 'FPS').")
    n_samples: int = Field(
        ..., gt=0, description="Number of structures to sample from the trajectory."
    )


class FullConfig(BaseModel):
    """The full, top-level configuration for the MLIP-AutoPipe workflow."""

    model_config = ConfigDict(extra="forbid")

    system: SystemConfig
    exploration: MDConfig
    sampling: SamplingConfig
    db_path: str = Field("mlip_autopipec_data.db", description="Path to the output ASE database.")
