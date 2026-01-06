# src/mlip_autopipec/common/pydantic_models.py
"""
Pydantic models for configuration validation.
Following the schema-first development approach, these models define the expected
structure and types for the user's configuration file.
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


class SystemConfig(BaseModel):
    """Configuration for the physical system."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    supercell_size: list[int] = Field(..., min_length=3, max_length=3)

    @field_validator("composition")
    def check_composition_sum(cls, v: dict[str, float]) -> dict[str, float]:  # noqa: N805
        """Validate that the composition fractions sum to 1.0."""
        if not abs(sum(v.values()) - 1.0) < 1e-6:
            error_msg = "Composition fractions must sum to 1.0"
            raise ValueError(error_msg)
        return v


class MDConfig(BaseModel):
    """Configuration for the Molecular Dynamics exploration."""

    model_config = ConfigDict(extra="forbid")

    temperature_k: float = Field(..., gt=0)
    pressure_gpa: float = Field(..., ge=0)
    timestep_fs: float = Field(0.5, gt=0)
    n_steps: int = Field(1000, gt=0)


class SamplingConfig(BaseModel):
    """Configuration for the sampling process."""

    model_config = ConfigDict(extra="forbid")

    method: Literal["random", "fps"] = "random"
    n_samples: int = Field(100, gt=0)


class FullConfig(BaseModel):
    """The root configuration model."""

    model_config = ConfigDict(extra="forbid")

    system: SystemConfig
    exploration: MDConfig
    sampling: SamplingConfig
    db_path: str = "mlip_data.db"
