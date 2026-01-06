"""Pydantic models for configuration."""
from typing import Literal

from pydantic import BaseModel, Field, model_validator, ConfigDict
import numpy as np


class SystemConfig(BaseModel):
    """Configuration for the physical system."""
    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = Field(..., gt=0)

    @model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Validate the composition."""
        if not np.isclose(sum(self.composition.values()), 1.0):
            raise ValueError("Composition values must sum to 1.0")
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError("Elements and composition keys must match")
        return self


class ExplorationConfig(BaseModel):
    """Configuration for the exploration stage."""
    model_config = ConfigDict(extra="forbid")
    temperature: float = Field(..., gt=0)


class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""
    model_config = ConfigDict(extra="forbid")
    method: Literal["random"]
    fraction: float = Field(..., gt=0, le=1)


class FullConfig(BaseModel):
    """Top-level configuration model."""
    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    project_name: str
