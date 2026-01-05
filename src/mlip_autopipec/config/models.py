"""Pydantic models for configuration."""
from typing import Dict, List, Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator


class SystemConfig(BaseModel):
    """Configuration for the physical system."""
    model_config = ConfigDict(extra="forbid")

    elements: List[str] = Field(..., min_length=1)
    composition: Dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = Field(..., gt=0)

    @model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Validate that the composition sums to 1.0 and elements match."""
        if not np.isclose(sum(self.composition.values()), 1.0):
            raise ValueError("Composition probabilities must sum to 1.0.")
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError("Elements and composition keys must match.")
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
    """Root configuration model."""
    model_config = ConfigDict(extra="forbid")

    project_name: str
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
