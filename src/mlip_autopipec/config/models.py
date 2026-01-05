"""Pydantic models for configuration management."""

import math
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SystemConfig(BaseModel):
    """Configuration for the physical system."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = Field(..., gt=0)

    @model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Validate that composition sums to 1.0 and elements match."""
        composition_sum_error = "Composition fractions must sum to 1.0."
        if not math.isclose(sum(self.composition.values()), 1.0):
            raise ValueError(composition_sum_error)

        element_mismatch_error = "Elements in 'elements' list and 'composition' dict must match."
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError(element_mismatch_error)
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
