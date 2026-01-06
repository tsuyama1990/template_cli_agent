# -*- coding: utf-8 -*-
"""Pydantic models for configuration, defining the schema."""
from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator


class SystemConfig(BaseModel):
    """Configuration for the physical system."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = Field(..., gt=0)

    @model_validator(mode="after")
    def validate_system(self) -> "SystemConfig":
        """Validate the system configuration."""
        composition_sum = sum(self.composition.values())
        if not np.isclose(composition_sum, 1.0):
            msg = f"Composition must sum to 1.0, but sums to {composition_sum}"
            raise ValueError(msg)
        if set(self.elements) != set(self.composition.keys()):
            msg = "Elements and composition keys must match."
            raise ValueError(msg)
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
    """Top-level Pydantic model for the entire configuration."""

    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    project_name: str = Field(..., min_length=1)
