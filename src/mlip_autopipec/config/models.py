"""Pydantic models for the MLIP-AutoPipe configuration.

This module defines the data structures for the user-provided YAML configuration
file. Pydantic is used to enforce a strict schema, ensuring that all
configuration parameters are of the correct type and within valid ranges
before any computation begins. This schema-first approach is a cornerstone of
the application's robustness.
"""

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator


class SystemConfig(BaseModel):
    """Defines the physical material system to be generated.

    This model captures all parameters related to the atomic structure,
    including the elements, composition, and crystal lattice.
    """

    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1, description="A list of chemical symbols.")
    composition: dict[str, float] = Field(
        ..., description="The fractional composition of each element."
    )
    lattice: Literal["fcc", "bcc", "hcp"] = Field(..., description="The crystal lattice type.")
    num_structures: int = Field(
        ..., gt=0, description="The number of initial seed structures to generate."
    )

    @model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Ensure the composition is consistent and sums to 1.0."""
        if not np.isclose(sum(self.composition.values()), 1.0):
            raise ValueError("Composition values must sum to 1.0")
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError("Elements and composition keys must match")
        return self


class ExplorationConfig(BaseModel):
    """Defines the parameters for the (placeholder) exploration stage.

    In Cycle 1, this is a simple placeholder, but the schema is defined for
    forward compatibility.
    """

    model_config = ConfigDict(extra="forbid")
    temperature: float = Field(..., gt=0, description="The simulation temperature in Kelvin.")


class SamplingConfig(BaseModel):
    """Defines the parameters for the sampling stage."""

    model_config = ConfigDict(extra="forbid")
    method: Literal["random"] = Field(..., description="The sampling method to use.")
    fraction: float = Field(..., gt=0, le=1, description="The fraction of structures to select.")


class FullConfig(BaseModel):
    """The root model for the entire pipeline configuration.

    This model aggregates all the other configuration components into a single,
    validated object.
    """

    model_config = ConfigDict(extra="forbid")
    system: SystemConfig = Field(..., description="The physical system configuration.")
    exploration: ExplorationConfig = Field(..., description="The exploration stage configuration.")
    sampling: SamplingConfig = Field(..., description="The sampling stage configuration.")
    project_name: str = Field(..., description="The name of the project, used for output files.")
