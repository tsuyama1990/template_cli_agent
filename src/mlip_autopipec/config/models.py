"""Pydantic models for configuration, defining the schema for the project."""
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(
        ..., min_length=1, description="List of chemical symbols."
    )
    composition: dict[str, float] = Field(
        ..., description="Mapping from element to its fractional composition."
    )
    lattice: Literal["fcc", "bcc", "hcp"] = Field(
        ..., description="Crystal lattice type."
    )
    num_structures: int = Field(
        ..., gt=0, description="Number of initial seed structures to generate."
    )

    @model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Validate that the composition is consistent with the elements list."""
        import math

        composition_sum = sum(self.composition.values())
        if not math.isclose(composition_sum, 1.0):
            raise ValueError(
                f"The sum of compositions must be 1.0, but it is {composition_sum}"
            )

        if set(self.composition.keys()) != set(self.elements):
            raise ValueError(
                "The elements in the composition dictionary must exactly match the elements list."
            )

        return self


class ExplorationConfig(BaseModel):
    """Configuration for the exploration stage."""

    model_config = ConfigDict(extra="forbid")

    temperature: float = Field(..., gt=0, description="Simulation temperature in Kelvin.")


class SamplingConfig(BaseModel):
    """Configuration for the sampling stage."""

    model_config = ConfigDict(extra="forbid")

    method: Literal["random"] = Field(..., description="Sampling method.")
    fraction: float = Field(
        ..., gt=0, le=1, description="Fraction of structures to select."
    )


class FullConfig(BaseModel):
    """Top-level configuration model for the entire pipeline."""

    model_config = ConfigDict(extra="forbid")

    project_name: str = Field(
        ..., min_length=1, description="Name of the project or run."
    )
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
