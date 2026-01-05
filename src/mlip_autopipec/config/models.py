"""Pydantic models for configuration, defining the schema."""
from typing import Literal

import pydantic


class SystemConfig(pydantic.BaseModel):
    """Configuration for the physical system."""

    elements: list[str]
    composition: dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = pydantic.Field(..., gt=0)

    @pydantic.model_validator(mode="after")
    def validate_composition(self) -> "SystemConfig":
        """Validate the composition."""
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError("Elements and composition keys do not match.")
        if not pydantic.VERSION.startswith("1."):
            # Pydantic v2
            import math
            if not math.isclose(sum(self.composition.values()), 1.0):
                raise ValueError("Composition must sum to 1.0.")
        else:
            # Pydantic v1
            if abs(sum(self.composition.values()) - 1.0) > 1e-9:
                raise ValueError("Composition must sum to 1.0.")
        return self


class ExplorationConfig(pydantic.BaseModel):
    """Configuration for the exploration stage."""

    temperature: float = pydantic.Field(..., gt=0)


class SamplingConfig(pydantic.BaseModel):
    """Configuration for the sampling stage."""

    method: Literal["random"]
    fraction: float = pydantic.Field(..., gt=0, le=1)


class FullConfig(pydantic.BaseModel):
    """Top-level configuration model."""

    model_config = pydantic.ConfigDict(extra="forbid")

    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    project_name: str
