from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    model_validator,
)


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    lattice: Literal["fcc", "bcc", "hcp"]
    num_structures: int = Field(..., gt=0)

    @model_validator(mode="after")
    def validate_system(self) -> "SystemConfig":
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError("Elements and composition keys do not match.")
        if not abs(sum(self.composition.values()) - 1.0) < 1e-6:
            raise ValueError("Composition values must sum to 1.0.")
        return self


class ExplorationConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    temperature: float = Field(..., gt=0)


class SamplingConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    method: Literal["random"]
    fraction: float = Field(..., gt=0, le=1)


class FullConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
    project_name: str
