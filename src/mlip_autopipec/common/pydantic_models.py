from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    supercell_size: list[int] = Field(..., min_items=3, max_items=3)

class MDConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    temperature_k: float = Field(..., gt=0)
    pressure_gpa: float
    time_step_fs: float = Field(..., gt=0)
    total_steps: int = Field(..., ge=10)

class SamplingConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    method: Literal["random", "fps"] = "random"
    num_samples: int = Field(..., ge=1)

class FullConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    exploration: MDConfig
    sampling: SamplingConfig
