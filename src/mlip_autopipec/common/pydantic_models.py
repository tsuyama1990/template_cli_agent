from pydantic import BaseModel, ConfigDict, Field


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    supercell_size: tuple[int, int, int]


class MDConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    temperature_k: float = Field(..., gt=0)
    pressure_gpa: float = Field(..., ge=0)
    time_step_fs: float = Field(..., gt=0)
    total_steps: int = Field(..., gt=0)


class SamplingConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    method: str = "random"
    number_of_samples: int = Field(..., gt=0)


class FullConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    exploration: MDConfig
    sampling: SamplingConfig
