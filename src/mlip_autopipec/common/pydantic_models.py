from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, field_validator


class SystemConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    elements: list[str]
    composition: dict[str, float]
    supercell_size: list[int] = Field(..., min_length=3, max_length=3)
    num_initial_structures: int = Field(..., gt=0)

    @field_validator("composition")
    def composition_sum_must_be_one(cls, v):
        if not abs(sum(v.values()) - 1.0) < 1e-6:
            raise ValueError("Composition values must sum to 1.0")
        return v

    @field_validator("composition")
    def composition_keys_must_be_in_elements(cls, v, info):
        if "elements" in info.data:
            if not set(v.keys()).issubset(set(info.data["elements"])):
                raise ValueError("Composition keys must be a subset of elements")
        return v


class GenerationConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    rattle_std_dev: float = Field(..., ge=0)
    volumetric_strain: float
    min_atomic_distance: float = Field(..., gt=0)


class MDEnsemble(str, Enum):
    nvt = "nvt"
    npt = "npt"


class ExplorationConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    md_calculator: str
    ensemble: MDEnsemble
    temperature_k: float = Field(..., gt=0)
    time_step_fs: float = Field(..., gt=0)
    num_steps: int = Field(..., gt=0)


class SamplingMethod(str, Enum):
    random = "random"


class SamplingConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    method: SamplingMethod
    num_samples: int = Field(..., gt=0)


class FullConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    generation: GenerationConfig
    exploration: ExplorationConfig
    sampling: SamplingConfig
