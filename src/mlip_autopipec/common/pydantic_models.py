from pydantic import BaseModel, ConfigDict, Field


class SystemConfig(BaseModel):
    """Configuration for the physical system to be generated."""

    model_config = ConfigDict(extra="forbid")

    elements: list[str] = Field(..., min_length=1)
    composition: dict[str, float]
    supercell_size: list[int] = Field(..., min_length=3, max_length=3)


class MDConfig(BaseModel):
    """Configuration for the Molecular Dynamics exploration."""

    model_config = ConfigDict(extra="forbid")

    temperature_k: float = Field(..., gt=0)
    pressure_gpa: float


class FullConfig(BaseModel):
    """The full, top-level configuration for the MLIP-AutoPipe workflow."""

    model_config = ConfigDict(extra="forbid")

    system: SystemConfig
    exploration: MDConfig
