from typing import Literal

from pydantic import BaseModel, ConfigDict, FilePath, field_validator


class DFTCompute(BaseModel):
    """Configuration for DFT computations."""

    model_config = ConfigDict(extra="forbid")

    code: Literal["quantum_espresso"]
    command: str
    pseudopotentials: str
    ecutwfc: float
    ecutrho: float
    kpoints_density: float

    @field_validator("ecutrho")
    def ecutrho_ge_ecutwfc(cls, v, values):
        if "ecutwfc" in values.data and v < values.data["ecutwfc"]:
            raise ValueError("ecutrho must be greater than or equal to ecutwfc")
        return v


class MLIPTraining(BaseModel):
    """Configuration for MLIP training."""

    model_config = ConfigDict(extra="forbid")

    model_type: Literal["ace"]
    r_cut: float
    delta_learning: bool
    base_potential: str | None = None
    loss_weights: dict[str, float]

    @field_validator("base_potential")
    def base_potential_required_for_delta_learning(cls, v, values):
        if (
            "delta_learning" in values.data
            and values.data["delta_learning"]
            and v is None
        ):
            raise ValueError("base_potential must be specified for delta_learning")
        return v


class Cycle01Config(BaseModel):
    """Top-level configuration for Cycle 01."""

    model_config = ConfigDict(extra="forbid")

    dft_compute: DFTCompute
    mlip_training: MLIPTraining
    database_path: FilePath
