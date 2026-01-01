# src/mlip_autopipec/configs/models.py

from typing import Any

from pydantic import BaseModel, ConfigDict


class DFTComputeConfig(BaseModel):
    """Configuration for the DFT compute settings."""

    model_config = ConfigDict(extra="forbid")

    code: str
    command: str
    pseudopotentials: str
    ecutwfc: float
    ecutrho: float
    kpoints_density: float
    smearing: str
    degauss: float


class MLIPTrainingConfig(BaseModel):
    """Configuration for the MLIP training process."""

    model_config = ConfigDict(extra="forbid")

    model_type: str
    r_cut: float
    delta_learning: bool
    base_potential: str
    loss_weights: dict[str, float]


class MainConfig(BaseModel):
    """Top-level configuration model."""

    model_config = ConfigDict(extra="forbid")

    system: dict[str, Any]
    dft_compute: DFTComputeConfig
    mlip_training: MLIPTrainingConfig
