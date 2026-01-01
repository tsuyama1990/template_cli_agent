from typing import Any, Dict

from pydantic import BaseModel, ConfigDict


class DFTComputeConfig(BaseModel):
    """Configuration for DFT compute settings."""

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
    """Configuration for MLIP training."""

    model_config = ConfigDict(extra="forbid")

    model_type: str
    r_cut: float
    delta_learning: bool
    base_potential: str
    loss_weights: Dict[str, float]


class MainConfig(BaseModel):
    """Main configuration model."""

    model_config = ConfigDict(extra="forbid")

    system: Dict[str, Any]
    dft_compute: DFTComputeConfig
    mlip_training: MLIPTrainingConfig
