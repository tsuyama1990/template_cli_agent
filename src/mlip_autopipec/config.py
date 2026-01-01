from enum import Enum
from typing import Any

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
)


class ModelType(str, Enum):
    ACE = "ACE"


class DFTInputConfig(BaseModel):
    """Pydantic model for DFT input parameters."""

    model_config = ConfigDict(extra="forbid")

    pseudopotentials: dict[str, str] = Field(
        ..., description="Mapping of element to pseudopotential filename."
    )
    kpoints: tuple[int, int, int] = Field(..., description="K-points grid.")
    ecutwfc: int = Field(..., description="Wavefunction cutoff energy.")
    control: dict[str, Any] = Field(
        ..., description="Control parameters for the DFT calculation."
    )


class DFTResult(BaseModel):
    """Pydantic model for storing DFT calculation results."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    energy: float = Field(..., description="Total energy from the calculation.")
    forces: np.ndarray = Field(..., description="Forces on each atom.")
    stress: np.ndarray = Field(..., description="Stress tensor.")

    @field_validator("forces", "stress", mode="before")
    @classmethod
    def validate_numpy_array(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        if not isinstance(v, np.ndarray):
            raise ValueError("Value must be a numpy array or a list.")
        return v

    @field_serializer("forces", "stress")
    def serialize_numpy_array(self, v: np.ndarray) -> list:
        return v.tolist()


class MLIPTrainingConfig(BaseModel):
    """Pydantic model for MLIP training configuration."""

    model_config = ConfigDict(extra="forbid")

    model_type: ModelType = Field(
        ..., description="Type of the MLIP model to be trained."
    )
    r_cut: float = Field(..., description="Cutoff radius for the potential.")
    loss_weights: dict[str, float] = Field(
        ..., description="Weights for the loss function components."
    )


class Settings(BaseModel):
    """Manages application-wide settings."""

    qe_command: str = Field(
        "pw.x", description="The command to execute Quantum Espresso."
    )
