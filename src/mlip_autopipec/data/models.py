from typing import Any, Literal

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    FilePath,
    field_validator,
)


class DFTCompute(BaseModel):
    """Configuration for the DFT compute engine."""
    model_config = ConfigDict(extra="forbid")

    code: Literal["quantum_espresso"]
    command: str
    pseudopotentials: dict[str, str]
    ecutwfc: float
    ecutrho: float
    kpoints_density: float

    @field_validator('ecutrho')
    def ecutrho_ge_ecutwfc(cls, v, values):
        if 'ecutwfc' in values.data and v < values.data['ecutwfc']:
            raise ValueError('ecutrho must be greater than or equal to ecutwfc')
        return v


class MLIPTraining(BaseModel):
    """Configuration for the MLIP training engine."""
    model_config = ConfigDict(extra="forbid")

    model_type: Literal["ace"]
    r_cut: float
    delta_learning: bool
    base_potential: str | None = None
    loss_weights: dict[str, float]

    @field_validator('base_potential')
    def base_potential_required_for_delta_learning(cls, v, values):
        if 'delta_learning' in values.data and values.data['delta_learning'] and v is None:
            raise ValueError('base_potential must be specified for delta_learning')
        return v


class Cycle01Config(BaseModel):
    """Top-level configuration for Cycle 01."""
    model_config = ConfigDict(extra="forbid")

    dft_compute: DFTCompute
    mlip_training: MLIPTraining
    database_path: FilePath


class DFTResults(BaseModel):
    """Data model for storing results from a DFT calculation."""
    model_config = ConfigDict(arbitrary_types_allowed=True)

    energy: float
    forces: np.ndarray
    stress: np.ndarray

    @field_validator('forces', 'stress', mode='before')
    def convert_list_to_numpy_array(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        return v
