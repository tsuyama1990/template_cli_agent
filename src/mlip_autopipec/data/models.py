from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    field_validator,
    model_validator,
)


class DFTCompute(BaseModel):
    """Configuration for the DFT compute engine."""

    model_config = ConfigDict(extra="forbid")

    code: Literal["quantum_espresso"]
    command: str
    pseudopotentials: Dict[str, str]
    ecutwfc: float
    ecutrho: float
    kpoints_density: float

    @model_validator(mode="after")
    def ecutrho_ge_ecutwfc(self):
        if self.ecutrho < self.ecutwfc:
            raise ValueError("ecutrho must be greater than or equal to ecutwfc")
        return self


class MLIPTraining(BaseModel):
    """Configuration for the MLIP training engine."""

    model_config = ConfigDict(extra="forbid")

    model_type: Literal["mace", "ace"]
    r_cut: float
    delta_learning: bool
    base_potential: Optional[str] = None
    loss_weights: Dict[str, float]

    @model_validator(mode="after")
    def base_potential_required_for_delta_learning(self):
        if self.delta_learning and self.base_potential is None:
            raise ValueError("base_potential must be specified for delta_learning")
        return self


class MinimalSystem(BaseModel):
    """Represents the user's minimal input in `input.yaml`."""

    model_config = ConfigDict(extra="forbid")

    elements: List[str]
    composition: Union[str, Dict[str, int]]
    simulation_temperature: Optional[List[float]] = None


class System(BaseModel):
    """Represents the fully defined system after heuristic expansion."""

    model_config = ConfigDict(extra="forbid")
    elements: List[str]
    composition: Dict[str, int]
    structure_type: Literal["alloy", "ionic", "covalent", "molecular"]


class StructureGeneration(BaseModel):
    """Contains settings for Module A: Structure Generator."""

    model_config = ConfigDict(extra="forbid")

    generation_strategy: Literal["sqs", "random"]
    supercell_size: Union[int, List[int]]
    strains: List[float]


class FullConfig(BaseModel):
    """Top-level configuration, fully expanded by the heuristic engine."""

    model_config = ConfigDict(extra="forbid")

    system: System
    generation: StructureGeneration
    dft_compute: DFTCompute
    mlip_training: MLIPTraining


class DFTResults(BaseModel):
    """Data model for storing results from a DFT calculation."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    energy: float
    forces: np.ndarray
    stress: np.ndarray

    @field_validator("forces", "stress", mode="before")
    def convert_list_to_numpy_array(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        return v
