import numpy as np
from pydantic import BaseModel, ConfigDict, field_validator, field_serializer
from typing import Any, Dict, Tuple, List

class DFTInputConfig(BaseModel):
    """Configuration for a Quantum Espresso calculation."""
    model_config = ConfigDict(extra="forbid")
    pseudopotentials: Dict[str, str]
    kpoints: Tuple[int, int, int]
    ecutwfc: int
    control: Dict[str, Any]

class DFTResult(BaseModel):
    """Stores the output of a QE calculation."""
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)
    energy: float
    forces: np.ndarray
    stress: np.ndarray

    @field_validator('forces', 'stress', mode='before')
    def _parse_np_array(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        if isinstance(v, np.ndarray):
            return v
        raise ValueError('Value must be a list or a numpy array')

    @field_serializer('forces', 'stress')
    def _serialize_np_array(self, v: np.ndarray) -> List[float]:
        return v.tolist()

class MLIPTrainingConfig(BaseModel):
    """Configuration for the MLIP training job."""
    model_config = ConfigDict(extra="forbid")
    model_type: str = "ACE"
    r_cut: float
    loss_weights: Dict[str, float]
