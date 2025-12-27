"""Pydantic models for data structures."""
from pydantic import BaseModel, ConfigDict
from typing import List, Optional


class DFTResult(BaseModel):
    """Represents the results of a DFT calculation."""
    model_config = ConfigDict(extra='forbid')

    total_energy_ev: Optional[float] = None
    forces: Optional[List[List[float]]] = None
    stress: Optional[List[List[float]]] = None
    was_successful: bool
    error_message: Optional[str] = None


class TrainingConfig(BaseModel):
    """Represents the configuration for an MLIP training run."""
    model_config = ConfigDict(extra='forbid')

    model_type: str
    learning_rate: float
    num_epochs: int
    r_cut: float
    delta_learn: bool
    baseline_potential: str
