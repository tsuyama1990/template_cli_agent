from pydantic import BaseModel, ConfigDict
from typing import List, Optional

class DFTResult(BaseModel):
    """
    A Pydantic model to store the results of a DFT calculation.
    It enforces a strict schema for the data returned by the LabellingEngine.
    """
    model_config = ConfigDict(extra='forbid')

    total_energy_ev: Optional[float] = None
    forces: Optional[List[List[float]]] = None
    stress: Optional[List[List[float]]] = None
    was_successful: bool
    error_message: Optional[str] = None

class TrainingConfig(BaseModel):
    """
    A Pydantic model for holding the configuration for the TrainingEngine.
    This ensures that the training setup is always valid and complete.
    """
    model_config = ConfigDict(extra='forbid')

    model_type: str
    learning_rate: float
    num_epochs: int
    r_cut: float
    delta_learn: bool
    baseline_potential: str
